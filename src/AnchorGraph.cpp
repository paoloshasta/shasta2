// Shasta.
#include "AnchorGraph.hpp"
#include "Anchor.hpp"
#include "AnchorPair.hpp"
#include "approximateTopologicalSort.hpp"
#include "deduplicate.hpp"
#include "dominatorTree.hpp"
#include "findLinearChains.hpp"
#include "findReachableVertices.hpp"
#include "graphvizToHtml.hpp"
#include "Journeys.hpp"
#include "Markers.hpp"
#include "memoryInformation.hpp"
#include "orderPairs.hpp"
#include "partialTransitiveReduction.hpp"
#include "performanceLog.hpp"
#include "ReadId.hpp"
#include "timestamp.hpp"
#include "tmpDirectory.hpp"
#include "transitiveReduction.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "fstream.hpp"
#include <queue>
#include "tuple.hpp"

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<AnchorGraph>;



// Construct the AnchorGraph from the Journeys.
// Only include edges with at least the specified minCoverage.
AnchorGraph::AnchorGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    uint64_t minEdgeCoverage,
    double minEdgeCoverageFraction) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<AnchorGraph>(*this)
{
    writeMemoryStatistics("AnchorGraph construction begins");
    AnchorGraph& anchorGraph = *this;

    // Create the vertices, one for each AnchorId.
    // In the AnchorGraph, vertex_descriptors are AnchorIds.
    const uint64_t anchorCount = anchors.size();
    for(AnchorId anchorId=0; anchorId<anchorCount; anchorId++) {
        add_vertex(anchorGraph);
    }



    // Loop over possible source vertices to create edges.
    nextEdgeId = 0;
    vector<AnchorPair> anchorPairs;
    for(AnchorId anchorIdA=0; anchorIdA<anchorCount; anchorIdA++) {
        const uint64_t coverageA = anchors[anchorIdA].coverage();
        AnchorPair::createChildren(anchors, journeys, anchorIdA, 0, anchorPairs);

        for(const AnchorPair& anchorPair: anchorPairs) {
            const uint64_t edgeCoverage = anchorPair.size();

            // The edge must satisfy both minEdgeCoverage and minEdgeCoverageFraction.
            if(edgeCoverage >= minEdgeCoverage) {
                // minEdgeCoverage is satisfied.
                // We must check if minEdgeCoverageFraction is also satisfied.

                const AnchorId anchorIdB = anchorPair.anchorIdB;
                const uint64_t coverageB = anchors[anchorIdB].coverage();

                if(double(edgeCoverage) >= minEdgeCoverageFraction * double(min(coverageA, coverageB))) {
                    addEdge(anchorIdA, anchorIdB, anchorPair.orientedReadIds, true);
                }
            }
        }
    }

    cout << "The anchor graph has " << num_vertices(*this) <<
        " vertices and " << num_edges(*this) << " edges." << endl;

    writeMemoryStatistics("AnchorGraph construction ends");
}



// Constructor from binary data.
AnchorGraph::AnchorGraph(const MappedMemoryOwner& mappedMemoryOwner, const string& name) :
    MappedMemoryOwner(mappedMemoryOwner),
    MultithreadedObject<AnchorGraph>(*this)
{
    load(name);
}



void AnchorGraph::save(ostream& s) const
{
    boost::archive::binary_oarchive archive(s);
    archive << *this;
}



void AnchorGraph::load(istream& s)
{
    boost::archive::binary_iarchive archive(s);
    archive >> *this;
}



void AnchorGraph::save(const string& name) const
{
    // If not using persistent binary data, do nothing.
    if(largeDataFileNamePrefix.empty()) {
        return;
    }

    // First save to a string.
    std::ostringstream s;
    save(s);
    const string dataString = s.str();

    // Now save the string to binary data.
    MemoryMapped::Vector<char> data;
    data.createNew(largeDataName(name), largeDataPageSize);
    data.resize(dataString.size());
    const char* begin = dataString.data();
    const char* end = begin + dataString.size();
    copy(begin, end, data.begin());
}



void AnchorGraph::load(const string& name)
{
    // Access the binary data.
    MemoryMapped::Vector<char> data;
    try {
        data.accessExistingReadOnly(name);
    } catch (std::exception&) {
        throw runtime_error(name + " is not available.");
    }
    const string dataString(data.begin(), data.size());

    // Load it from here.
    std::istringstream s(dataString);
    try {
        load(s);
    } catch(std::exception& e) {
        throw runtime_error(string("Error reading " + name + ": ") + e.what());
    }
}



void AnchorGraph::transitiveReduction(
    uint64_t transitiveReductionMaxEdgeCoverage,
    uint64_t transitiveReductionMaxDistance)
{
    AnchorGraph& anchorGraph = *this;
    performanceLog << timestamp << "AnchorGraph transitive reduction begins." << endl;

    // Initially make sure all edges are flag as "useForAssembly".
    // The transitive reduction process sets useForAssembly to false
    // for edges removed by transitive reduction.
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        anchorGraph[e].useForAssembly = true;
    }

    // Loop over edge coverage.
    // At each iteration we only consider edges with this coverage.
    vector<edge_descriptor> edgesToProcess;
    vector<edge_descriptor> edgesToRemove;
    for(uint64_t edgeCoverage=1; edgeCoverage<=transitiveReductionMaxEdgeCoverage; edgeCoverage++) {

        // Gather edges with this coverage.
        edgesToProcess.clear();
        BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
            if(anchorGraph[e].coverage() == edgeCoverage) {
                edgesToProcess.push_back(e);
            }
        }

        // If there are none, there is nothing to do.
        if(edgesToProcess.empty()) {
            continue;
        }

        // Loop over all edges with this coverage.
        // This can be multithreaded.
        edgesToRemove.clear();
        for(const edge_descriptor e: edgesToProcess) {
            if(transitiveReductionCanRemove(e, transitiveReductionMaxDistance)) {
                edgesToRemove.push_back(e);
            }
        }

        // Turn off the useForAssembly flag for edges removed at this iteration over coverage.
        for(const edge_descriptor e: edgesToRemove) {
            anchorGraph[e].useForAssembly = false;
        }
        cout << "Edge coverage " << edgeCoverage <<
            ": processed " << edgesToProcess.size() <<
            " edges and flagged " << edgesToRemove.size() << endl;
    }

    uint64_t useForAssemblyCount = 0;
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        if(anchorGraph[e].useForAssembly) {
            ++useForAssemblyCount;
        }
    }
    cout << useForAssemblyCount << " anchor graph edges flagged for use in assembly out of " <<
        num_edges(anchorGraph) << " total." << endl;

    performanceLog << timestamp << "AnchorGraph transitive reduction ends." << endl;
}



bool AnchorGraph::transitiveReductionCanRemove(
    edge_descriptor e,
    uint64_t transitiveReductionMaxDistance) const
{
    const AnchorGraph& anchorGraph = *this;
    const uint64_t edgeCoverage = anchorGraph[e].coverage();

    const vertex_descriptor v0 = source(e, anchorGraph);
    const vertex_descriptor v1 = target(e, anchorGraph);

    const bool debug = ((anchorIdToString(v0) == "45549+") and (anchorIdToString(v1) == "78505-"));

    // Do a forward BFS starting at v0, using edges
    // still marked as "use for assembly"
    // with coverage greater than edgeCoverage
    // and with maximum distance (number of edges)
    // equal to transitiveReductionMaxDistance.
    // If we encounter v1, return true.
    std::queue<vertex_descriptor> q;
    q.push(v0);

    // A map to store vertices already encountered and their distance from v0.
    std::map<vertex_descriptor, uint64_t> m;
    m.insert(make_pair(v0, 0));



    // Main BFS loop.
    while(not q.empty()) {

        // Dequeue a vertex.
        const vertex_descriptor vA = q.front();
        q.pop();
        const auto itA = m.find(vA);
        SHASTA2_ASSERT(itA != m.end());
        const uint64_t distanceA = itA->second;
        const uint64_t distanceB = distanceA + 1;

        // Loop over its out-edges still marked as useForAssembly
        // and with sufficient coverage.
        BGL_FORALL_OUTEDGES(vA, eAB, anchorGraph, AnchorGraph) {
            const AnchorGraphEdge& edgeAB = anchorGraph[eAB];
            if(not edgeAB.useForAssembly) {
                continue;
            }

            // Only use edges with higher coverage for the BFS,
            if(edgeAB.coverage() <= edgeCoverage) {
                continue;
            }

            // If we reached v1, return true;
            const vertex_descriptor vB = target(eAB, anchorGraph);
            if(vB == v1) {
                if(debug) {
                    cout << "Edge " << anchorIdToString(v0) << " " << anchorIdToString(v1) <<
                        " flagged by transitive reduction." << endl;
                }
                return true;
            }

            // If we already encountered vB, don't do anything.
            if(m.contains(vB)) {
                continue;
            }

            if(distanceB < transitiveReductionMaxDistance) {
                q.push(vB);
                m.insert(make_pair(vB, distanceB));
            }
        }
    }

    // If getting here we did not encounter v1 in the BFS loop.
    if(debug) {
        cout << "Edge " << anchorIdToString(v0) << " " << anchorIdToString(v1) <<
            " not flagged by transitive reduction." << endl;
    }
    return false;
}



// Return the reverse complement of a vertex.
// In the AnchorGraph, vertex_descriptors are AnchorIds.
AnchorGraph::vertex_descriptor AnchorGraph::reverseComplement(vertex_descriptor v) const
{
    const AnchorId anchorId = v;
    const AnchorId anchorIdRc = reverseComplementAnchorId(anchorId);
    const vertex_descriptor vRc = anchorIdRc;
    SHASTA2_ASSERT(vRc != v);
    return vRc;
}



// Return the reverse complement of an edge.
AnchorGraph::edge_descriptor AnchorGraph::reverseComplement(edge_descriptor e) const
{
    const AnchorGraph& anchorGraph = *this;

    const vertex_descriptor v0 = source(e, anchorGraph);
    const vertex_descriptor v1 = target(e, anchorGraph);

    const vertex_descriptor v0Rc = reverseComplement(v0);
    const vertex_descriptor v1Rc = reverseComplement(v1);

    auto[eRc, edgeExists] = boost::edge(v1Rc, v0Rc, anchorGraph);
    SHASTA2_ASSERT(edgeExists);
    SHASTA2_ASSERT(eRc != e);

    return eRc;
}
