// Shasta.
#include "AssemblyGraph.hpp"
#include "AnchorGraph.hpp"
#include "DisjointSets.hpp"
#include "findLinearChains.hpp"
#include "findReachableVertices.hpp"
#include "Journeys.hpp"
#include "LocalAssembly6.hpp"
#include "LocalAssembly7.hpp"
#include "memoryInformation.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
#include "ReadFollowing4.hpp"
#include "RestrictedAnchorGraph.hpp"
#include "SegmentStepSupport.hpp"
#include "Tangle.hpp"
#include "TangleMatrix.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "chrono.hpp"
#include "fstream.hpp"
#include <ranges>



// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<AssemblyGraph>;



// Initial construction from the AnchorGraph.
AssemblyGraph::AssemblyGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    const AnchorGraph& anchorGraph,
    const Options& options,
    shared_ptr<const HomopolymerModel> homopolymerModelPointer) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<AssemblyGraph>(*this),
    anchors(anchors),
    journeys(journeys),
    options(options),
    homopolymerModelPointer(homopolymerModelPointer),
    orderById(*this),
    edgesToBeAssembledA(orderById),
    edgesToBeAssembledB(orderById)
{
    AssemblyGraph& assemblyGraph = *this;

    // Create a filtered AnchorGraph containing only the edges marked as "useForAssembly".
    class EdgePredicate {
    public:
        bool operator()(const AnchorGraph::edge_descriptor& e) const
        {
            return (*anchorGraph)[e].useForAssembly;
        }
        EdgePredicate(const AnchorGraph& anchorGraph) : anchorGraph(&anchorGraph) {}
        EdgePredicate() : anchorGraph(0) {}
        const AnchorGraph* anchorGraph;
    };
    using FilteredAnchorGraph = boost::filtered_graph<AnchorGraph, EdgePredicate>;
    FilteredAnchorGraph filteredAnchorGraph(anchorGraph, EdgePredicate(anchorGraph));



    // Find linear chains of edges in the FilteredAnchorGraph.
    vector< std::list<AnchorGraph::edge_descriptor> > chains;
    findLinearChains(filteredAnchorGraph, 1, chains);



    // The AnchorGraph is guaranteed to be exactly strand-symmetric
    // and we want the AssemblyGraph to also be strand-symmetric.
    // For each vertex the reverse conplemented vertex will be stored in vRc.
    // For each edge the reverse complemented edge will be stored in eRc.



    // Before we look for pairs of reverse complemented chains, we need to
    // normalize any isolated circular chains so they begin and end at the lowest
    // numbered AnchorId in the chain. This happens rarely.
    for(std::list<AnchorGraph::edge_descriptor>& chain: chains)
    {
        // In a circular chain, the first and last vertex must coincide.
        const AnchorGraph::edge_descriptor e0 = chain.front();
        const AnchorGraph::edge_descriptor e1 = chain.back();
        const AnchorGraph::vertex_descriptor v0 = source(e0, anchorGraph);
        const AnchorGraph::vertex_descriptor v1 = target(e1, anchorGraph);
        if(v0 != v1) {
            // This is not a circular chain. Do nothing.
            continue;
        }

        // In a circular chain, the first vertex, which is also the last,
        // must have in-degree and out-degree 1.
        if(in_degree(v0, anchorGraph) != 1) {
            // This is not an isolated circular chain. Do nothing.
            continue;
        }
        if(out_degree(v0, anchorGraph) != 1) {
            // This is not an isolated circular chain. Do nothing.
            continue;
        }

        // This is an isolated circular chain. Gather its vertices.
        // In the AnchorGraph, vertex descriptors are AnchorIds.
        vector<AnchorGraph::vertex_descriptor> chainVertices;
        for(const AnchorGraph::edge_descriptor e: chain) {
            chainVertices.push_back(source(e, anchorGraph));
        }

        // This happens rarely and is pathological, so let's write out the chain.
        cout << "Found the following isolated circular anchor chain." << endl;
        cout << "Before normalization:" << endl;
        for(const AnchorId anchorId: chainVertices) {
            cout << anchorIdToString(anchorId) << " ";
        }
        cout << endl;

        // Find the position of the lowest numbered vertex_descriptor
        // (which is also an AnchorId).
        const auto it = std::min_element(chainVertices.begin(), chainVertices.end());

        // Rotate the chain so the lowest numbered vertex_descriptor
        // (which is also an AnchorId) is at the beginning.
        std::rotate(chainVertices.begin(), it, chainVertices.end());

        cout << "After normalization:" << endl;
        for(const AnchorId anchorId: chainVertices) {
            cout << anchorIdToString(anchorId) << " ";
        }
        cout << endl;

        // Use the rotated vertices to created the normalized (rotated) chain.
        chainVertices.push_back(chainVertices.front());
        std::list<AnchorGraph::edge_descriptor> rotatedChain;
        for(uint64_t i1=1; i1<chainVertices.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const AnchorGraph::vertex_descriptor v0 = chainVertices[i0];
            const AnchorGraph::vertex_descriptor v1 = chainVertices[i1];
            auto [e, edgeExists] = boost::edge(v0, v1, anchorGraph);
            SHASTA2_ASSERT(edgeExists);
            rotatedChain.push_back(e);
        }

        // Replace the initial chain with the rotated chain.
        chain.swap(rotatedChain);
    }



    // Index the chains by their first AnchorGraph::edge_descriptor.
    std::map<AnchorGraph::edge_descriptor, uint64_t> chainMap;
    for(uint64_t chainId=0; chainId<chains.size(); chainId++) {
        const std::list<AnchorGraph::edge_descriptor>& chain = chains[chainId];
        const AnchorGraph::edge_descriptor eFirst = chain.front();
        SHASTA2_ASSERT(not chainMap.contains(eFirst));
        chainMap.insert(make_pair(eFirst, chainId));
    }



    // Now we can find pairs of reverse complemented chains.
    vector< pair<uint64_t, uint64_t> > chainPairs;
    for(uint64_t chainId=0; chainId<chains.size(); chainId++) {
        const std::list<AnchorGraph::edge_descriptor>& chain = chains[chainId];
        const AnchorGraph::edge_descriptor e0 = chain.front();
        const AnchorGraph::edge_descriptor e1 = chain.back();
        const AnchorGraph::edge_descriptor e0Rc = anchorGraph.reverseComplement(e0);
        const AnchorGraph::edge_descriptor e1Rc = anchorGraph.reverseComplement(e1);

        // The reverse complemented chain begins at e1Rc.
        const uint64_t chainIdRc = chainMap.at(e1Rc);
        SHASTA2_ASSERT(chainIdRc != chainId);
        const std::list<AnchorGraph::edge_descriptor>& chainRc = chains[chainIdRc];
        SHASTA2_ASSERT(chainRc.front() == e1Rc);
        SHASTA2_ASSERT(chainRc.back() == e0Rc);

        // Only store each pair once.
        if(chainId < chainIdRc) {
            chainPairs.push_back(make_pair(chainId, chainIdRc));
        }
    }
    SHASTA2_ASSERT(2 * chainPairs.size() == chains.size());




    // Generate vertices.
    // At this stage there is a vertex for each AnchorGraph vertex
    // that is at the beginning or end of a linear chain,
    // so there is only one vertex for a given AnchorId.
    // However, in later stages of the AssemblyGraph there can be more than one vertex
    // with a given AnchorId. So the vertexMap is only used in this constructor.
    std::map<AnchorId, vertex_descriptor> vertexMap;
    for(const auto& chain: chains) {
        const AnchorId anchorId0 = anchorGraph[chain.front()].anchorIdA;
        const AnchorId anchorId1 = anchorGraph[chain.back()].anchorIdB;

        if(not vertexMap.contains(anchorId0)) {
            const vertex_descriptor v0 = add_vertex(AssemblyGraphVertex(anchorId0, nextVertexId++), assemblyGraph);
            vertexMap.insert(make_pair(anchorId0, v0));
        }

        if(not vertexMap.contains(anchorId1)) {
            const vertex_descriptor v1 = add_vertex(AssemblyGraphVertex(anchorId1, nextVertexId++), assemblyGraph);
            vertexMap.insert(make_pair(anchorId1, v1));
        }
    }
    SHASTA2_ASSERT(vertexMap.size() == num_vertices(assemblyGraph));



    // Store the reverse complement of each vertex.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        AssemblyGraphVertex& vertex = assemblyGraph[v];
        const AnchorId anchorId = vertex.anchorId;
        const AnchorId anchorIdRc = reverseComplementAnchorId(anchorId);
        if(anchorId < anchorIdRc) {
            const vertex_descriptor vRc = vertexMap.at(anchorIdRc);
            AssemblyGraphVertex& vertexRc = assemblyGraph[vRc];
            vertex.vRc = vRc;
            vertexRc.vRc = v;
        }
    }



    // Each pair of reverse complement chains generates a pair
    // reverse complement of AssembleGraph edges.
    for(const auto&[chainId, chainIdRc]: chainPairs) {
        const std::list<AnchorGraph::edge_descriptor>& chain = chains[chainId];
        const std::list<AnchorGraph::edge_descriptor>& chainRc = chains[chainIdRc];

        // Get the first and last AnchorId of these two chains.
        const AnchorId anchorId0 = anchorGraph[chain.front()].anchorIdA;
        const AnchorId anchorId1 = anchorGraph[chain.back()].anchorIdB;
        const AnchorId anchorId0Rc = anchorGraph[chainRc.front()].anchorIdA;
        const AnchorId anchorId1Rc = anchorGraph[chainRc.back()].anchorIdB;
        SHASTA2_ASSERT(anchorId1Rc == reverseComplementAnchorId(anchorId0));
        SHASTA2_ASSERT(anchorId0Rc == reverseComplementAnchorId(anchorId1));

        // Get the corresponding AssemblyGraph vertices.
        const vertex_descriptor v0 = vertexMap.at(anchorId0);
        const vertex_descriptor v1 = vertexMap.at(anchorId1);
        const vertex_descriptor v0Rc = vertexMap.at(anchorId0Rc);
        const vertex_descriptor v1Rc = vertexMap.at(anchorId1Rc);
        SHASTA2_ASSERT(assemblyGraph[v0].vRc == v1Rc);
        SHASTA2_ASSERT(assemblyGraph[v1].vRc == v0Rc);

        // Now we can generate the two AssemblyGraph edges
        // corresponding to these chains.
        // Each AnchorGraph edge in a chain contributes a step to this AssemblyGraph edge.

        auto[e, wasAdded] = add_edge(v0, v1, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
        SHASTA2_ASSERT(wasAdded);
        AssemblyGraphEdge& edge = assemblyGraph[e];
        for(const AnchorGraph::edge_descriptor eA: chain) {
            const AnchorPair anchorPair = anchorGraph.getAnchorPair(eA);
            edge.emplace_back(anchorPair, anchorPair.getAverageOffset(anchors));
        }

        auto[eRc, wasAddedRc] = add_edge(v0Rc, v1Rc, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
        SHASTA2_ASSERT(wasAddedRc);
        AssemblyGraphEdge& edgeRc = assemblyGraph[eRc];
        for(const AnchorGraph::edge_descriptor eA: chainRc) {
            const AnchorPair anchorPair = anchorGraph.getAnchorPair(eA);
            edgeRc.emplace_back(anchorPair, anchorPair.getAverageOffset(anchors));
        }

        // Store the eRc fields.
        edge.eRc = eRc;
        edgeRc.eRc = e;
    }


#if 0
    // OLD CODE TO GENERATE EDGES WITHOUT SETTING THE eRc FIELDS.
    // Generate the edges. There is an edge for each linear chain.
    for(const auto& chain: chains) {
        const AnchorId anchorId0 = anchorGraph[chain.front()].anchorIdA;
        const AnchorId anchorId1 = anchorGraph[chain.back()].anchorIdB;

        const vertex_descriptor v0 = vertexMap.at(anchorId0);
        const vertex_descriptor v1 = vertexMap.at(anchorId1);

        edge_descriptor e;
        bool edgeWasAdded;
        tie(e, edgeWasAdded) = add_edge(v0, v1, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
        AssemblyGraphEdge& edge = assemblyGraph[e];

        // Each AnchorGraph edge in the chain contributes a step to this AssemblyGraph edge.
        for(const AnchorGraph::edge_descriptor eA: chain) {
            const AnchorPair anchorPair = anchorGraph.getAnchorPair(eA);
            edge.emplace_back(anchorPair, anchorPair.getAverageOffset(anchors));
        }
    }
#endif

}



// Deserialize constructor.
AssemblyGraph::AssemblyGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    const Options& options,
    shared_ptr<const HomopolymerModel> homopolymerModelPointer,
    const string& stage) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<AssemblyGraph>(*this),
    anchors(anchors),
    journeys(journeys),
    options(options),
    homopolymerModelPointer(homopolymerModelPointer),
    orderById(*this),
    edgesToBeAssembledA(orderById),
    edgesToBeAssembledB(orderById)
{
    load(stage);
}



// Top level assembly function.
void AssemblyGraph::simplifyAndAssemble()
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t maxIterationCount = 6;

    writeMemoryStatistics("AssemblyGraph::simplifyAndAssemble begins");

    // Initial output.
    writeIntermediateStageIfRequested("A");

    // Iterate detangling, read following, phasing.
    for(uint64_t iteration=0; iteration<maxIterationCount; ++iteration) {
        if(not simplifyIteration(iteration)) {
            cout << "No assembly graph changes at simplify iteration " << iteration << "." << endl;
            break;
        }
    }
    writeIntermediateStageIfRequested("B");

    // Make the AssemblyGraph single-stranded.
    makeSingleStranded();
    write("Final-NoSequence");

    // Sequence assembly.
    assembleAllStrandSymmetric();
    write("Final");
    writeFasta("Final");

    writeMemoryStatistics("AssemblyGraph::simplifyAndAssemble ends");
}



// This runs one iteration of phasing, detangling, read following.
// It returns true if any changes in the AssemblyGraph were made.
bool AssemblyGraph::simplifyIteration(uint64_t iteration)
{
    const string iterationString = to_string(iteration);
    const bool debug = false;

    const uint64_t oldNextEdgeId = nextEdgeId;

    if(debug) write(iterationString + "A");

    // Remove small connected components and strict pruning.
    removeSmallComponents();
    strictPrune();
    if(debug) write(iterationString + "B");
    strandSymmetricCompress();
    if(debug) write(iterationString + "B-Compressed");

    // Remove or simplify bubbles likely caused by errors.
    bubblePairCleanup();
    if(debug) write(iterationString + "C");
    strandSymmetricCompress();
    if(debug) write(iterationString + "C-Compressed");

    // Phasing.
    strandSymmetricPhaseSuperbubbleChains();
    if(debug) write(iterationString + "D");
    strandSymmetricCompress();
    if(debug) write(iterationString + "D-Compressed");

    // Vertex detangling.
    detangleVertices(iterationString);
    if(debug) write(iterationString + "E");
    strandSymmetricCompress();
    if(debug) write(iterationString + "E-Compressed");

    // Superbubble detangling and read following.
    detangleAndReadFollowingSuperbubbles(iterationString);
    if(debug) write(iterationString + "F");
    strandSymmetricCompress();
    if(debug) write(iterationString + "F-Compressed");

    // Edge detangling.
    detangleEdges(iterationString);
    if(debug) write(iterationString + "G");
    strandSymmetricCompress();
    if(debug) write(iterationString + "G-Compressed");

    // Strand separation for self-complementary tangles.
    separateStrands1(iterationString);
    if(debug) write(iterationString + "H");
    strandSymmetricCompress();
    if(debug) write(iterationString + "H-Compressed");

    const bool changesWereMade = (nextEdgeId > oldNextEdgeId);
    return changesWereMade;
}



void AssemblyGraph::check(bool writeDetails) const
{
    cout << timestamp << "AssemblyGraph::check begins." << endl;

    const AssemblyGraph& assemblyGraph = *this;

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        if(edge.empty()) {
            continue;
        }

        if(writeDetails) {
            cout << "Checking edge " << edge.id << " with " << edge.size() << " steps." << endl;
        }

        // Check that the first/last AnchorIds of this edge are consistent
        // with the ones in the source/target vertices.
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);

        const AnchorId anchorId0 = assemblyGraph[v0].anchorId;
        const AnchorId anchorId1 = assemblyGraph[v1].anchorId;

        if(writeDetails) {
            cout << "Source AnchorId for this edge is " << anchorIdToString(anchorId0) << endl;
            cout << "Target AnchorId for this edge is " << anchorIdToString(anchorId1) << endl;
            cout << "AnchorIds for each step:" << endl;
            for(const AssemblyGraphEdgeStep& step: edge) {
                cout << anchorIdToString(step.anchorPair.anchorIdA) << " " <<
                    anchorIdToString(step.anchorPair.anchorIdB) << endl;
            }
        }

        SHASTA2_ASSERT(edge.front().anchorPair.anchorIdA == anchorId0);
        SHASTA2_ASSERT(edge.back().anchorPair.anchorIdB == anchorId1);

        // Check that AnchorPairs in this edge are adjacent to each other.
        for(uint64_t i1=1; i1<edge.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            SHASTA2_ASSERT(edge[i0].anchorPair.anchorIdB == edge[i1].anchorPair.anchorIdA);
        }
    }



    // Check consistency of the vRc fields in the vertices.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphVertex& vertex = assemblyGraph[v];
        const vertex_descriptor vRc = vertex.vRc;
        SHASTA2_ASSERT(vRc != null_vertex());
        SHASTA2_ASSERT(vRc != v);
        const AssemblyGraphVertex& vertexRc = assemblyGraph[vRc];
        SHASTA2_ASSERT(vertexRc.vRc != null_vertex());
        SHASTA2_ASSERT(vertexRc.vRc != vRc);
        if(vertexRc.vRc != v) {
            cout << "Check fails at vertices " << assemblyGraph[v].id << " " << assemblyGraph[vRc].id << endl;
            cout << "v " << v << endl;
            cout << "assemblyGraph[v].id " << assemblyGraph[v].id << endl;
            cout << "assemblyGraph[v].vRc " << assemblyGraph[v].vRc << endl;
            cout << "vRc " << vRc << endl;
            cout << "assemblyGraph[vRc].id " << assemblyGraph[vRc].id << endl;
            cout << "assemblyGraph[vRc].vRc " << assemblyGraph[vRc].vRc << endl;

        }
        SHASTA2_ASSERT(vertexRc.vRc == v);
    }



    // Check consistency of the eRc fields in the edges.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        if(writeDetails) {
            cout << "Checking " << edge.id << endl;
        }
        const edge_descriptor eRc = edge.eRc;
        if(writeDetails) {
            if(eRc == assemblyGraphNullEdge) {
                cout << "Segment " << assemblyGraph[eRc].id <<
                    " has null reverse complement." << endl;
            }
        }
        SHASTA2_ASSERT(eRc != assemblyGraphNullEdge);
        SHASTA2_ASSERT(eRc != e);
        const AssemblyGraphEdge& edgeRc = assemblyGraph[eRc];
        if(writeDetails) {
            cout << "Checking reverse complemented edge pair " <<
                edge.id << " " << edgeRc.id << endl;
        }
        SHASTA2_ASSERT(edgeRc.eRc != assemblyGraphNullEdge);
        SHASTA2_ASSERT(edgeRc.eRc != eRc);
        SHASTA2_ASSERT(edgeRc.eRc == e);

        // Also check the steps of these two edges.
        const uint64_t n = edge.size();
        SHASTA2_ASSERT(edgeRc.size() == n);
        for(uint64_t i=0; i<n; i++) {
            if(writeDetails) {
                cout << "Checking step " << i << " of " << n << endl;
            }
            const AssemblyGraphEdgeStep& step = edge[i];
            const AssemblyGraphEdgeStep& stepRc = edgeRc[n - 1 - i];
            SHASTA2_ASSERT(step.offset == stepRc.offset);
            const AnchorPair& anchorPair = step.anchorPair;
            const AnchorPair& anchorPairRc = stepRc.anchorPair;
            const vector<OrientedReadId>& orientedReadIds = anchorPair.orientedReadIds;
            const vector<OrientedReadId>& orientedReadIdsRc = anchorPairRc.orientedReadIds;
            const uint64_t coverage = orientedReadIds.size();
            SHASTA2_ASSERT(orientedReadIdsRc.size() == coverage);

            for(uint64_t j=0; j<coverage; j++) {
                const OrientedReadId orientedReadId = orientedReadIds[j];
                const OrientedReadId orientedReadIdRc = orientedReadIdsRc[j];
                SHASTA2_ASSERT(orientedReadId.getReadId() == orientedReadIdRc.getReadId());
                SHASTA2_ASSERT(orientedReadId.getStrand() == 1 - orientedReadIdRc.getStrand());
            }
        }
    }

    cout << timestamp << "AssemblyGraph::check ends." << endl;
}



uint64_t AssemblyGraphEdge::offset() const
{
    uint64_t sum = 0;
    for(const auto& step: *this) {
        sum += step.offset;
    }
    return sum;
}



void AssemblyGraph::write(const string& stage)
{
    performanceLog << timestamp << "AssemblyGraph::write begins for stage " << stage << endl;

    cout << "Stage " << stage << ": " <<
        num_vertices(*this) << " vertices, " <<
        num_edges(*this) << " segments. Next segment id is " << nextEdgeId << "." << endl;

    const uint64_t zeroLengthSegmentCount = countZeroLengthSegments();
    if(zeroLengthSegmentCount) {
        cout << "There are " << zeroLengthSegmentCount << " zero length segments." << endl;
    }

    if((options.memoryMode == "filesystem") and options.keepBinaryData) {
        save(stage);
    }
    writeGfa("Assembly-" + stage + ".gfa");
    writeGraphviz("Assembly-" + stage + ".dot");
    writeCsv("Assembly-" + stage + ".csv");
    writeSequenceLengthByCoverageCsv("Assembly-SequenceLengthByCoverage-" + stage + ".csv");

    if(options.writeAssemblyDetails) {
        writeDetailsCsv("AssemblyDetails-" + stage + ".csv");
    }

    performanceLog << timestamp << "AssemblyGraph::write ends for stage " << stage << endl;
}



void AssemblyGraph::writeIntermediateStageIfRequested(const string& name)
{
    if(options.writeIntermediateAssemblyStages) {
        write(name);
    }
}



void AssemblyGraph::writeFasta(const string& stage) const
{
    const AssemblyGraph& assemblyGraph = *this;;

    ofstream fasta("Assembly-" + stage + ".fasta");

    vector<shasta2::Base> sequence;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        edge.getSequence(sequence);

        fasta << ">" << edge.id << "\n";
        copy(sequence.begin(), sequence.end(), ostream_iterator<shasta2::Base>(fasta));
        fasta << "\n";
    }
}



void AssemblyGraph::writeGfa(const string& fileName) const
{
    ofstream gfa(fileName);
    writeGfa(gfa);
}



void AssemblyGraph::writeGfa(ostream& gfa) const
{
    const AssemblyGraph& assemblyGraph = *this;

    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";

    // Each edge generates a gfa segment.
    vector<shasta2::Base> sequence;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        const double coverage = edge.lengthWeightedAverageCoverage();

        // Record type.
        gfa << "S\t";

        // Name.
        gfa << edge.id << "\t";

        // Sequence.
        if(edge.wasAssembled) {
            edge.getSequence(sequence);
            copy(sequence.begin(), sequence.end(), ostream_iterator<shasta2::Base>(gfa));
            const uint64_t length = sequence.size();
            gfa << "\tLN:i:" << length;
            gfa << "\tRC:i:" << uint64_t(std::round(coverage * double(length)));
            gfa << "\n";

        } else {
            if(edge.empty()) {
                gfa << "*\tLN:i:0\n";

            } else {
                const uint64_t offset = edge.offset();
                gfa << "*\tLN:i:" << offset;
                gfa << "\tRC:i:" << uint64_t(std::round(coverage * double(offset)));
                gfa << "\n";
            }
        }
    }



    // For each vertex, generate a link between each pair of
    // incoming/outgoing edges.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        BGL_FORALL_INEDGES(v, e0, assemblyGraph, AssemblyGraph) {
            const uint64_t id0 = assemblyGraph[e0].id;
            BGL_FORALL_OUTEDGES(v, e1, assemblyGraph, AssemblyGraph) {
                const uint64_t id1 = assemblyGraph[e1].id;

                gfa <<
                    "L\t" <<
                    id0 << "\t+\t" <<
                    id1 << "\t+\t*\n";
            }
        }
    }



    // If two or more segments end at a vertex with out-degree 0,
    // add a dummy segment so the two segments are joined at the end.
    // Similarly in the opposite direction.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        if((out_degree(v, assemblyGraph) == 0) and (in_degree(v, assemblyGraph) > 1)) {
            const string dummyEdgeName = "Vertex" + to_string(id(v));
            gfa << dummyEdgeName << "\t*\tLN:i:0\n";
            BGL_FORALL_INEDGES(v, segment, assemblyGraph, AssemblyGraph) {
                gfa <<
                    "L\t" << id(segment) << "\t+\t" <<
                    dummyEdgeName << "\t+\t*\n";
            }
        }
        if((in_degree(v, assemblyGraph) == 0) and (out_degree(v, assemblyGraph) > 1)) {
            const string dummyEdgeName = "Vertex" + to_string(id(v));
            gfa << dummyEdgeName << "\t*\tLN:i:0\n";
            BGL_FORALL_OUTEDGES(v, segment, assemblyGraph, AssemblyGraph) {
                gfa <<
                    "L\t" << dummyEdgeName << "\t+\t" <<
                    id(segment) << "\t+\t*\n";
            }
        }
    }
}



// GFA output of a subset of Segments.
// The Segments must be sorted by id.
void AssemblyGraph::writeGfa(const string& fileName, const vector<Segment>& segments) const
{
    ofstream gfa(fileName);
    writeGfa(gfa, segments);
}



void AssemblyGraph::writeGfa(ostream& gfa, const vector<Segment>& segments) const
{
    const AssemblyGraph& assemblyGraph = *this;

    SHASTA2_ASSERT(std::is_sorted(segments.begin(), segments.end(), orderById));

    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";

    vector<shasta2::Base> sequence;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        if(not std::binary_search(segments.begin(), segments.end(), e, orderById)) {
            continue;
        }
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        if(not std::binary_search(segments.begin(), segments.end(), e)) {
            continue;
        }
        const double coverage = edge.lengthWeightedAverageCoverage();

        // Record type.
        gfa << "S\t";

        // Name.
        gfa << edge.id << "\t";

        // Sequence.
        if(edge.wasAssembled) {
            edge.getSequence(sequence);
            copy(sequence.begin(), sequence.end(), ostream_iterator<shasta2::Base>(gfa));
            const uint64_t length = sequence.size();
            gfa << "\tLN:i:" << length;
            gfa << "\tRC:i:" << uint64_t(std::round(coverage * double(length)));
            gfa << "\n";

        } else {
            if(edge.empty()) {
                gfa << "*\tLN:i:0\n";

            } else {
                const uint64_t offset = edge.offset();
                gfa << "*\tLN:i:" << offset;
                gfa << "\tRC:i:" << uint64_t(std::round(coverage * double(offset)));
                gfa << "\n";
            }
        }
    }



    // For each vertex, generate a link between each pair of
    // incoming/outgoing edges.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        BGL_FORALL_INEDGES(v, e0, assemblyGraph, AssemblyGraph) {
            if(not std::binary_search(segments.begin(), segments.end(), e0, orderById)) {
                continue;
            }
            const uint64_t id0 = assemblyGraph[e0].id;
            BGL_FORALL_OUTEDGES(v, e1, assemblyGraph, AssemblyGraph) {
                if(not std::binary_search(segments.begin(), segments.end(), e1, orderById)) {
                    continue;
                }
                const uint64_t id1 = assemblyGraph[e1].id;

                gfa <<
                    "L\t" <<
                    id0 << "\t+\t" <<
                    id1 << "\t+\t*\n";
            }
        }
    }
}



void AssemblyGraph::writeGraphviz(const string& fileName) const
{
    ofstream dot(fileName);
    writeGraphviz(dot);
}



void AssemblyGraph::writeGraphviz(ostream& dot) const
{
    const AssemblyGraph& assemblyGraph = *this;

    dot << "digraph AssemblyGraph {\n";



    // Write the vertices.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
    	const AssemblyGraphVertex& vertex = assemblyGraph[v];
    	dot <<
    		vertex.id <<
    		" [label=\"" << anchorIdToString(vertex.anchorId) << "\\n" << vertex.id;
    	if(vertex.vRc != null_vertex()) {
    	    dot << "\\n(" << assemblyGraph[vertex.vRc].id << ")";
    	}
    	dot << "\"]"
    	    ";\n";
    }



    // Write the edges.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
    	const vertex_descriptor v0 = source(e, assemblyGraph);
    	const vertex_descriptor v1 = target(e, assemblyGraph);
    	const AssemblyGraphVertex& vertex0 = assemblyGraph[v0];
    	const AssemblyGraphVertex& vertex1 = assemblyGraph[v1];
    	dot <<
    	    vertex0.id << "->" <<
    	    vertex1.id <<
    	    " [label=\"" << edge.id << "\\n" <<
    	    (edge.wasAssembled ? edge.sequenceLength() : edge.offset()) <<
    	    "\\n" << edge.size() <<
    	    "\"]"
    	    ";\n";
    }

    dot << "}\n";
}



// Assemble sequence for all edges.
void AssemblyGraph::assembleAll()
{
    writeMemoryStatistics("AssemblyGraph::assembleAll begins");

    const AssemblyGraph& assemblyGraph = *this;

    vector<edge_descriptor> edgesToBeAssembled;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        edgesToBeAssembled.push_back(e);
    }
    assemble(edgesToBeAssembled);

    writeMemoryStatistics("AssemblyGraph::assembleAll ends");
}



// Assemble sequence for the specified edge.
void AssemblyGraph::assemble(edge_descriptor e)
{
    vector<edge_descriptor> edgesToBeAssembled;
    edgesToBeAssembled.push_back(e);
    assemble(edgesToBeAssembled);
}



// Assemble sequence for step i of the specified edge.
// This is the lowest level sequence assembly function and is not multithreaded.
// It runs a LocalAssembly2 on the AnchorPair for that step.
void AssemblyGraph::assembleStep(edge_descriptor e, uint64_t i)
{
    AssemblyGraph& assemblyGraph = *this;
    AssemblyGraphEdge& edge = assemblyGraph[e];
    AssemblyGraphEdgeStep& step = edge[i];

    // cout << timestamp << " " << getPeakMemoryUsage() <<
    //     " Begin local assembly for edge " << edge.id << " " << " step " << i << endl;

    if(step.anchorPair.anchorIdA == step.anchorPair.anchorIdB) {
        step.sequence.clear();
        return;
    }



    // Let the local assembly use OrientedReadIds from the previous and next step.
    vector<OrientedReadId> additionalOrientedReadIds;
    if(i > 0) {
        std::ranges::copy(edge[i - 1].anchorPair.orientedReadIds, back_inserter(additionalOrientedReadIds));
    }
    if(i < edge.size() - 1) {
        std::ranges::copy(edge[i + 1].anchorPair.orientedReadIds, back_inserter(additionalOrientedReadIds));
    }
    deduplicate(additionalOrientedReadIds);

    ostream html(0);
    try {

        // Combine the Oriented reads in the AnchorPair
        // and the additional OrientedReadIds.
        vector<OrientedReadId> orientedReadIds = additionalOrientedReadIds;
        const AnchorPair& anchorPair = edge[i].anchorPair;
        std::ranges::copy(anchorPair.orientedReadIds, back_inserter(orientedReadIds));
        deduplicate(orientedReadIds);

        if(false) {
            LocalAssembly6 localAssembly(
                anchors,
                anchorPair.anchorIdA,
                anchorPair.anchorIdB,
                html,
                orientedReadIds);
            step.sequence = localAssembly.sequence;

        } else {

            LocalAssembly7::Options localAssembly7Options;
            localAssembly7Options.useMsaRepair = options.useMsaRepair;
            LocalAssembly7 localAssembly(
                localAssembly7Options,
                anchors,
                homopolymerModelPointer,
                anchorPair.anchorIdA,
                anchorPair.anchorIdB,
                html,
                orientedReadIds);
            step.sequence = localAssembly.sequence;

            if(not localAssembly.success) {
                std::lock_guard<std::mutex> lock(mutex);
                throw runtime_error("Local assembly for segment " + to_string(edge.id) +
                    " step " + to_string(i) + " failed.");
            }
        }

        if(step.sequence.empty()) {
            std::lock_guard<std::mutex> lock(mutex);
            throw runtime_error("Local assembly for segment " + to_string(edge.id) +
                " step " + to_string(i) + ": empty sequence assembled.");
        }

    } catch(const std::exception&) {
        std::lock_guard<std::mutex> lock(mutex);
        cout << "Error occurred assembling segment " <<
            assemblyGraph[e].id << " step " << i << endl;
        throw;
    }
}



// Assemble sequence for all edges in the edgesToBeAssembled vector.
// This fills in the stepsToBeAssembled with all steps of those edges,
// then assembles each of the steps in parallel.
void AssemblyGraph::assemble(const vector<edge_descriptor>& edgesToBeAssembled)
{
    const bool debug = false;
    AssemblyGraph& assemblyGraph = *this;

    // The homopolymer model must be consistent with the Options: present if
    // and only if the Options specify one (see
    // Assembler::createHomopolymerModel). All sequence assembly goes through
    // here.
    if(options.homopolymerModelName.empty()) {
        SHASTA2_ASSERT(not homopolymerModelPointer);
    } else {
        SHASTA2_ASSERT(homopolymerModelPointer);
    }

    performanceLog << timestamp << "Sequence assembly begins for " << edgesToBeAssembled.size() <<
        " assembly graph edges." << endl;
    if(debug) {
        cout << "Edges to be assembled:" << endl;
        for(const Segment segment: edgesToBeAssembled) {
            cout << "Segment  " << id(segment) << ": " << assemblyGraph[segment].size() << " steps." << endl;
        }
    }


    stepsToBeAssembled.clear();
    for(const edge_descriptor e: edgesToBeAssembled) {
        AssemblyGraphEdge& edge = assemblyGraph[e];
        for(uint64_t i=0; i<edge.size(); i++) {
            stepsToBeAssembled.push_back(make_pair(e, i));
        }
    }

    const uint64_t batchCount = 1;
    setupLoadBalancing(stepsToBeAssembled.size(), batchCount);
    runThreads(&AssemblyGraph::assembleThreadFunction, options.actualThreadCount());

    // Mark them as assembled.
    for(const edge_descriptor e: edgesToBeAssembled) {
        assemblyGraph[e].wasAssembled = true;
    }

    stepsToBeAssembled.clear();
    stepsToBeAssembled.shrink_to_fit();

    performanceLog << timestamp << "Sequence assembly ends." << endl;
}



void AssemblyGraph::assembleThreadFunction([[maybe_unused]] uint64_t threadId)
{
    AssemblyGraph& assemblyGraph = *this;

    // ofstream out("assembleThreadFunction-" + to_string(threadId) + ".csv");

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all assembly steps assigned to this batch.
        for(uint64_t j=begin; j!=end; j++) {

            if((j % 1000000) == 0) {
                std::lock_guard<std::mutex> lock(mutex);
                performanceLog << timestamp << "Starting sequence assembly step " << j << " of " <<
                    stepsToBeAssembled.size() << endl;
            }

            const auto& p = stepsToBeAssembled[j];
            const edge_descriptor e = p.first;
            const uint64_t i = p.second;
            AssemblyGraphEdge& edge = assemblyGraph[e];
            SHASTA2_ASSERT(i < edge.size());

            try {
                const auto t0 = steady_clock::now();
                // out << "Begin segment " << assemblyGraph[e].id << " step " << i << endl;
                assembleStep(e, i);
                // out << "End segment " << assemblyGraph[e].id << " step " << i << endl;
                const auto t1 = steady_clock::now();
                const double t01 = seconds(t1-t0);
                if(t01 > 10.) {
                    std::lock_guard<std::mutex> lock(mutex);
                    performanceLog << "Slow assembly step: segment " <<
                        assemblyGraph[e].id << " step " << i <<
                        " " << t01 << " s." << endl;

                }
            } catch(std::exception&) {
                std::lock_guard<std::mutex> lock(mutex);
                cout << "Error occurred during local assembly step: segment " <<
                    assemblyGraph[e].id << " step " << i << endl;
                throw;
            }
         }
    }
}



// Clear sequence from all steps of all edges.
void AssemblyGraph::clearSequence()
{
    AssemblyGraph& assemblyGraph = *this;

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        AssemblyGraphEdge& edge = assemblyGraph[e];
        edge.wasAssembled= false;
        for(AssemblyGraphEdgeStep& step: edge) {
            step.sequence.clear();
            step.sequence.shrink_to_fit();
        }
    }
}



void AssemblyGraphEdge::getSequence(vector<Base>& sequence) const
{
    sequence.clear();
    for(const auto& step: *this) {
        copy(step.sequence.begin(), step.sequence.end(), back_inserter(sequence));
    }
}



uint64_t AssemblyGraphEdge::sequenceLength() const
{
    SHASTA2_ASSERT(wasAssembled);

    uint64_t length = 0;
    for(const auto& step: *this) {
        length += step.sequence.size();
    }
    return length;
}



uint64_t AssemblyGraphEdge::length() const
{
    if(wasAssembled) {
        return sequenceLength();
    } else {
        return offset();
    }
}



// Compress linear chains of edges into a single edge.
uint64_t AssemblyGraph::compress()
{
    AssemblyGraph& assemblyGraph = *this;

    // Find linear chains of 2 or more edges.
    vector< std::list<edge_descriptor> > chains;
    findLinearChains(assemblyGraph, 2, chains);

    for(const auto& chain: chains) {
        compressLinearChain(chain);
    }

    return chains.size();
}



uint64_t AssemblyGraph::strandSymmetricCompress()
{
    performanceLog << timestamp << "AssemblyGraph::strandSymmetricCompress begins." << endl;
    const bool debug = false;

    AssemblyGraph& assemblyGraph = *this;

    // Find linear chains of 2 or more edges.
    vector< std::list<edge_descriptor> > chains;
    findLinearChains(assemblyGraph, 2, chains);

    // To find pairs of reverse complement chais,
    // index the edges by the chain they belong to.
    // This is very strict.
    std::map<edge_descriptor, uint64_t> chainMap;
    for(uint64_t chainId=0; chainId<chains.size(); chainId++) {
        const std::list<edge_descriptor>& chain = chains[chainId];
        for(const edge_descriptor e: chain) {
            SHASTA2_ASSERT(not chainMap.contains(e));
            chainMap.insert(make_pair(e, chainId));
        }
    }

    // Now for each chain we can find its reverse complement.
    vector<uint64_t> chainTable(chains.size(), invalid<uint64_t>);
    vector<uint64_t> v;
    for(uint64_t chainId=0; chainId<chains.size(); chainId++) {
        // cout << "Working on chain " << chainId << endl;
        const std::list<edge_descriptor>& chain = chains[chainId];
        v.clear();
        for(const edge_descriptor e: chain) {
            const edge_descriptor eRc = assemblyGraph[e].eRc;
            SHASTA2_ASSERT(eRc != assemblyGraphNullEdge);
            const uint64_t chainIdRc = chainMap.at(eRc);
            v.push_back(chainIdRc);
        }
        deduplicate(v);
        SHASTA2_ASSERT(v.size() == 1);
        chainTable[chainId] = v.front();
    }

    // Sanity check.
    for(uint64_t chainId=0; chainId<chains.size(); chainId++) {
        const uint64_t chainIdRc = chainTable[chainId];
        SHASTA2_ASSERT(chainIdRc != chainId);
        SHASTA2_ASSERT(chainTable[chainIdRc] == chainId);
    }

    for(uint64_t chainId=0; chainId<chains.size(); chainId++) {
        const uint64_t chainIdRc = chainTable[chainId];
        if(chainId < chainIdRc) {
            auto& chain = chains[chainId];
            auto& chainRc = chains[chainIdRc];

            if(debug) {
                cout << "Working on chain";
                for(const edge_descriptor e: chain) {
                    cout << " " << assemblyGraph[e].id;
                }
                cout << endl;
                cout << "Reverse complemented chain is ";
                for(const edge_descriptor e: chainRc) {
                    cout << " " << assemblyGraph[e].id;
                }
                cout << endl;
            }

            // For the special case of an isolated circular chain,
            // we need to make sure that the reverse complemented
            // chain begins/ends at a vertex which is the
            // reverse complement of the begin/end of the chain.
            // We do this by rotating the reverse complemented chain.
            {
                const edge_descriptor e0 = chain.front();
                const edge_descriptor e1 = chain.back();
                const vertex_descriptor v0 = source(e0, assemblyGraph);
                const vertex_descriptor v1 = target(e1, assemblyGraph);
                if(v0 == v1) {
                    // It is a circular chain.
                    if((in_degree(v0, assemblyGraph) == 1) and (out_degree(v0, assemblyGraph) == 1)) {
                        // It is isolated.
                        const edge_descriptor e0Rc = assemblyGraph[e0].eRc;
                        const auto it = find(chainRc.rbegin(), chainRc.rend(), e0Rc);
                        SHASTA2_ASSERT(it != chainRc.rend());
                        std::rotate(chainRc.rbegin(), it, chainRc.rend());

                        if(debug) {
                            cout << "After rotation, reverse complemented chain is ";
                            for(const edge_descriptor e: chainRc) {
                                cout << " " << assemblyGraph[e].id;
                            }
                            cout << endl;
                        }
                    }
                }
            }


            // Sanity check on the reverse complemented chain.
            const uint64_t n = chain.size();
            SHASTA2_ASSERT(chainRc.size() == n);
            auto it = chain.begin();
            auto itRc = chainRc.rbegin();
            for(uint64_t i=0; i<n; i++) {
                const edge_descriptor e = *it++;
                const edge_descriptor eRc = *(itRc++);
                SHASTA2_ASSERT(assemblyGraph[e].eRc == eRc);
                SHASTA2_ASSERT(assemblyGraph[eRc].eRc == e);
            }


            const edge_descriptor e = compressLinearChain(chains[chainId]);
            if(debug) {
                cout << "Created new edge " << assemblyGraph[e].id <<
                    ", length " << assemblyGraph[e].length() <<
                    ", source " << assemblyGraph[source(e, assemblyGraph)].id <<
                    ", target " << assemblyGraph[target(e, assemblyGraph)].id <<
                    endl;
            }
            const edge_descriptor eRc = createReverseComplementEdge(e);
            if(debug) {
                cout << "Created new reverse complement edge " << assemblyGraph[eRc].id <<
                    ", length " << assemblyGraph[eRc].length() <<
                    ", source " << assemblyGraph[source(eRc, assemblyGraph)].id <<
                    ", target " << assemblyGraph[target(eRc, assemblyGraph)].id <<
                    endl;
            }


            // The call to compressLinearChain removed the edges and
            // the internal vertices of chainId. Do the same for chainIdRc.
            bool isFirst = true;
            for(const edge_descriptor e: chains[chainIdRc]) {
                if(isFirst) {
                    isFirst = false;
                } else {
                    const vertex_descriptor v = source(e, assemblyGraph);
                    boost::clear_vertex(v, assemblyGraph);
                    boost::remove_vertex(v, assemblyGraph);
                }
            }
        }
    }

    performanceLog << timestamp << "AssemblyGraph::strandSymmetricCompress ends." << endl;

    return chains.size();
}



#if 0
// OLD VERSION
// Compress linear chains of edges into a single edge.
uint64_t AssemblyGraph::compress()
{
    AssemblyGraph& assemblyGraph = *this;
    uint64_t compressCount = 0;

    // Find linear chains of 2 or more edges.
    vector< std::list<edge_descriptor> > chains;
    findLinearChains(assemblyGraph, 2, chains);

    for(const auto& chain: chains) {
        SHASTA2_ASSERT(chain.size() > 1);

        // Get the first and last edge of this chain.
        const edge_descriptor e0 = chain.front();
        const edge_descriptor e1 = chain.back();

        // Get the first and last edge of this chain.
        const vertex_descriptor v0 = source(e0, assemblyGraph);
        const vertex_descriptor v1 = target(e1, assemblyGraph);

        // Add the new edge.
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(v0, v1, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
        AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];

        // Concatenate the steps of all the edges in the chain.
        for(const edge_descriptor e: chain) {
            const AssemblyGraphEdge& edge = assemblyGraph[e];
            copy(edge.begin(), edge.end(), back_inserter(edgeNew));
        }



        // Minimal debug output.
        if(compressDebugLevel >= 1) {
            cout << "Compress " << assemblyGraph[chain.front()].id << "..." <<
                assemblyGraph[chain.back()].id <<
                " into " << edgeNew.id << endl;
        }



        // Compact debug output.
        if(compressDebugLevel >= 2) {
            cout << "Compress";
            for(const edge_descriptor e: chain) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << " into " << edgeNew.id << endl;
        }



        // Detailed debug output.
        if(compressDebugLevel >= 3) {
            uint64_t stepCount = 0;
            for(const edge_descriptor e: chain) {
                const AssemblyGraphEdge& edge = assemblyGraph[e];
                const uint64_t stepBegin = stepCount;
                const uint64_t stepEnd = stepBegin + edge.size();
                cout << edge.id << " " << edge.size() << " steps become " <<
                    edgeNew.id << " steps " << stepBegin << "-" << stepEnd << endl;
                stepCount = stepEnd;
            }
        }



        // Now we can remove the edges of the chain and its internal vertices.
        bool isFirst = true;
        for(const edge_descriptor e: chain) {
            if(isFirst) {
                isFirst = false;
            } else {
                const vertex_descriptor v = source(e, assemblyGraph);
                boost::clear_vertex(v, assemblyGraph);
                boost::remove_vertex(v, assemblyGraph);
            }
        }

        ++compressCount;
    }

    return compressCount;
}
#endif



AssemblyGraph::edge_descriptor AssemblyGraph::compressLinearChain(const std::list<edge_descriptor>& chain)
{
    AssemblyGraph& assemblyGraph = *this;
    SHASTA2_ASSERT(chain.size() > 1);

    // Get the first and last edge of this chain.
    const edge_descriptor e0 = chain.front();
    const edge_descriptor e1 = chain.back();

    // Get the first and last vertex of this chain.
    const vertex_descriptor v0 = source(e0, assemblyGraph);
    const vertex_descriptor v1 = target(e1, assemblyGraph);

    // Add the new edge.
    edge_descriptor eNew;
    tie(eNew, ignore) = add_edge(v0, v1, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
    AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];

    // Concatenate the steps of all the edges in the chain.
    for(const edge_descriptor e: chain) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        copy(edge.begin(), edge.end(), back_inserter(edgeNew));
    }



    // Minimal debug output.
    if(compressDebugLevel >= 1) {
        cout << "Compress " << assemblyGraph[chain.front()].id << "..." <<
            assemblyGraph[chain.back()].id <<
            " into " << edgeNew.id << endl;
    }



    // Compact debug output.
    if(compressDebugLevel >= 2) {
        cout << "Compress";
        for(const edge_descriptor e: chain) {
            cout << " " << assemblyGraph[e].id;
        }
        cout << " into " << edgeNew.id << endl;
    }



    // Detailed debug output.
    if(compressDebugLevel >= 3) {
        uint64_t stepCount = 0;
        for(const edge_descriptor e: chain) {
            const AssemblyGraphEdge& edge = assemblyGraph[e];
            const uint64_t stepBegin = stepCount;
            const uint64_t stepEnd = stepBegin + edge.size();
            cout << edge.id << " " << edge.size() << " steps become " <<
                edgeNew.id << " steps " << stepBegin << "-" << stepEnd << endl;
            stepCount = stepEnd;
        }
    }



    // Now we can remove the edges of the chain and its internal vertices.
    bool isFirst = true;
    for(const edge_descriptor e: chain) {
        if(isFirst) {
            isFirst = false;
        } else {
            const vertex_descriptor v = source(e, assemblyGraph);
            boost::clear_vertex(v, assemblyGraph);
            boost::remove_vertex(v, assemblyGraph);
        }
    }

    return eNew;
}



void AssemblyGraph::save(ostream& s) const
{
    boost::archive::binary_oarchive archive(s);
    archive << *this;
}



void AssemblyGraph::load(istream& s)
{
    boost::archive::binary_iarchive archive(s);
    archive >> *this;
}



void AssemblyGraph::save(const string& stage) const
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
    const string name = largeDataName("AssemblyGraph-" + stage);
    MemoryMapped::Vector<char> data;
    data.createNew(name, largeDataPageSize);
    data.resize(dataString.size());
    const char* begin = dataString.data();
    const char* end = begin + dataString.size();
    copy(begin, end, data.begin());
}



void AssemblyGraph::load(const string& assemblyStage)
{
    // Access the binary data.
    MemoryMapped::Vector<char> data;
    try {
        const string name = largeDataName("AssemblyGraph-" + assemblyStage);
        data.accessExistingReadOnly(name);
    } catch (std::exception&) {
        throw runtime_error("Assembly graph at stage " + assemblyStage +
            " is not available.");
    }
    const string dataString(data.begin(), data.size());

    // Load it from here.
    std::istringstream s(dataString);
    try {
        load(s);
    } catch(std::exception& e) {
        throw runtime_error("Error reading assembly graph at stage " + assemblyStage +
            ": " + e.what());
    }
}



// Write a csv file that can be loaded in Bandage.
void AssemblyGraph::writeCsv(const string& fileName) const
{
    ofstream csv(fileName);
    writeCsv(csv);
}



// Write a csv file that can be loaded in Bandage.
void AssemblyGraph::writeCsv(ostream& csv) const
{
    const AssemblyGraph& assemblyGraph = *this;

    csv << "Segment,Reverse complement segment,Number of steps,Average coverage,Estimated length,Actual length,Annotation,\n";
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        const uint64_t coverage = uint64_t(std::round(edge.averageCoverage()));
        csv << edge.id << ",";
        if(edge.eRc != assemblyGraphNullEdge) {
            csv << assemblyGraph[edge.eRc].id;
        }
        csv <<
            "," <<
            edge.size() << "," <<
            coverage << "," <<
            edge.offset() << ",";
        if(edge.wasAssembled) {
            csv << edge.sequenceLength();
        }
        csv << "," << edge.annotation << ",\n";
    }
}



void AssemblyGraph::writeSequenceLengthByCoverageCsv(const string& fileName) const
{
    ofstream csv(fileName);
    writeSequenceLengthByCoverageCsv(csv);
}



void AssemblyGraph::writeSequenceLengthByCoverageCsv(ostream& csv) const
{
    const AssemblyGraph& assemblyGraph = *this;

    // Compute the total sequence length for each coverage value.
    vector<uint64_t> v; // Sequence
    vector<uint64_t> n; // Number of segments
    uint64_t minCoverage = std::numeric_limits<uint64_t>::max();
    uint64_t maxCoverage =0;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        const uint64_t length = edge.length();
        if(length == 0) {
            continue;
        }
        const uint64_t coverage = uint64_t(std::round(edge.averageCoverage()));
        minCoverage = min(minCoverage, coverage);
        maxCoverage = max(maxCoverage, coverage);

        if(v.size() <= coverage) {
            n.resize(coverage + 1, 0);
            v.resize(coverage + 1, 0);
        }
        ++(n[coverage]);
        v[coverage] += length;
    }

    // Write it out.
    csv << "Coverage,Number of segments,Sequence length,"
        "Cunmulative number of segments,Cumulative sequence length\n";
    uint64_t cumulativeNumberOfSegments = 0;
    uint64_t cumulativeSequenceLength = 0;
    for(uint64_t coverage=minCoverage; coverage<=maxCoverage; coverage++) {
        cumulativeNumberOfSegments += n[coverage];
        cumulativeSequenceLength += v[coverage];
        csv <<
            coverage << "," <<
            n[coverage] << "," <<
            v[coverage] << "," <<
            cumulativeNumberOfSegments << "," <<
            cumulativeSequenceLength << "," <<
            "\n";
    }

}



void AssemblyGraph::writeDetailsCsv(const string& fileName) const
{
    ofstream csv(fileName);
    writeDetailsCsv(csv);
}



void AssemblyGraph::writeDetailsCsv(ostream& csv) const
{
    const AssemblyGraph& assemblyGraph = *this;

    // Write a csv header.
    csv << "SegmentId,Step,AnchorIdA,AnchorIdB,Coverage,Estimated length,Actual length,Sequence begin,Sequence end,\n";

    // Loop over all segments (AssemblyGraph edges).
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        uint64_t begin = 0;

        // Loop over all steps of this segment.
        for(uint64_t stepId=0; stepId<edge.size(); stepId++) {
            const AssemblyGraphEdgeStep& step = edge[stepId];
            const uint64_t end = begin + (edge.wasAssembled ? step.sequence.size() : step.offset);

            csv <<
                edge.id << "," <<
                stepId << "," <<
                anchorIdToString(step.anchorPair.anchorIdA) << "," <<
                anchorIdToString(step.anchorPair.anchorIdB) << "," <<
                step.anchorPair.orientedReadIds.size() << "," <<
                step.offset << ",";
            if(edge.wasAssembled) {
                csv << step.sequence.size();
            }
            csv <<
                "," <<
                begin << "," <<
                end << "," <<
                "\n";

            begin = end;
        }
    }
}



// Find the non-trivial strongly connected components.
// Each component is stored with vertices sorted to permit binary searches.
void AssemblyGraph::findStrongComponents(
    vector< vector<vertex_descriptor> >& strongComponents) const
{
    const AssemblyGraph& assemblyGraph = *this;

    // Map the vertices to integers.
    uint64_t vertexIndex = 0;
    std::map<vertex_descriptor, uint64_t> vertexMap;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexMap.insert({v, vertexIndex++});
    }

    // Compute strong components.
    std::map<vertex_descriptor, uint64_t> componentMap;
    boost::strong_components(
        assemblyGraph,
        boost::make_assoc_property_map(componentMap),
        boost::vertex_index_map(boost::make_assoc_property_map(vertexMap)));

    // Gather the vertices in each strong component.
    std::map<uint64_t, vector<vertex_descriptor> > componentVertices;
    for(const auto& p: componentMap) {
        componentVertices[p.second].push_back(p.first);
    }

    strongComponents.clear();
    for(const auto& p: componentVertices) {
        const vector<vertex_descriptor>& component = p.second;
        if(component.size() > 1) {
            strongComponents.push_back(component);
            sort(strongComponents.back().begin(), strongComponents.back().end());
        }
    }


}



// This creates a csv file that can be loaded in bandage to see
// the strongly connected components.
void AssemblyGraph::colorStrongComponents() const
{
    const AssemblyGraph& assemblyGraph = *this;

    vector< vector<vertex_descriptor> > strongComponents;
    findStrongComponents(strongComponents);

    ofstream csv("StrongComponents.csv");
    csv << "Id,Color,Component\n";

    // Loop over the non-trivial strongly connected components.
    for(uint64_t i=0; i<strongComponents.size(); i++) {
        const vector<vertex_descriptor>& strongComponent = strongComponents[i];

        for(const vertex_descriptor v0: strongComponent) {
            BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
                const vertex_descriptor v1 = target(e, assemblyGraph);
                if(std::binary_search(strongComponent.begin(), strongComponent.end(), v1)) {
                    csv << assemblyGraph[e].id << ",Green," << i << "\n";
                }
            }
        }

    }


}




double AssemblyGraphEdge::averageCoverage() const
{
    if(empty()) {
        return 0.;
    }

    uint64_t sum = 0;
    for(const AssemblyGraphEdgeStep& step: *this) {
        sum += step.anchorPair.orientedReadIds.size();
    }

    return double(sum) / double(size());
}



double AssemblyGraphEdge::lengthWeightedAverageCoverage() const
{
    uint64_t sum0 = 0;
    uint64_t sum1 = 0;
    for(const AssemblyGraphEdgeStep& step: *this) {
        const uint64_t length = wasAssembled ? step.sequence.size() : step.offset;
        const uint64_t coverage = step.anchorPair.orientedReadIds.size();
        sum0 += length;
        sum1 += length * coverage;
    }

    return double(sum1) / double(sum0);
}



// Compute oriented read journeys in the AssemblyGraph.
void AssemblyGraph::computeJourneys()
{
    const uint64_t orientedReadCount = journeys.size();

    AssemblyGraph& assemblyGraph = *this;


    // Compute uncompressed journeys for all oriented reads.
    class AssemblyGraphJourneyEntry {
    public:
        edge_descriptor e;
        uint64_t stepId;

        // These are positions in the standard Journeys on the Anchors.
        uint64_t positionInJourneyA;
        uint64_t positionInJourneyB;

        bool operator<(const AssemblyGraphJourneyEntry& that) const {
            return positionInJourneyA + positionInJourneyB <
                that.positionInJourneyA + that.positionInJourneyB;
        }
    };
    vector< vector<AssemblyGraphJourneyEntry> > assemblyGraphJourneys(orientedReadCount);

    // Loop over AssemblyGraph edges.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        // Loop over steps of this edge.
        for(uint64_t stepId=0; stepId<edge.size(); stepId++) {
            const AssemblyGraphEdgeStep& step = edge[stepId];

            // Locate the anchors for this step.
            const AnchorPair& anchorPair = step.anchorPair;
            const AnchorId anchorIdA = anchorPair.anchorIdA;
            const AnchorId anchorIdB = anchorPair.anchorIdB;
            const Anchor anchorA = anchors[anchorIdA];
            const Anchor anchorB = anchors[anchorIdB];

            // Loop over OrientedReadIds of this step.
            auto itA = anchorA.begin();
            auto itB = anchorB.begin();
            for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {

                // Locate this OrientedReadId in the two anchors.
                for(; (itA != anchorA.end()) and (itA->orientedReadId != orientedReadId); ++itA) {}
                SHASTA2_ASSERT(itA != anchorA.end());
                SHASTA2_ASSERT(itA->orientedReadId == orientedReadId);
                const AnchorMarkerInfo& infoA = *itA;
                for(; (itB != anchorB.end()) and (itB->orientedReadId != orientedReadId); ++itB) {}
                SHASTA2_ASSERT(itB != anchorB.end());
                const AnchorMarkerInfo& infoB = *itB;
                SHASTA2_ASSERT(itB->orientedReadId == orientedReadId);

                const uint32_t positionInJourneyA = infoA.positionInJourney;
                const uint32_t positionInJourneyB = infoB.positionInJourney;

                AssemblyGraphJourneyEntry entry;
                entry.e = e;
                entry.stepId = stepId;
                entry.positionInJourneyA = positionInJourneyA;
                entry.positionInJourneyB = positionInJourneyB;

                assemblyGraphJourneys[orientedReadId.getValue()].push_back(entry);

            }

        }
    }


    // Sort the journeys.
    for(vector<AssemblyGraphJourneyEntry>& v: assemblyGraphJourneys) {
        sort(v.begin(), v.end());
    }



    // Create the compressed journeys.
    // Here, we only consider transitions between AssemblyGraph edges.
    // Compressed journeys consisting of only one edge are considered empty.
    compressedJourneys.clear();
    compressedJourneys.resize(orientedReadCount);
    for(ReadId orientedReadIdValue=0; orientedReadIdValue<orientedReadCount; orientedReadIdValue++) {
        const vector<AssemblyGraphJourneyEntry>& assemblyGraphJourney = assemblyGraphJourneys[orientedReadIdValue];
        vector<edge_descriptor>& compressedJourney = compressedJourneys[orientedReadIdValue];

        for(uint64_t i1=0; i1<assemblyGraphJourney.size(); i1++) {
            const AssemblyGraphJourneyEntry& entry1 = assemblyGraphJourney[i1];
            const edge_descriptor e1 = entry1.e;

            if(i1 == 0) {
                compressedJourney.push_back(e1);
            } else {
                const uint64_t i0 = i1 - 1;
                const AssemblyGraphJourneyEntry& entry0 = assemblyGraphJourney[i0];
                const edge_descriptor e0 = entry0.e;
                if(e1 != e0) {
                    compressedJourney.push_back(e1);
                }
            }
        }
        if(compressedJourney.size() == 1) {
            compressedJourney.clear();
        }
    }



    // Write out the compressed journeys.
    {
        ofstream csv("AssemblyGraphCompressedJourneys.csv");
        for(ReadId orientedReadIdValue=0; orientedReadIdValue<orientedReadCount; orientedReadIdValue++) {
            const OrientedReadId orientedReadId = OrientedReadId::fromValue(orientedReadIdValue);
            const vector<edge_descriptor>& compressedJourney = compressedJourneys[orientedReadIdValue];
            csv << orientedReadId << ",";
            for(const edge_descriptor e: compressedJourney) {
                csv << assemblyGraph[e].id << ",";
            }
            csv << "\n";
        }
    }

    {
        ofstream csv("AssemblyGraphJourneys.csv");
        csv << "OrientedReadId,Segment,Step,PositionInJourneyA,PositionInJourneyB\n";
        for(ReadId orientedReadIdValue=0; orientedReadIdValue<orientedReadCount; orientedReadIdValue++) {
            const OrientedReadId orientedReadId = OrientedReadId::fromValue(orientedReadIdValue);
            const vector<AssemblyGraphJourneyEntry>& assemblyGraphJourney = assemblyGraphJourneys[orientedReadIdValue];
            for(const AssemblyGraphJourneyEntry& entry: assemblyGraphJourney) {
                csv << orientedReadId << ",";
                csv << assemblyGraph[entry.e].id << ",";
                csv << entry.stepId << ",";
                csv << entry.positionInJourneyA << ",";
                csv << entry.positionInJourneyB << "\n";
            }
        }
    }
}



// The detangling process can generate empty edges (edges without steps).
// This removes them by collapsing the vertices they join.
void AssemblyGraph::removeEmptyEdges()
{
    AssemblyGraph& assemblyGraph = *this;
    const bool debug = false;

    // We need to find groups of vertices that need to be collapsed together.
    // Usually it will be just two vertices to be collapsed together,
    // but it could be larger groups if there are adjacent empty edges.
    // So we need to find the connected components generated by the empty edges.

    // Map vertices to integer.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    vector<vertex_descriptor> vertexTable;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexIndexMap.insert(make_pair(v, vertexIndex++));
        vertexTable.push_back(v);
    }

    // Initialize the disjoint sets data structure.
    const uint64_t n = vertexIndexMap.size();
    // Initialize the disjoint sets data structure.
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }

    // Loop over the empty edges.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        if(assemblyGraph[e].empty()) {
            const vertex_descriptor v0 = source(e, assemblyGraph);
            const vertex_descriptor v1 = target(e, assemblyGraph);
            const uint64_t i0 = vertexIndexMap[v0];
            const uint64_t i1 = vertexIndexMap[v1];
            disjointSets.union_set(i0, i1);
            edgesToBeRemoved.push_back(e);
        }
    }

    // Gather the vertices in each connected component.
    vector< vector<uint64_t> > groups(n);
    for(uint64_t i=0; i<n; i++) {
        groups[disjointSets.find_set(i)].push_back(i);
    }

    // Collapse each group into a single vertex.
    for(const vector<uint64_t>& group: groups) {
        if(group.size() < 2) {
            continue;
        }
        if(debug) {
            cout << "Found a group of " << group.size() << " vertices to be collapsed:";
            for(const uint64_t vertexIndex: group) {
                const vertex_descriptor v = vertexTable[vertexIndex];
                cout << " " << assemblyGraph[v].id << "/" << anchorIdToString(assemblyGraph[v].anchorId);
            }
            cout << endl;
        }

        // Create the collapsed vertex.
        const AnchorId anchorId = assemblyGraph[vertexTable[group.front()]].anchorId;
        const vertex_descriptor vNew = add_vertex(
             AssemblyGraphVertex(anchorId, assemblyGraph.nextVertexId++), assemblyGraph);



        // For each vertex in this group, reroute all incoming/outgoing non-empty edges
        // to/from the new vertex.
        for(const uint64_t vertexIndex: group) {
            const vertex_descriptor v = vertexTable[vertexIndex];

            // Loop over non-empty incoming edges.
            BGL_FORALL_INEDGES(v, eOld, assemblyGraph, AssemblyGraph) {
                AssemblyGraphEdge& oldEdge = assemblyGraph[eOld];
                if(oldEdge.empty()) {
                    continue;
                }
                const vertex_descriptor u = source(eOld, assemblyGraph);

                // Create the new edge, with target vNew.
                AssemblyGraph::edge_descriptor eNew;
                tie(eNew, ignore) = add_edge(u, vNew, assemblyGraph);
                AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
                newEdge.id = oldEdge.id;
                newEdge.swapSteps(oldEdge);

                edgesToBeRemoved.push_back(eOld);
            }


            // Loop over non-empty outgoing edges.
            BGL_FORALL_OUTEDGES(v, eOld, assemblyGraph, AssemblyGraph) {
                AssemblyGraphEdge& oldEdge = assemblyGraph[eOld];
                if(oldEdge.empty()) {
                    continue;
                }
                const vertex_descriptor u = target(eOld, assemblyGraph);

                // Create the new edge, with source vNew.
                AssemblyGraph::edge_descriptor eNew;
                tie(eNew, ignore) = add_edge(vNew, u, assemblyGraph);
                AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
                newEdge.id = oldEdge.id;
                newEdge.swapSteps(oldEdge);

                edgesToBeRemoved.push_back(eOld);
            }
        }
    }

    deduplicate(edgesToBeRemoved);
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, assemblyGraph);
    }

    // Remove any vertices that were left isolated.
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        if((in_degree(v, assemblyGraph) == 0) and (out_degree(v, assemblyGraph) == 0)) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        boost::remove_vertex(v, assemblyGraph);
    }
}



// Return true if the two specified edges can be connected for assembly.
// If necessary, this constructs the RestrictedAnchorGraph between the two edges
// and checks that all is good.
bool AssemblyGraph::canConnect(edge_descriptor e0, edge_descriptor e1, bool deep) const
{
    const AssemblyGraph& assemblyGraph = *this;

    // If the last AnchorId of e0 is the same as the first AnchorId of v1,
    // They can be trivially connected.
    const AssemblyGraph::vertex_descriptor v0 = target(e0, assemblyGraph);
    const AssemblyGraph::vertex_descriptor v1 = source(e1, assemblyGraph);
    if(v0 == v1) {
        return true;
    }
    const AnchorId anchorId0 = assemblyGraph[v0].anchorId;
    const AnchorId anchorId1 = assemblyGraph[v1].anchorId;
    if(anchorId0 == anchorId1) {
        return true;
    }



    // In the general case, we need to construct the RestrictedAnchorGraph
    // between e0 and e1. The RestrictedAnchorGraph requires a TangleMatrix,
    // so we create a 1x1 tangle matrix with e0 as the entrance and e1 as the exit.
    ostream html(0);
    TangleMatrix tangleMatrix(assemblyGraph, {e0}, {e1}, html);



    // Create the RestrictedAnchorGraph. If this fails, return false.
    try {
        ostream html(0);

        // Create the RestrictedAnchorGraph.
        RestrictedAnchorGraph restrictedAnchorGraph(anchors, journeys, tangleMatrix, 0, 0, html, deep);
        vector<RestrictedAnchorGraph::edge_descriptor> longestPath;

        // Check that we have a path.
        restrictedAnchorGraph.findOptimalPath(anchorId0, anchorId1, longestPath);

        // Check coverage.
        uint64_t minCoverage = std::numeric_limits<uint64_t>::max();
        for(const RestrictedAnchorGraph::edge_descriptor e: longestPath) {
            const RestrictedAnchorGraphEdge& edge = restrictedAnchorGraph[e];
            minCoverage = min(minCoverage, edge.anchorPair.size());
        }
        if(minCoverage == 0) {
            return false;
        }

    } catch (const std::exception& e) {
        return false;
    }



    // If we get here, all is good.
    return true;
}



void AssemblyGraph::removeIsolatedVertices()
{
    AssemblyGraph& assemblyGraph = *this;


    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        if((out_degree(v, assemblyGraph) == 0) and (in_degree(v, assemblyGraph) == 0)) {
            verticesToBeRemoved.push_back(v);
        }
    }

    for(const vertex_descriptor v: verticesToBeRemoved) {
        boost::remove_vertex(v, assemblyGraph);
    }

}



uint64_t AssemblyGraph::countIsolatedVertices() const
{
    const AssemblyGraph& assemblyGraph = *this;

    uint64_t n = 0;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        if((out_degree(v, assemblyGraph) == 0) and (in_degree(v, assemblyGraph) == 0)) {
            ++n;
        }
    }
    return n;
}


// Remove connected components with a low N50.
void AssemblyGraph::removeLowN50Components()
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minComponentN50 = 100000;

    AssemblyGraph& assemblyGraph = *this;
    const bool debug = false;

    // Map the vertices to integers.
    uint64_t vertexIndex = 0;
    vector<vertex_descriptor> vertexTable;                  // Map integers to vertices
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;   // Map vertices to integers
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexTable.push_back(v);
        vertexIndexMap.insert({v, vertexIndex++});
    }



    // Compute connected components.
    // We can't use boost::connected_components because that only
    // supports undirected graphs.
    DisjointSets disjointSets(vertexIndexMap.size());
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);
        disjointSets.unionSet(vertexIndexMap[v0], vertexIndexMap[v1]);
    }
    vector< vector<uint64_t> > components;
    disjointSets.gatherComponents(1, components);
    if(debug) {
        cout << "Found " << components.size() << " connected components of the AssemblyGraph." << endl;
    }



    // Loop over connected components.
    for(const vector<uint64_t>& component: components) {

        // Gather the lengths of all edges in this component.
        vector<uint64_t> lengths;
        for(uint64_t vIndex0: component) {
            const vertex_descriptor v0 = vertexTable[vIndex0];
            BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
                const vertex_descriptor v1 = target(e, assemblyGraph);
                const uint64_t vIndex1 = vertexIndexMap[v1];
                SHASTA2_ASSERT(disjointSets.findSet(vIndex0) == disjointSets.findSet(vIndex1));
                lengths.push_back(assemblyGraph[e].length());
            }
        }
        const uint64_t totalLength = std::accumulate(lengths.begin(), lengths.end(), 0);

        // Find the N50.
        std::ranges::sort(lengths, std::greater<uint64_t>());
        uint64_t cumulativeLength = 0;
        uint64_t n50 = 0;
        for(const uint64_t length: lengths) {
            cumulativeLength += length;
            if(2 * cumulativeLength >= totalLength) {
                n50 = length;
                break;
            }
        }


        // Decide if we keep this component.
        const bool keep = (n50 >= minComponentN50);
        if(debug) {
            if(keep) {
                cout << "Keeping";
            } else {
                cout << "Discarding";
            }
            cout << " a connected component with " << lengths.size() <<
                " segments, total length " << totalLength <<
                ", N50 " << n50 << endl;
        }

        // If not keeping this component, remove its vertices and edges.
        if(not keep) {
            for(uint64_t vIndex: component) {
                const vertex_descriptor v = vertexTable[vIndex];
                boost::clear_vertex(v, assemblyGraph);
                boost::remove_vertex(v, assemblyGraph);
            }

        }

    }

}



// Remove connected components with small total length.
void AssemblyGraph::removeSmallComponents()
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minComponentTotalLength = 10000;

    AssemblyGraph& assemblyGraph = *this;
    const bool debug = false;

    // Map the vertices to integers.
    uint64_t vertexIndex = 0;
    vector<vertex_descriptor> vertexTable;                  // Map integers to vertices
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;   // Map vertices to integers
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexTable.push_back(v);
        vertexIndexMap.insert({v, vertexIndex++});
    }



    // Compute connected components.
    // We can't use boost::connected_components because that only
    // supports undirected graphs.
    DisjointSets disjointSets(vertexIndexMap.size());
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);
        disjointSets.unionSet(vertexIndexMap[v0], vertexIndexMap[v1]);
    }
    vector< vector<uint64_t> > components;
    disjointSets.gatherComponents(1, components);
    if(debug) {
        cout << "Found " << components.size() << " connected components of the AssemblyGraph." << endl;
    }



    // Loop over connected components.
    for(const vector<uint64_t>& component: components) {

        // Compute the total length of all edges in this component.
        uint64_t totalLength = 0;
        for(uint64_t vIndex0: component) {
            const vertex_descriptor v0 = vertexTable[vIndex0];
            BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
                const vertex_descriptor v1 = target(e, assemblyGraph);
                const uint64_t vIndex1 = vertexIndexMap[v1];
                SHASTA2_ASSERT(disjointSets.findSet(vIndex0) == disjointSets.findSet(vIndex1));
                totalLength += assemblyGraph[e].length();
            }
        }

        // Decide if we keep this component.
        const bool keep = (totalLength >= minComponentTotalLength);

        // If not keeping this component, remove its vertices and edges.
        if(not keep) {
            for(uint64_t vIndex: component) {
                const vertex_descriptor v = vertexTable[vIndex];
                boost::clear_vertex(v, assemblyGraph);
                boost::remove_vertex(v, assemblyGraph);
            }

        }

    }

}



void AssemblyGraph::prune()
{
    while(pruneIteration());
}



// Prune can be used on an uncompressed AssemblyGraph.
// So we find linear chains that are leafs, and prune
// the ones that are too short.
uint64_t AssemblyGraph::pruneIteration()
{
    AssemblyGraph& assemblyGraph = *this;

    // Find linear chains.
    using Chain = std::list<edge_descriptor>;
    vector<Chain> chains;
    findLinearChains(assemblyGraph, 1, chains);



    // Find the ones to be removed.
    vector<uint64_t> chainsToBeRemoved;
    for(uint64_t chainId=0; chainId<chains.size(); chainId++) {
        const Chain& chain = chains[chainId];
        SHASTA2_ASSERT(not chain.empty());

        // Find the first and last vertex.
        const vertex_descriptor v0 = source(chain.front(), assemblyGraph);
        const vertex_descriptor v1 = target(chain.back(), assemblyGraph);

        // If not a leaf, skip it.
        const bool isLeaf =
            (in_degree(v0, assemblyGraph) == 0) or
            (out_degree(v1, assemblyGraph) == 0);
        if(not isLeaf) {
            continue;
        }

        // If not short, skip it.
        uint64_t length = 0;
        for(const edge_descriptor e: chain) {
            length += assemblyGraph[e].length();
        }
        if(length > options.pruneLength) {
            continue;
        }

        // This chain is a short leaf, so we will remove it.
        chainsToBeRemoved.push_back(chainId);
    }


    // Remove them.
    for(const uint64_t chainId: chainsToBeRemoved) {
        const Chain& chain = chains[chainId];
        for(const edge_descriptor e: chain) {
            boost::remove_edge(e, assemblyGraph);
        }
    }

    return chainsToBeRemoved.size();
}



void AssemblyGraph::strictPrune()
{
    AssemblyGraph& assemblyGraph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);

        // If out-degree(v1)=0, check for
        // an edge v0->v2 such that out_degree(v2)>0.
        if(out_degree(v1, assemblyGraph) == 0) {
            bool found = false;
            BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
                const vertex_descriptor v2 = target(e, assemblyGraph);
                if(out_degree(v2, assemblyGraph) > 0) {
                    found = true;
                    break;
                }
            }
            if(found) {
                if(assemblyGraph[e].length() <= options.strictPruneLength) {
                    edgesToBeRemoved.push_back(e);
                }
                continue;
            }
        }

        // If in-degree(v0)=0, check for
        // an edge v2->v0 such that in_degree(v2)>0.
        if(in_degree(v0, assemblyGraph) == 0) {
            bool found = false;
            BGL_FORALL_INEDGES(v1, e, assemblyGraph, AssemblyGraph) {
                const vertex_descriptor v2 = source(e, assemblyGraph);
                if(in_degree(v2, assemblyGraph) > 0) {
                    found = true;
                    break;
                }
            }
            if(found) {
                if(assemblyGraph[e].length() <= options.strictPruneLength) {
                    edgesToBeRemoved.push_back(e);
                }
                continue;
            }
        }

    }


#if 0
    // EXPENSIVE CHECK.
    for(const edge_descriptor e: edgesToBeRemoved) {
        if(not std::ranges::contains(edgesToBeRemoved, assemblyGraph[e].eRc)) {
            cout << assemblyGraph[e].id << " will be removed but " <<
                assemblyGraph[assemblyGraph[e].eRc].id << " will not." << endl;
        }
        SHASTA2_ASSERT(std::ranges::contains(edgesToBeRemoved, assemblyGraph[e].eRc));
    }
#endif


    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, assemblyGraph);
    }
    removeIsolatedVertices();

}



// This cleans up linear chains by removing edges that have low
// corrected Jaccard similarity with nearby edges and
// replacing with new edges, constructed by connecting
// the remaining edges.
void AssemblyGraph::cleanupLinearChains()
{
    // Find linear chains of 3 or more edges.
    vector< vector<edge_descriptor> > chains;
    findLinearChains(*this, 3, chains);

    // Process each one of them separately.
    for(const auto& chain: chains) {
        cleanupLinearChain(chain);
    }
}



void AssemblyGraph::cleanupLinearChain(const vector<edge_descriptor>& chain)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minCommonCount = 3;
    const double correctedJaccardThresholdStart = 0.8;
    const double correctedJaccardThresholdDelta = 0.05;

    const bool debug = false;
    AssemblyGraph& assemblyGraph = *this;
    ostream html(0);
    const uint32_t representativeRegionStepCount =  uint32_t(options.representativeRegionStepCount);

    if(debug) {
        cout << "AssemblyGraph::cleanupLinearChain begins for chain:" << endl;
        for(const edge_descriptor e: chain) {
            cout << assemblyGraph[e].id << " ";
        }
        cout << endl;
    }

    // If the chain does not have at least 3 edges, do nothing.
    if(chain.size() < 3) {
        return;
    }

    // Sanity check: all vertices internal to the chain must have
    // in_degree 1 and out_degree 1.
    for(uint64_t i=1; i<chain.size(); i++) {
        const edge_descriptor e = chain[i];
        const vertex_descriptor v = source(e, assemblyGraph);
        SHASTA2_ASSERT(in_degree(v, assemblyGraph) == 1);
        SHASTA2_ASSERT(out_degree(v, assemblyGraph) == 1);
    }



    // Class to store information about a pair of edges in the chain.
    class EdgePairInfo {
    public:
        uint64_t i0;
        uint64_t i1;
        uint64_t commonCount;
        double correctedJaccard;
    };



    // Check all pairs of edges in the chain.
    vector<EdgePairInfo> edgePairInfos;
    for(uint64_t i0=0; i0<chain.size()-1; i0++) {
        const edge_descriptor e0 = chain[i0];
        for(uint64_t i1=i0+1; i1<chain.size(); i1++) {
            const edge_descriptor e1 = chain[i1];
            if(canConnect(e0, e1)) {
                const SegmentPairInformation segmentPairInformation =  SegmentStepSupport::analyzeSegmentPair(
                    html, assemblyGraph, e0, e1, representativeRegionStepCount);

                if(segmentPairInformation.commonCount < minCommonCount) {
                    continue;
                }

                if(false) {
                    cout << assemblyGraph[e0].id << " " << assemblyGraph[e1].id << " " <<
                        segmentPairInformation.commonCount << " " << segmentPairInformation.correctedJaccard << endl;
                }

                EdgePairInfo& edgePairInfo = edgePairInfos.emplace_back();
                edgePairInfo.i0 = i0;
                edgePairInfo.i1 = i1;
                edgePairInfo.commonCount = segmentPairInformation.commonCount;
                edgePairInfo.correctedJaccard = segmentPairInformation.correctedJaccard;
            }
        }
    }



    class Vertex{
    public:
        edge_descriptor e;
        // The length of the longest path ending here.
        uint64_t pathLength = invalid<uint64_t>;
        bool isOnLongestPath = false;
    };
    class Edge {
    public:
        uint64_t commonCount = 0;
        double correctedJaccard = 0.;
    };
    using Graph = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::bidirectionalS,
        Vertex,
        Edge>;
    Graph graph(chain.size());



    // Main loop over decreasing values of correctedJaccardThreshold.
    double correctedJaccardThreshold = correctedJaccardThresholdStart;
    while(true) {

        // Add edges not already present.
        for(const EdgePairInfo& edgePairInfo: edgePairInfos) {
            if(edgePairInfo.correctedJaccard < correctedJaccardThreshold) {
                continue;
            }

            auto[ignore, edgeExists] = boost::edge(edgePairInfo.i0, edgePairInfo.i1, graph);
            if(edgeExists) {
                continue;
            }

            auto[e, edgeWasAdded] = boost::add_edge(edgePairInfo.i0, edgePairInfo.i1, graph);
            Edge& edge = graph[e];
            edge.commonCount = edgePairInfo.commonCount;
            edge.correctedJaccard = edgePairInfo.correctedJaccard;
        }

        // Check reachability.
        if(isReachable(graph, 0, chain.size()-1, 0)) {
            break;
        }

        // If no reachability at any threshold, do nothing.
        if(correctedJaccardThreshold < 0.) {
            if(debug) {
                cout << "No reachability, giving up on this chain." << endl;
            }
            return;
        }

        correctedJaccardThreshold -= correctedJaccardThresholdDelta;
    }
    if(debug) {
        cout << "Reachability achieved at corrected Jaccard threshold " << correctedJaccardThreshold << endl;
    }



    // We want to find the longest path between vFirst and vLast.
    // To do that we first compute the length of the longest path
    // starting at vFirst and ending at each vertex.
    graph[0].pathLength = 0;
    for(Graph::vertex_descriptor v0=1; v0<chain.size(); v0++) {
        uint64_t maximumLength = 0;
        bool hasReachableParents = false;
        BGL_FORALL_INEDGES_T(v0, e, graph, Graph) {
            const Graph::vertex_descriptor v1 = source(e, graph);
            const uint64_t pathLength1 = graph[v1].pathLength;
            if(pathLength1 != invalid<uint64_t>) {
                maximumLength = max(maximumLength, pathLength1);
                hasReachableParents = true;
            }
        }
        if(hasReachableParents) {
            graph[v0].pathLength = maximumLength + 1;
        }
    }
    SHASTA2_ASSERT(graph[chain.size() - 1].pathLength != invalid<uint64_t>);


    // Now to compute the longest path between vFirst and vLast we start
    // at vLast and move backward, choosing at each step the parent with the
    // greatest pathLength.
    vector<Graph::vertex_descriptor> longestPath;
    longestPath.push_back(chain.size() - 1);
    while(true) {
        const Graph::vertex_descriptor v0 = longestPath.back();
        uint64_t maximumLength = 0;
        Graph::vertex_descriptor v1Best = Graph::null_vertex();
        BGL_FORALL_INEDGES(v0, e, graph, Graph) {
            const Graph::vertex_descriptor v1 = source(e, graph);
            const uint64_t pathLength1 = graph[v1].pathLength;
            if(pathLength1 != invalid<uint64_t>) {
                maximumLength = max(maximumLength, pathLength1);
                v1Best = v1;
            }
        }
        SHASTA2_ASSERT(v1Best != Graph::null_vertex());
        longestPath.push_back(v1Best);
        if(v1Best == 0) {
            break;
        }
    }
    std::ranges::reverse(longestPath);
    SHASTA2_ASSERT(longestPath.front() == 0);
    SHASTA2_ASSERT(longestPath.back() == chain.size() - 1);
    for(const Graph::vertex_descriptor v: longestPath) {
        graph[v].isOnLongestPath = true;
    }



    if(debug) {
        ofstream dot("CleanupLinearChain.dot");
        dot << std::fixed << std::setprecision(2);
        dot << "digraph cleanupLinearChain {\n";
        BGL_FORALL_VERTICES(v, graph, Graph) {

            dot << v << " [label=\"" << v << "\\n" <<
                assemblyGraph[chain[v]].id;
            if(graph[v].pathLength != invalid<uint64_t>) {
                dot << "\\n" <<
                graph[v].pathLength;
            }
            dot << "\"";
            if(graph[v].isOnLongestPath) {
                dot << " style=filled fillcolor=pink";
            }
            dot << "];\n";
        }
        BGL_FORALL_EDGES(e, graph, Graph) {
            const Edge& edge = graph[e];
            const Graph::vertex_descriptor v0 = source(e, graph);
            const Graph::vertex_descriptor v1 = target(e, graph);
            dot << v0 << "->" << v1 <<
                "[label=\"" << edge.commonCount << "/" << edge.correctedJaccard << "\"];\n";
        }
        dot << "}\n";
    }



    // At places where the longest path skips AssemblyGraph edges,
    // we have to remove the AssemblyGraph edges that are skipped
    // and replace them with a new edge to connect the skip.
    for(uint64_t i1=1; i1<longestPath.size(); i1++) {
        const uint64_t i0 = i1 - 1;
        const Graph::vertex_descriptor v0 = longestPath[i0];
        const Graph::vertex_descriptor v1 = longestPath[i1];
        if(v1 != v0+1) {
            const edge_descriptor e0 = chain[v0];
            const edge_descriptor e1 = chain[v1];
            if(debug) {
                cout << "Skipping between " <<
                    assemblyGraph[e0].id << " and " <<
                    assemblyGraph[e1].id << ":";
                for(Graph::vertex_descriptor u=v0+1; u<v1; u++) {
                    cout << " " << assemblyGraph[chain[u]].id;
                }
                cout << endl;
            }

            // Create a new edge to connect e0 to e1.
            const vertex_descriptor v0 = target(e0, assemblyGraph);
            const vertex_descriptor v1 = source(e1, assemblyGraph);

            const AnchorId anchorId0 = assemblyGraph[v0].anchorId;
            const AnchorId anchorId1 = assemblyGraph[v1].anchorId;

            // Create the new edge.
            // If the two anchors are the same, leave it empty without any steps.
            // Otherwise use the same process in Tangle1::addConnectPair.
            edge_descriptor eNew;
            tie(eNew, ignore) = add_edge(v0, v1, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
            AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
            if(anchorId0 != anchorId1) {

                // Create the RestrictedAnchorGraph, then:
                // - Remove vertices not accessible from anchorId0 and anchorId1.
                // - Remove cycles.
                // - Find the longest path.
                // - Add one step for each edge of the longest path of the RestrictedAnchorGraph.

                ostream html(0);
                const TangleMatrix tangleMatrix(
                    assemblyGraph,
                    vector<edge_descriptor>(1, e0),
                    vector<edge_descriptor>(1, e1),
                    html);

                RestrictedAnchorGraph restrictedAnchorGraph(anchors, journeys, tangleMatrix, 0, 0, html);
                vector<RestrictedAnchorGraph::edge_descriptor> longestPath;
                // restrictedAnchorGraph.findLongestPath(longestPath);
                restrictedAnchorGraph.findOptimalPath(anchorId0, anchorId1, longestPath);

                for(const RestrictedAnchorGraph::edge_descriptor re: longestPath) {
                    const auto& rEdge = restrictedAnchorGraph[re];
                    if(rEdge.anchorPair.size() == 0) {
                        newEdge.clear();
                        SHASTA2_ASSERT(0);
                    }
                    newEdge.push_back(AssemblyGraphEdgeStep(rEdge.anchorPair, rEdge.offset));
                }
            }
        }
        // Now remove the edges we skipped.
        for(Graph::vertex_descriptor u=v0+1; u<v1; u++) {
            const edge_descriptor e = chain[u];
            boost::remove_edge(e, assemblyGraph);
        }
    }


}



void AssemblyGraph::setAnnotation(edge_descriptor e, const string& annotation)
{
    AssemblyGraph& assemblyGraph = *this;
    assemblyGraph[e].annotation = annotation;
}



void AssemblyGraph::readFollowing()
{
    ReadFollowing4::ReadFollower readFollower(*this);
    // readFollower.updateAssemblyGraph(*this);
    readFollower.updateAssemblyGraphStrandSymmetric(*this);
}



void AssemblyGraph::connectDanglingSegments()
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minLength = 50000;
    const uint64_t minReadCount = 1;

    AssemblyGraph& assemblyGraph = *this;
    using Segment = edge_descriptor;
    const bool debug = false;


    // A graph to store information about dangling segments.
    class Vertex {
    public:
        Segment segment;
        bool isSource;
        bool isTarget;
    };
    using Graph = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        Vertex>;
    Graph graph;



    // To generate Graph vertices, loop over AssemblyGraph edges (segments).
    std::map<Segment, Graph::vertex_descriptor> vertexMap;
    BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {

        // If not dangling, skip it.
        const vertex_descriptor v0 = source(segment, assemblyGraph);
        const vertex_descriptor v1 = target(segment, assemblyGraph);
        const bool isSource = (out_degree(v1, assemblyGraph) == 0);
        const bool isTarget = (in_degree(v0, assemblyGraph) == 0);
        if(not (isSource or isTarget)) {
            continue;
        }

        // If too short, skip it.
        if(assemblyGraph[segment].length() < minLength) {
            continue;
        }

        // Create a vertex.
        const Graph::vertex_descriptor v = add_vertex(Vertex({segment, isSource, isTarget}), graph);
        vertexMap.insert(make_pair(segment, v));
        if(false) {
            cout << "Dangling segment " << assemblyGraph[segment].id;
            if(isSource) {
                cout << " source";
            }
            if(isTarget) {
                cout << " target";
            }
            cout << endl;
        }

    }



    // To generate edges, loop over pairs of vertices in which the first is
    // a source and the second is a target.
    BGL_FORALL_VERTICES(v0, graph, Graph) {
        if(not graph[v0].isSource) {
            continue;
        }
        BGL_FORALL_VERTICES(v1, graph, Graph) {
            if(not graph[v1].isTarget) {
                continue;
            }

            const Segment e0 = graph[v0].segment;
            const Segment e1 = graph[v1].segment;

            // Skip if the segments are the same
            // or the reverse complement of each other.
            if(e0 == e1) {
                continue;
            }
            if(assemblyGraph[e0].eRc == e1) {
                continue;
            }



            // Create the tangle matrix.
            std::ostream html(0);
            const TangleMatrix tangleMatrix(
                assemblyGraph,
                vector<edge_descriptor>(1, e0),
                vector<edge_descriptor>(1, e1),
                html);
            const double m = tangleMatrix.tangleMatrix[0][0];

            if(debug and (m > 0.)) {
                cout << "Dangling segments pair " <<
                    assemblyGraph[e0].id << " " << assemblyGraph[e1].id <<
                    ", " << m << " common oriented reads." << endl;
            }

            if(m >= double(minReadCount)) {
                add_edge(v0, v1, graph);
                if(debug) {
                    cout << "Dangling segments pair " <<
                        assemblyGraph[e0].id << " " << assemblyGraph[e1].id <<
                        " connected, " << m << " common reads." << endl;
                }
            }
        }
    }



    if(debug) {
        ofstream dot("ConnectDanglingSegments.dot");
        dot << "digraph ConnectDanglingSegments {\n";

        BGL_FORALL_VERTICES(v, graph, Graph) {
            const Vertex& vertex = graph[v];
            const Segment e = vertex.segment;
            const Segment eRc = assemblyGraph[e].eRc;
            SHASTA2_ASSERT(eRc != assemblyGraphNullEdge);
            const uint64_t id = assemblyGraph[e].id;
            const uint64_t idRc = assemblyGraph[eRc].id;
            dot << id << "[label=\"" <<
                id << "\\n(" << idRc << ")\"";
            string color = "Black";
            if(vertex.isSource) {
                if(vertex.isTarget) {
                    color = "Blue";
                } else {
                    color = "Green";
                }
            } if(vertex.isTarget) {
                if(vertex.isSource) {
                    color = "Blue";
                } else {
                    color = "Red";
                }
            }
            dot << " style=filled fillcolor=" << color << "];\n";
        }

        BGL_FORALL_EDGES(e, graph, Graph){
            const Graph::vertex_descriptor v0 = source(e, graph);
            const Graph::vertex_descriptor v1 = target(e, graph);
            const Vertex& vertex0 = graph[v0];
            const Vertex& vertex1 = graph[v1];
            dot << assemblyGraph[vertex0.segment].id <<
                "->" << assemblyGraph[vertex1.segment].id << endl;
        }
        dot << "}\n";
    }



    // Check that the graph is strand-symmetric.
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Segment segment = graph[v].segment;
        const Segment segmentRc = assemblyGraph[segment].eRc;
        SHASTA2_ASSERT(segmentRc != assemblyGraphNullEdge);
        SHASTA2_ASSERT(segmentRc != segment);
        const Graph::vertex_descriptor vRc = vertexMap.at(segmentRc);
        const Vertex& vertex = graph[v];
        const Vertex& vertexRc = graph[vRc];
        SHASTA2_ASSERT(vertex.isSource == vertexRc.isTarget);
        SHASTA2_ASSERT(vertex.isTarget == vertexRc.isSource);
    }
    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const Segment segment0 = graph[v0].segment;
        const Segment segment1 = graph[v1].segment;
        const Segment segment0Rc = assemblyGraph[segment0].eRc;
        const Segment segment1Rc = assemblyGraph[segment1].eRc;
        SHASTA2_ASSERT(segment0Rc != assemblyGraphNullEdge);
        SHASTA2_ASSERT(segment1Rc != assemblyGraphNullEdge);
        const Graph::vertex_descriptor v0Rc = vertexMap.at(segment0Rc);
        const Graph::vertex_descriptor v1Rc = vertexMap.at(segment1Rc);
        auto[eRc, edgeExists] = boost::edge(v1Rc, v0Rc, graph);
        SHASTA2_ASSERT(edgeExists);
        SHASTA2_ASSERT(eRc != e);
    }



    // Now we loop over all edges v0->v1 such that out_degree(v0)==1
    // and in_degree(v1)==1 and connected the segments corresponding to v0 and v1.
    // To keep the AssemblyGraph strand symmetric, when we make a connection
    // we also make the reverse complement connection.
    // We keep track of connections already made.
    std::set< pair<Segment, Segment> > connectionsMade;
    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        if((out_degree(v0, graph) == 1) and (in_degree(v1, graph) == 1)) {
            const edge_descriptor e0 = graph[v0].segment;
            const edge_descriptor e1 = graph[v1].segment;
            if(connectionsMade.contains({e0, e1})) {
                continue;
            }
            if(canConnect(e0, e1)) {
                if(debug) {
                    cout << "Connecting dangling segments " <<
                        assemblyGraph[e0].id << " " <<
                        assemblyGraph[e1].id << endl;
                }

                const AssemblyGraph::vertex_descriptor u0 = target(e0, assemblyGraph);
                const AssemblyGraph::vertex_descriptor u1 = source(e1, assemblyGraph);

                const AnchorId anchorId0 = assemblyGraph[u0].anchorId;
                const AnchorId anchorId1 = assemblyGraph[u1].anchorId;


                // Create the new edge.
                // If the two anchors are the same, leave it empty without any steps.
                // Otherwise use the same process in Tangle1::addConnectPair.
                edge_descriptor newSegment;
                tie(newSegment, ignore) = add_edge(u0, u1, AssemblyGraphEdge(assemblyGraph.nextEdgeId++), assemblyGraph);
                AssemblyGraphEdge& newEdge = assemblyGraph[newSegment];
                if(anchorId0 != anchorId1) {

                    // Create the RestrictedAnchorGraph, then:
                    // - Remove vertices not accessible from anchorId0 and anchorId1.
                    // - Remove cycles.
                    // - Find the longest path.
                    // - Add one step for each edge of the longest path of the RestrictedAnchorGraph.

                    ostream html(0);
                    const TangleMatrix tangleMatrix(
                        assemblyGraph,
                        vector<edge_descriptor>(1, e0),
                        vector<edge_descriptor>(1, e1),
                        html);

                    try {
                        RestrictedAnchorGraph restrictedAnchorGraph(assemblyGraph.anchors, assemblyGraph.journeys, tangleMatrix, 0, 0, html);
                        vector<RestrictedAnchorGraph::edge_descriptor> longestPath;
                        // restrictedAnchorGraph.findLongestPath(longestPath);
                        restrictedAnchorGraph.findOptimalPath(anchorId0, anchorId1, longestPath);

                        for(const RestrictedAnchorGraph::edge_descriptor re: longestPath) {
                            const auto& rEdge = restrictedAnchorGraph[re];
                            if(rEdge.anchorPair.size() == 0) {
                                newEdge.clear();
                                SHASTA2_ASSERT(0);
                            }
                            newEdge.push_back(AssemblyGraphEdgeStep(rEdge.anchorPair, rEdge.offset));

                        }
                    } catch(RestrictedAnchorGraph::NoTransitions&) {
                        cout << "Could not connect " << assemblyGraph[e0].id <<
                            " with " << assemblyGraph[e1].id << endl;
                        SHASTA2_ASSERT(0);
                    }

                }
                connectionsMade.insert({e0, e1});
                const Segment newSegmentRc = createReverseComplementEdge(newSegment);
                connectionsMade.insert({assemblyGraph[e1].eRc, assemblyGraph[e0].eRc});

                if(debug) {
                    cout << "Created " << assemblyGraph[newSegment].id <<
                        " connecting " << assemblyGraph[e0].id << " with " <<
                        assemblyGraph[e1].id << endl;
                    cout << "Reverse complement: created " << assemblyGraph[newSegmentRc].id <<
                        " connecting " << assemblyGraph[assemblyGraph[e1].eRc].id << " with " <<
                        assemblyGraph[assemblyGraph[e0].eRc].id << endl;

                }

            }
        }
    }

    strandSymmetricCompress();
}



// Connect two segments (edges) using the RestrictedAnchorGraph.
// Return the newly created edge.
AssemblyGraph::edge_descriptor AssemblyGraph::connect(edge_descriptor e0, edge_descriptor e1, bool deep)
{
    AssemblyGraph& assemblyGraph = *this;

    const AssemblyGraph::vertex_descriptor v0 = target(e0, assemblyGraph);
    const AssemblyGraph::vertex_descriptor v1 = source(e1, assemblyGraph);

    const AnchorId anchorId0 = assemblyGraph[v0].anchorId;
    const AnchorId anchorId1 = assemblyGraph[v1].anchorId;


    // Create the new edge.
    // If the two anchors are the same, leave it empty without any steps.
    // Otherwise use the RestrictedAnchorGraph.
    edge_descriptor newSegment;
    tie(newSegment, ignore) = add_edge(v0, v1, AssemblyGraphEdge(assemblyGraph.nextEdgeId++), assemblyGraph);
    AssemblyGraphEdge& newEdge = assemblyGraph[newSegment];
    if(anchorId0 != anchorId1) {

        // Create the RestrictedAnchorGraph, then:
        // - Remove vertices not accessible from anchorId0 and anchorId1.
        // - Remove cycles.
        // - Find the longest path.
        // - Add one step for each edge of the longest path of the RestrictedAnchorGraph.

        ostream html(0);
        const TangleMatrix tangleMatrix(
            assemblyGraph,
            vector<edge_descriptor>(1, e0),
            vector<edge_descriptor>(1, e1),
            html);

        try {
            RestrictedAnchorGraph restrictedAnchorGraph(assemblyGraph.anchors, assemblyGraph.journeys, tangleMatrix, 0, 0, html, deep);
            vector<RestrictedAnchorGraph::edge_descriptor> longestPath;
            // restrictedAnchorGraph.findLongestPath(longestPath);
            restrictedAnchorGraph.findOptimalPath(anchorId0, anchorId1, longestPath);

            for(const RestrictedAnchorGraph::edge_descriptor re: longestPath) {
                const auto& rEdge = restrictedAnchorGraph[re];
                if(rEdge.anchorPair.size() == 0) {
                    newEdge.clear();
                    SHASTA2_ASSERT(0);
                }
                newEdge.push_back(AssemblyGraphEdgeStep(rEdge.anchorPair, rEdge.offset));

            }
        } catch(std::exception& e) {
            cout << "Could not connect " << assemblyGraph[e0].id <<
                " with " << assemblyGraph[e1].id << endl;
            cout << e.what() << endl;
            assemblyGraph.write("AtAssertion");
            SHASTA2_ASSERT(0);
        }

    }

    return newSegment;
}



// Simple connection of two segments (edges) without using
// the RestrictedAnchorGraph.
void AssemblyGraph::simpleConnect(edge_descriptor e0, edge_descriptor e1)
{
    AssemblyGraph& assemblyGraph = *this;

    // Access the target vertex of e0 and the source vertex of e1.
    const vertex_descriptor v0 = target(e0, assemblyGraph);
    const vertex_descriptor v1 = source(e1, assemblyGraph);

    // If the two vertices are the same, e0 and e1 are already connected,
    // and we don't have to do anything.
    if(v0 == v1) {
        return;
    }

    // Get the corresponding AnchorIds.
    const AnchorId anchorId0 = assemblyGraph[v0].anchorId;
    const AnchorId anchorId1 = assemblyGraph[v1].anchorId;

    // If the anchorIds are the same, connect v0 and v1 using an empty
    // AssemblyGraphEdge.
    if(anchorId0 == anchorId1) {
        add_edge(v0, v1, AssemblyGraphEdge(assemblyGraph.nextEdgeId++), assemblyGraph);
        return;
    }

    // Otherwise, create a new AssemblyGraphEdge consisting of a single
    // AssemblyGraphEdgeStep, using the common reads between the last
    // step of e0 and the first step of e1.
    vector<OrientedReadId> orientedReadIds;
    findOrientedReadIdsForSimpleConnect(e0, e1, orientedReadIds);
    SHASTA2_ASSERT(not orientedReadIds.empty());
    const auto [eNew, ignore] = add_edge(v0, v1, AssemblyGraphEdge(assemblyGraph.nextEdgeId++), assemblyGraph);
    AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];
    AnchorPair anchorPair;
    anchorPair.anchorIdA = anchorId0;
    anchorPair.anchorIdB = anchorId1;
    anchorPair.orientedReadIds = orientedReadIds;
    const uint32_t offset = anchorPair.getAverageOffset(anchors);
    edgeNew.push_back(AssemblyGraphEdgeStep(anchorPair, offset));

}



bool AssemblyGraph::canSimpleConnect(edge_descriptor e0, edge_descriptor e1)
{
    vector<OrientedReadId> orientedReadIds;
    findOrientedReadIdsForSimpleConnect(e0, e1, orientedReadIds);
    return not orientedReadIds.empty();
}



void AssemblyGraph::findOrientedReadIdsForSimpleConnect(
    edge_descriptor e0,
    edge_descriptor e1,
    vector<OrientedReadId>& orientedReadIds) const
{
    const AssemblyGraph& assemblyGraph = *this;
    orientedReadIds.clear();

    // Access the two edges.
    const AssemblyGraphEdge& edge0 = assemblyGraph[e0];
    const AssemblyGraphEdge& edge1 = assemblyGraph[e1];

    if(edge0.empty() or edge1.empty()) {
        return;
    }

    // Access the last step of e0 and the first step of e1.
    const AssemblyGraphEdgeStep& step0 = edge0.back();
    const AssemblyGraphEdgeStep& step1 = edge1.front();

    const AnchorPair& anchorPair0 = step0.anchorPair;
    const AnchorPair& anchorPair1 = step1.anchorPair;

    // Find common oriented reads between step0 and step1.
    std::ranges::set_intersection(
        anchorPair0.orientedReadIds,
        anchorPair1.orientedReadIds,
        back_inserter(orientedReadIds));
}



// Clear reverse complement information from all vertices and edges.
// This needs to be done before operations that don't maintain
// the vertices/edges vRc/eRc field.
void AssemblyGraph::clearReverseComplementInformation()
{
    AssemblyGraph& assemblyGraph = *this;

    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        assemblyGraph[v].vRc = null_vertex();
    }

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        assemblyGraph[e].eRc = assemblyGraphNullEdge;
    }
}



// Remove zero length segments by collapsing their vertices.
void AssemblyGraph::removeZeroLengthSegments()
{
    AssemblyGraph& assemblyGraph = *this;

    // Map vertices to integers.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    vector<vertex_descriptor> vertexTable;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexTable.push_back(v);
        vertexIndexMap.insert(make_pair(v, vertexIndex++));
    }

    // Compute connected components using only zero-length edges,
    // that is, edges with zero steps.
    // Each non-trivial connected component will be collapsed
    // into a single vertex.
    DisjointSets disjointSets(vertexIndexMap.size());
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        if(assemblyGraph[e].empty()) {
            const vertex_descriptor v0 = source(e, assemblyGraph);
            const vertex_descriptor v1 = target(e, assemblyGraph);
            const uint64_t vertexIndex0 = vertexIndexMap.at(v0);
            const uint64_t vertexIndex1 = vertexIndexMap.at(v1);
            disjointSets.unionSet(vertexIndex0, vertexIndex1);
        }
    }

    // Get the non-trivial connected components.
    vector< vector<uint64_t> > components;
    disjointSets.gatherComponents(2, components);
    // cout << "Found " << components.size() << " groups of vertices to be collapsed." << endl;

    // Collapse the vertices in each component.
    for(const vector<uint64_t>& component: components) {
        vector<vertex_descriptor> componentVertices;
        for(const uint64_t vertexIndex: component) {
            componentVertices.push_back(vertexTable[vertexIndex]);
        }
        collapseVertices(componentVertices);
    }
}



uint64_t AssemblyGraph::countZeroLengthSegments() const
{
    const AssemblyGraph& assemblyGraph = *this;

    uint64_t n = 0;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        if(assemblyGraph[e].empty()) {
            ++n;
        }
    }
    return n;
}



void AssemblyGraph::removeZeroLengthSegmentsStrandSymmetric()
{
    performanceLog << timestamp << "AssemblyGraph::removeZeroLengthSegmentsStrandSymmetric begins." << endl;
    const bool debug = false;

    AssemblyGraph& assemblyGraph = *this;

    // Map vertices to integers.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    vector<vertex_descriptor> vertexTable;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexTable.push_back(v);
        vertexIndexMap.insert(make_pair(v, vertexIndex++));
    }
    const uint64_t n = vertexTable.size();
    SHASTA2_ASSERT(vertexIndexMap.size() == n);


    // Compute connected components using only zero-length edges,
    // that is, edges with zero steps.
    // Each non-trivial connected component will be collapsed
    // into a single vertex.
    // Only get the connected components corresponding to positive anchors.
    // The ones with negative anchors will be processed below as the
    // reverse complement component.
    DisjointSets disjointSets(n);
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        if(assemblyGraph[e].empty()) {
            const vertex_descriptor v0 = source(e, assemblyGraph);
            const vertex_descriptor v1 = target(e, assemblyGraph);
            const AnchorId anchorId0 = assemblyGraph[v0].anchorId;
            const AnchorId anchorId1 = assemblyGraph[v1].anchorId;
            SHASTA2_ASSERT(anchorId0 == anchorId1);
            if(anchorId0 & 1) {
                continue;
            }
            const uint64_t vertexIndex0 = vertexIndexMap.at(v0);
            const uint64_t vertexIndex1 = vertexIndexMap.at(v1);
            SHASTA2_ASSERT(vertexIndex0 < n);
            SHASTA2_ASSERT(vertexIndex1 < n);
            disjointSets.unionSet(vertexIndex0, vertexIndex1);

            if(debug) {
                cout << "Empty edge " << assemblyGraph[e].id <<
                    ", source " << assemblyGraph[v0].id << "/" << vertexIndex0 << "/" << anchorIdToString(assemblyGraph[v0].anchorId) <<
                    ", target " << assemblyGraph[v1].id << "/" << vertexIndex1 << "/" << anchorIdToString(assemblyGraph[v0].anchorId) << endl;
            }
        }
    }

    // Get the non-trivial connected components.
    vector< vector<uint64_t> > components;
    disjointSets.gatherComponents(2, components);
    // cout << "Found " << components.size() << " groups of vertices to be collapsed." << endl;



    // Collapse the vertices in each component and its reverse complement.
    for(const vector<uint64_t>& component: components) {
        vector<vertex_descriptor> componentVertices;
        for(const uint64_t vertexIndex: component) {
            SHASTA2_ASSERT(vertexIndex < n);
            componentVertices.push_back(vertexTable[vertexIndex]);
        }

        if(debug) {
            cout << "Indexes of vertices to be collapsed:" << endl;
            for(const uint64_t vertexIndex: component) {
                cout << vertexIndex << endl;
            }
            cout << "Collapsing vertices:" << endl;;
            for(const vertex_descriptor v: componentVertices) {
                const AssemblyGraphVertex& vertex = assemblyGraph[v];
                cout << v << " " << vertex.id << " " << anchorIdToString(assemblyGraph[v].anchorId) << endl;
            }
        }

        // Sanity check: they must all have the same AnchorId.
        AnchorId anchorId = invalid<AnchorId>;
        for(const vertex_descriptor v: componentVertices) {
            if(anchorId == invalid<AnchorId>) {
                anchorId = assemblyGraph[v].anchorId;
            } else {
                SHASTA2_ASSERT(anchorId == assemblyGraph[v].anchorId);
            }
        }
        if(debug) {
            cout << "AnchorId " << anchorIdToString(anchorId) << " " << (anchorId & 1UL) << endl;
        }
        SHASTA2_ASSERT((anchorId & 1UL) == 0);

        collapseVerticesStrandSymmetric(componentVertices);
    }
    performanceLog << timestamp << "AssemblyGraph::removeZeroLengthSegmentsStrandSymmetric ends." << endl;
}



void AssemblyGraph::collapseVertices(const vector<vertex_descriptor>& verticesToBeCollapsed)
{
    SHASTA2_ASSERT(verticesToBeCollapsed.size() > 1);
    const bool debug = false;

    AssemblyGraph& assemblyGraph = *this;

    if(debug) {
        cout << "This group of vertices will be collapsed:";
        for(const vertex_descriptor v: verticesToBeCollapsed) {
            cout << " " << assemblyGraph[v].id;
        }
        cout << endl;
    }

    // Sanity check: they must all kave the same AnchorId.
    AnchorId anchorId = invalid<AnchorId>;
    for(const vertex_descriptor v: verticesToBeCollapsed) {
        if(anchorId == invalid<AnchorId>) {
            anchorId = assemblyGraph[v].anchorId;
        } else {
            SHASTA2_ASSERT(anchorId == assemblyGraph[v].anchorId);
        }
    }

    // Remove edges between these vertices.
    vector<edge_descriptor> edgesToBeRemoved;
    for(const vertex_descriptor v0: verticesToBeCollapsed) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if(std::ranges::find(verticesToBeCollapsed, v1) != verticesToBeCollapsed.end()) {
                edgesToBeRemoved.push_back(e);
            }
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        if(debug) {
            cout << "Removing edge " << assemblyGraph[e].id << endl;
        }
        boost::remove_edge(e, assemblyGraph);
    }

    // Create the new vertex.
    const vertex_descriptor vNew = add_vertex(AssemblyGraphVertex(anchorId, nextVertexId++), assemblyGraph);
    if(debug) {
        cout << "Created new vertex " << assemblyGraph[vNew].id << endl;
    }



    // Reroute incoming/outgoing edges to the new, collapsed vertex.
    edgesToBeRemoved.clear();
    for(const vertex_descriptor v0: verticesToBeCollapsed) {

        if(debug) {
            cout << "Rerouting out-edges of " << assemblyGraph[v0].id << endl;
        }
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            AssemblyGraphEdge& edge = assemblyGraph[e];
            const vertex_descriptor v1 = target(e, assemblyGraph);
            SHASTA2_ASSERT(std::ranges::find(verticesToBeCollapsed, v1) == verticesToBeCollapsed.end());

            auto[eNew, wasAdded] = add_edge(vNew, v1, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
            SHASTA2_ASSERT(wasAdded);
            AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
            if(debug) {
                cout << "Edge " << edge.id << " rerouted to collapsed vertex, becomes " << newEdge.id << endl;
            }
            newEdge.swapSteps(edge);

            // Activate this when read following becomes strand-symmetric.
            /*
            if(edge.eRc != assemblyGraphNullEdge) {
                SHASTA2_ASSERT(assemblyGraph[edge.eRc].eRc == e);
                assemblyGraph[edge.eRc].eRc = eNew;
                newEdge.eRc = edge.eRc;
            }
            */

            edgesToBeRemoved.push_back(e);
        }

        if(debug) {
            cout << "Rerouting in-edges of " << assemblyGraph[v0].id << endl;
        }
        BGL_FORALL_INEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            AssemblyGraphEdge& edge = assemblyGraph[e];
            const vertex_descriptor v1 = source(e, assemblyGraph);
            SHASTA2_ASSERT(std::ranges::find(verticesToBeCollapsed, v1) == verticesToBeCollapsed.end());

            auto[eNew, wasAdded] = add_edge(v1, vNew, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
            SHASTA2_ASSERT(wasAdded);
            AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
            if(debug) {
                cout << "Edge " << edge.id << " rerouted to collapsed vertex, becomes " << newEdge.id << endl;
            }
            newEdge.swapSteps(edge);

            // Activate this when read following becomes strand-symmetric.
            /*
            if(edge.eRc != assemblyGraphNullEdge) {
                SHASTA2_ASSERT(assemblyGraph[edge.eRc].eRc == e);
                assemblyGraph[edge.eRc].eRc = eNew;
                newEdge.eRc = edge.eRc;
            }
            */

            edgesToBeRemoved.push_back(e);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        if(debug) {
            cout << "Removing edge " << assemblyGraph[e].id << endl;
        }
        boost::remove_edge(e, assemblyGraph);
    }

    // Finally, we can remove the vertices that were collapsed.
    for(const vertex_descriptor v: verticesToBeCollapsed) {
        SHASTA2_ASSERT(in_degree(v, assemblyGraph) == 0);
        SHASTA2_ASSERT(out_degree(v, assemblyGraph) == 0);
        boost::remove_vertex(v, assemblyGraph);
    }
}



void AssemblyGraph::collapseVerticesStrandSymmetric(
    const vector<vertex_descriptor>& verticesToBeCollapsed)
{
    SHASTA2_ASSERT(verticesToBeCollapsed.size() > 1);
    const bool debug = false;

    AssemblyGraph& assemblyGraph = *this;

    // Sanity check: they must all kave the same AnchorId.
    AnchorId anchorId = invalid<AnchorId>;
    for(const vertex_descriptor v: verticesToBeCollapsed) {
        if(anchorId == invalid<AnchorId>) {
            anchorId = assemblyGraph[v].anchorId;
        } else {
            SHASTA2_ASSERT(anchorId == assemblyGraph[v].anchorId);
        }
    }

    // Remove edges between these vertices and their reverse complement.
    vector<edge_descriptor> edgesToBeRemoved;
    for(const vertex_descriptor v0: verticesToBeCollapsed) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if(std::ranges::find(verticesToBeCollapsed, v1) != verticesToBeCollapsed.end()) {
                edgesToBeRemoved.push_back(e);
                edgesToBeRemoved.push_back(assemblyGraph[e].eRc);
            }
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        if(debug) {
            cout << "Removing edge " << assemblyGraph[e].id << endl;
        }
        boost::remove_edge(e, assemblyGraph);
    }

    // Create the new vertex and its reverse complement.
    const vertex_descriptor vNew = add_vertex(AssemblyGraphVertex(anchorId, nextVertexId++), assemblyGraph);
    if(debug) {
        cout << "Created new vertex " << assemblyGraph[vNew].id << endl;
    }
    createReverseComplementVertex(vNew);



    // Reroute incoming/outgoing edges to the new, collapsed vertex.
    // Do the same operations for the reverse complements of these vertices and edges.
    edgesToBeRemoved.clear();
    for(const vertex_descriptor v0: verticesToBeCollapsed) {

        if(debug) {
            cout << "Rerouting out-edges of " << assemblyGraph[v0].id << endl;
        }
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            AssemblyGraphEdge& edge = assemblyGraph[e];
            const vertex_descriptor v1 = target(e, assemblyGraph);
            SHASTA2_ASSERT(std::ranges::find(verticesToBeCollapsed, v1) == verticesToBeCollapsed.end());

            auto[eNew, wasAdded] = add_edge(vNew, v1, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
            SHASTA2_ASSERT(wasAdded);
            AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
            if(debug) {
                cout << "Edge " << edge.id << " rerouted to collapsed vertex, becomes " << newEdge.id << endl;
            }
            newEdge.swapSteps(edge);
            createReverseComplementEdge(eNew);

            edgesToBeRemoved.push_back(e);
            edgesToBeRemoved.push_back(assemblyGraph[e].eRc);
        }

        if(debug) {
            cout << "Rerouting in-edges of " << assemblyGraph[v0].id << endl;
        }
        BGL_FORALL_INEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            AssemblyGraphEdge& edge = assemblyGraph[e];
            const vertex_descriptor v1 = source(e, assemblyGraph);
            SHASTA2_ASSERT(std::ranges::find(verticesToBeCollapsed, v1) == verticesToBeCollapsed.end());

            auto[eNew, wasAdded] = add_edge(v1, vNew, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
            SHASTA2_ASSERT(wasAdded);
            AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
            if(debug) {
                cout << "Edge " << edge.id << " rerouted to collapsed vertex, becomes " << newEdge.id << endl;
            }
            newEdge.swapSteps(edge);

            createReverseComplementEdge(eNew);

            edgesToBeRemoved.push_back(e);
            edgesToBeRemoved.push_back(assemblyGraph[e].eRc);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        if(debug) {
            cout << "Removing edge " << assemblyGraph[e].id << endl;
        }
        boost::remove_edge(e, assemblyGraph);
    }



    // Finally, we can remove the vertices that were collapsed and their reverse complements.
    for(const vertex_descriptor v: verticesToBeCollapsed) {
        SHASTA2_ASSERT(in_degree(v, assemblyGraph) == 0);
        SHASTA2_ASSERT(out_degree(v, assemblyGraph) == 0);
        boost::remove_vertex(assemblyGraph[v].vRc, assemblyGraph);
        boost::remove_vertex(v, assemblyGraph);
    }
}



void AssemblyGraph::clearAllSequence()
{
    AssemblyGraph& assemblyGraph = *this;

    // Loop over all edges.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {

        // Mark this edge as not assembled.
        AssemblyGraphEdge& edge = assemblyGraph[e];
        edge.wasAssembled = false;

        // Clear the sequence in all the steps of this edge.
        for(AssemblyGraphEdgeStep& step: edge) {
            step.sequence.clear();
            step.sequence.shrink_to_fit();
        }
    }
}



void AssemblyGraph::addEdgeToBeAssembledStrandSymmetric(edge_descriptor e)
{
    const AssemblyGraph& assemblyGraph = *this;

    // Get the reverse complement edge.
    const edge_descriptor eRc = assemblyGraph[e].eRc;

    // If there is no reverse complement edge,
    // we will assemble this edge.
    if(eRc == assemblyGraphNullEdge) {
        edgesToBeAssembledA.insert(e);
        return;
    }

    // If we already have the reverse complement in edgesToBeAssembledA,
    // store this edge in edgesToBeAssembledB.
    // Its sequence  will be obtained by reverse complementing
    // the sequence of its reverse complement.
    if(edgesToBeAssembledA.contains(eRc)) {
        edgesToBeAssembledB.insert(e);
    } else {
        edgesToBeAssembledA.insert(e);
    }
}



void AssemblyGraph::assembleStrandSymmetric()
{
    AssemblyGraph& assemblyGraph = *this;
    using shasta2::Base;

    cout << "Strand symmetric sequence assembly for " <<
        edgesToBeAssembledA.size() << " + " <<
        edgesToBeAssembledB.size() << " edges begins." << endl;

    // Assemble the edges stored in edgesToBeAssembledA.
    vector<edge_descriptor> edgesToBeAssembled(edgesToBeAssembledA.begin(), edgesToBeAssembledA.end());
    assemble(edgesToBeAssembled);



    // For the edges stored in edgesToBeAssembledB, obtain the sequence
    // by reverse complement.
    for(const edge_descriptor e: edgesToBeAssembledB) {
        AssemblyGraphEdge& edge = assemblyGraph[e];

        // Access the reverse complemented edge.
        const edge_descriptor eRc = edge.eRc;
        SHASTA2_ASSERT(eRc != assemblyGraphNullEdge);
        SHASTA2_ASSERT(eRc != e);
        const AssemblyGraphEdge& edgeRc = assemblyGraph[eRc];

        // Loop over all steps.
        const uint64_t n = edge.size();
        SHASTA2_ASSERT(edgeRc.size() == n);
        for(uint64_t i=0; i<n; i++) {
            AssemblyGraphEdgeStep& step = edge[i];
            const AssemblyGraphEdgeStep& stepRc = edgeRc[n - 1 - i];
            vector<Base>& sequence = step.sequence;
            const vector<Base>& sequenceRc = stepRc.sequence;
            for(const Base b: std::views::reverse(sequenceRc)) {
                sequence.push_back(b.complement());
            }
        }
        edge.wasAssembled = true;

    }

    cout << "Strand symmetric sequence assembly for " <<
        edgesToBeAssembledA.size() << " + " <<
        edgesToBeAssembledB.size() << " edges ends." << endl;
}



void AssemblyGraph::assembleAllStrandSymmetric()
{
    AssemblyGraph& assemblyGraph = *this;

    writeMemoryStatistics("AssemblyGraph::assembleAllStrandSymmetric begins");

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        addEdgeToBeAssembledStrandSymmetric(e);
    }
    assembleStrandSymmetric();

    writeMemoryStatistics("AssemblyGraph::assembleAllStrandSymmetric ends");
}



// Make the AssemblyGraph single-stranded.
// This requires to AssemblyGraph to be strand-symmetric in input.
// It removes half of the edges (segments).
// We iterate over increasing values of the maxLength parameter
// used to compute tangles.
// If we find a self-complementary tangle, we disconnect
// all of its segments.
// After doing that, if maxLength is greater than
// the longest segment length we terminate the loop.
// We then recompute the tangles (which are now connected components),
// and we are guaranteed that no
// more self-complementary tangles (connected component) are present.
// We then remove one tangle in each reverse complemented pair.
void AssemblyGraph::makeSingleStranded()
{
    AssemblyGraph& assemblyGraph = *this;
    check();

    // Find the maximum segment length.
    uint64_t maxSegmentLength = 0;
    BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {
        maxSegmentLength = max(maxSegmentLength, assemblyGraph[segment].length());
    }



    // The main loop. At each iteration we disconnect all of the
    // self-complementary tangles.
    uint64_t maxLength = 1000;
    vector< vector<vertex_descriptor> > tangles;
    vector<uint64_t> tanglesRc;
    while(true) {
        createTanglesBySegmentLength(maxLength, tangles, tanglesRc);

        // Disconnect all the segments in self-complementary tangles.
        // Make sure to do it keeping the AssemblyGraph strand symmetric.
        vector<Segment> segmentsToBeRemoved;
        for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
            if(tanglesRc[tangleId] == tangleId) {
                const Tangle tangle(assemblyGraph, tangles[tangleId]);
                for(const Segment oldSegment: tangle.tangleEdges) {
                    const Segment oldSegmentRc = assemblyGraph[oldSegment].eRc;
                    SHASTA2_ASSERT(oldSegmentRc != oldSegment);
                    if(id(oldSegment) < id(oldSegmentRc)) {
                        const Segment newSegment = createDisconnectedCopy(oldSegment);
                        createReverseComplementVertex(source(newSegment, assemblyGraph));
                        createReverseComplementVertex(target(newSegment, assemblyGraph));
                        createReverseComplementEdge(newSegment);
                        segmentsToBeRemoved.push_back(oldSegment);
                        segmentsToBeRemoved.push_back(oldSegmentRc);
                    }
                }
            }
        }
        // Only now we can remove the old Segments.
        for(const Segment segment: segmentsToBeRemoved) {
            boost::remove_edge(segment, assemblyGraph);
        }

        if(maxLength >= maxSegmentLength) {
            break;
        }
        maxLength *= 2;
    }
    check();



    // Now recompute the tangles at maxLength.
    // These are now connected components, and none of them
    // is self-complementary. Remove all the segments
    // in one tangle of each pair.
    createTanglesBySegmentLength(maxLength, tangles, tanglesRc);
    for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
        const uint64_t tangleIdRc = tanglesRc[tangleId];
        SHASTA2_ASSERT(tangleIdRc != tangleId);
        if(tangleId > tangleIdRc) {
            const Tangle tangle(assemblyGraph, tangles[tangleId]);
            for(const Segment segment: tangle.tangleEdges) {
                boost::remove_edge(segment, assemblyGraph);
            }
        }
    }

    // Now we clear the vRc and eRc fields from all vertices and edges, because
    // the AssemblyGraph is no longer strand-symmetric.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        assemblyGraph[v].vRc = null_vertex();
    }
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        assemblyGraph[e].eRc = assemblyGraphNullEdge;
    }

    // Finally, remove any vertices that are left isolated.
    removeIsolatedVertices();
}



// This creates an edge eB, identical to the reverse complement
// of edge eA. It sets the eRc fields in eA and eB
// and returns eB.
// The vertices of eB must already exist.
AssemblyGraph::edge_descriptor AssemblyGraph::createReverseComplementEdge(
    edge_descriptor eA,
    uint64_t id)
{
    AssemblyGraph& assemblyGraph = *this;
    const bool debug = false;

    // Access the eA edge and its vertices.
    AssemblyGraphEdge& edgeA = assemblyGraph[eA];
    const vertex_descriptor vA0 = source(eA, assemblyGraph);
    const vertex_descriptor vA1 = target(eA, assemblyGraph);
    const AssemblyGraphVertex& vertexA0 = assemblyGraph[vA0];
    const AssemblyGraphVertex& vertexA1 = assemblyGraph[vA1];
    // cout << "MMM " << vertexA0.id << " " << vertexA1.id << endl;

    // Create the new edge eB.
    const vertex_descriptor vB0 = vertexA1.vRc;
    const vertex_descriptor vB1 = vertexA0.vRc;
    // cout << "NNN " << assemblyGraph[vB0].id << " " << assemblyGraph[vB1].id << endl;
    if(id == invalid<uint64_t>) {
        id = nextEdgeId++;
    }
    auto[eB, wasAdded] = add_edge(vB0, vB1, AssemblyGraphEdge(id), assemblyGraph);
    // cout << "OOO" << endl;
    SHASTA2_ASSERT(wasAdded);
    AssemblyGraphEdge& edgeB = assemblyGraph[eB];

    // Update the eRc fields of both edges.
    edgeA.eRc = eB;
    edgeB.eRc = eA;

    // cout << "PPP " << edgeA.id << " " << edgeB.id << endl;
    // cout << "QQQ " << edgeA.size() << " " << edgeB.size() << endl;



    // Now add the steps of eA to eB, with reverse complement.
    for(const AssemblyGraphEdgeStep& stepA: edgeA) {
        AssemblyGraphEdgeStep& stepB = edgeB.emplace_back();
        stepB.offset = stepA.offset;
        stepB.anchorPair = stepA.anchorPair;

        // Swap and reverse complement the AnchorIds.
        std::swap(stepB.anchorPair.anchorIdA, stepB.anchorPair.anchorIdB);
        stepB.anchorPair.anchorIdA = reverseComplementAnchorId(stepB.anchorPair.anchorIdA);
        stepB.anchorPair.anchorIdB = reverseComplementAnchorId(stepB.anchorPair.anchorIdB);

        // Reverse complement the OrientedReadIds.
        for(OrientedReadId& orientedReadId: stepB.anchorPair.orientedReadIds) {
            orientedReadId.flipStrand();
        }
    }
    // Finally, reverse the steps of eB.
    std::ranges::reverse(edgeB);



    if(debug) {
        cout << "addReverseComplementEdge details:" << endl;

        cout << "edgeA: " << edgeA.id << " " <<
            anchorIdToString(assemblyGraph[source(eA, assemblyGraph)].anchorId) << " " <<
            anchorIdToString(assemblyGraph[target(eA, assemblyGraph)].anchorId) << endl;
        cout << "edgeA steps:" << endl;
        for(const auto& step: edgeA) {
            cout << anchorIdToString(step.anchorPair.anchorIdA) << " " <<
                anchorIdToString(step.anchorPair.anchorIdB) << endl;
        }

        cout << "edgeB: " << edgeB.id << " " <<
            anchorIdToString(assemblyGraph[source(eB, assemblyGraph)].anchorId) << " " <<
            anchorIdToString(assemblyGraph[target(eB, assemblyGraph)].anchorId) << endl;
        cout << "edgeB steps:" << endl;
        for(const auto& step: edgeB) {
            cout << anchorIdToString(step.anchorPair.anchorIdA) << " " <<
                anchorIdToString(step.anchorPair.anchorIdB) << endl;
        }
    }

    // Sanity check.
    if(not edgeB.empty()) {
        if(
            (edgeB.front().anchorPair.anchorIdA != assemblyGraph[vB0].anchorId) or
            (edgeB.back() .anchorPair.anchorIdB != assemblyGraph[vB1].anchorId)) {
            write("AtAssertion");
            cout << "Assertion failing with:" << endl;
            cout << "edgeB " << edgeB.id << endl;
            cout << "edgeB.front().anchorPair.anchorIdA " << edgeB.front().anchorPair.anchorIdA << endl;
            cout << "edgeB.back().anchorPair.anchorIdB "  << edgeB.back() .anchorPair.anchorIdB << endl;
            cout << " assemblyGraph[vB0].anchorId " << assemblyGraph[vB0].anchorId << endl;
            cout << " assemblyGraph[vB1].anchorId " << assemblyGraph[vB1].anchorId << endl;
            cout << "See AssemblyGraph stage AtAssertion." << endl;
        }
        SHASTA2_ASSERT(edgeB.front().anchorPair.anchorIdA == assemblyGraph[vB0].anchorId);
        SHASTA2_ASSERT(edgeB.back().anchorPair.anchorIdB == assemblyGraph[vB1].anchorId);
    }

    return eB;
}



// This creates a vertex vB, identical to the reverse complement
// of vertex vA. It sets the vRc fields in vA and vB
// and returns vB.
AssemblyGraph::vertex_descriptor AssemblyGraph::createReverseComplementVertex(vertex_descriptor vA)
{
    const bool debug = false;

    AssemblyGraph& assemblyGraph = *this;
    AssemblyGraphVertex& vertexA = assemblyGraph[vA];

    // Get the AnchorIds.
    const AnchorId anchorIdA = vertexA.anchorId;
    const AnchorId anchorIdB = reverseComplementAnchorId(anchorIdA);

    // Create the new vertex.
    const vertex_descriptor vB = add_vertex(AssemblyGraphVertex(anchorIdB, nextVertexId++), assemblyGraph);
    AssemblyGraphVertex& vertexB = assemblyGraph[vB];

    // Set the vRc fields.
    vertexA.vRc = vB;
    vertexB.vRc = vA;

    if(debug) {
        cout << "Created a reverse complemented copy of a vertex: " <<
            vertexA.id << " " << anchorIdToString(vertexA.anchorId) <<
            vertexB.id << " " << anchorIdToString(vertexB.anchorId) << endl;

    }

    return vB;
}



// This makes a copy of a Segment, disconnected at its end,
// then removes the original Segment.
// It also performs the same operation on the reverse complement of the segment,
// to keep the AssemblyGraph strand symmetric, and sets the eRc and vRc
// fields of the newly created edges and vertices.
// It returns the edge_descriptor for the newly created copy.
AssemblyGraph::edge_descriptor AssemblyGraph::disconnectAtEnd(edge_descriptor eOld)
{
    AssemblyGraph& assemblyGraph = *this;

    // Get some information we need below.
    const AssemblyGraphEdge& edgeOld = assemblyGraph[eOld];
    const edge_descriptor eOldRc = edgeOld.eRc;
    const vertex_descriptor v0 = source(eOld, assemblyGraph);
    const vertex_descriptor v1Old = target(eOld, assemblyGraph);
    const AnchorId anchorId1 = assemblyGraph[v1Old].anchorId;

    // Create the new target vertex of eOld and the new source vertex of
    // its reverse complement.
    const vertex_descriptor v1New = add_vertex(AssemblyGraphVertex(anchorId1, nextVertexId++), assemblyGraph);
    createReverseComplementVertex(v1New);

    // Create the edge for the new Segment, keeping the same id.
    edge_descriptor eNew;
    tie(eNew, ignore)= add_edge(v0, v1New, edgeOld, assemblyGraph);

    // Also create the new reverse complement segment.
    createReverseComplementEdge(eNew, id(eOldRc));

    // Remove eOld and its reverse complement.
    boost::remove_edge(eOld, assemblyGraph);
    boost::remove_edge(eOldRc, assemblyGraph);

    return eNew;
}




// Same as above, but disconnect at the beginning.
AssemblyGraph::edge_descriptor AssemblyGraph::disconnectAtBeginning(edge_descriptor eOld)
{
    AssemblyGraph& assemblyGraph = *this;

    // Get some information we need below.
    const AssemblyGraphEdge& edgeOld = assemblyGraph[eOld];
    const edge_descriptor eOldRc = edgeOld.eRc;
    const vertex_descriptor v0Old = source(eOld, assemblyGraph);
    const vertex_descriptor v1 = target(eOld, assemblyGraph);
    const AnchorId anchorId0 = assemblyGraph[v0Old].anchorId;

    // Create the new source vertex of eOld and the new target vertex of
    // its reverse complement.
    const vertex_descriptor v0New = add_vertex(AssemblyGraphVertex(anchorId0, nextVertexId++), assemblyGraph);
    createReverseComplementVertex(v0New);

    // Create the edge for the new Segment, keeping the same id.
    edge_descriptor eNew;
    tie(eNew, ignore)= add_edge(v0New, v1, edgeOld, assemblyGraph);

    // Also create the new reverse complement segment. Keep the same id.
    createReverseComplementEdge(eNew, id(eOldRc));

    // Remove eOld and its reverse complement.
    boost::remove_edge(eOld, assemblyGraph);
    boost::remove_edge(eOldRc, assemblyGraph);

    return eNew;
}



AssemblyGraph::edge_descriptor AssemblyGraph::createDisconnectedCopy(edge_descriptor oldSegment)
{
    AssemblyGraph& assemblyGraph = *this;

    // Get the assembly graph edge and vertices of the oldSegment.
    const AssemblyGraphEdge& oldEdge = assemblyGraph[oldSegment];
    const AssemblyGraph::vertex_descriptor oldV0 = source(oldSegment, assemblyGraph);
    const AssemblyGraph::vertex_descriptor oldV1 = target(oldSegment, assemblyGraph);
    const AssemblyGraphVertex& oldVertex0 = assemblyGraph[oldV0];
    const AssemblyGraphVertex& oldVertex1 = assemblyGraph[oldV1];

    // Create the new vertices.
    AssemblyGraphVertex newVertex0 = oldVertex0;
    AssemblyGraphVertex newVertex1 = oldVertex1;
    newVertex0.id = assemblyGraph.nextVertexId++;
    newVertex1.id = assemblyGraph.nextVertexId++;
    const AssemblyGraph::vertex_descriptor newV0 = add_vertex(newVertex0, assemblyGraph);
    const AssemblyGraph::vertex_descriptor newV1 = add_vertex(newVertex1, assemblyGraph);

    // Create the new Segment.
    Segment newSegment;
    tie(newSegment, ignore)= add_edge(newV0, newV1, oldEdge, assemblyGraph);
    assemblyGraph[newSegment].id = nextEdgeId++;

    return newSegment;

}



// BFS starting at a given vertex, in the specified direction
// (0=forward, 1=backward). The BFS cannot use the forbiddenSegments,
// and stops when a stopSegment is reached.
// It records the stopSegments that are reached.
void AssemblyGraph::bfs(
    vertex_descriptor v,
    uint64_t direction,
    const vector<Segment>& forbiddenSegments,
    const vector<Segment>& stopSegments,
    vector<Segment>& reachableStopSegments) const
{
    if(direction == 0) {
        forwardBfs(v, forbiddenSegments, stopSegments, reachableStopSegments);
    } else {
        backwardBfs(v, forbiddenSegments, stopSegments, reachableStopSegments);
    }
}



void AssemblyGraph::forwardBfs(
    vertex_descriptor vStart,
    const vector<Segment>& forbiddenSegments,
    const vector<Segment>& stopSegments,
    vector<Segment>& reachableStopSegments) const
{
    const AssemblyGraph& assemblyGraph = *this;
    SHASTA2_ASSERT(std::is_sorted(forbiddenSegments.begin(), forbiddenSegments.end(), orderById));
    SHASTA2_ASSERT(std::is_sorted(stopSegments.begin(), stopSegments.end(), orderById));

    reachableStopSegments.clear();

    // Initialize the BFS.
    std::queue<vertex_descriptor> q;
    q.push(vStart);
    std::set<vertex_descriptor> visited;
    visited.insert(vStart);

    // BFS loop.
    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();

        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            if(std::binary_search(stopSegments.begin(), stopSegments.end(), e, orderById)){
                reachableStopSegments.push_back(e);
                continue;
            }
            if(std::binary_search(forbiddenSegments.begin(), forbiddenSegments.end(), e, orderById)){
                continue;
            }
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if(not visited.contains(v1)) {
                q.push(v1);
                visited.insert(v1);
            }
        }
    }

    sort(reachableStopSegments.begin(), reachableStopSegments.end(), orderById);
}



void AssemblyGraph::backwardBfs(
    vertex_descriptor vStart,
    const vector<Segment>& forbiddenSegments,
    const vector<Segment>& stopSegments,
    vector<Segment>& reachableStopSegments) const
{
    const AssemblyGraph& assemblyGraph = *this;
    SHASTA2_ASSERT(std::is_sorted(forbiddenSegments.begin(), forbiddenSegments.end(), orderById));
    SHASTA2_ASSERT(std::is_sorted(stopSegments.begin(), stopSegments.end(), orderById));

    reachableStopSegments.clear();

    // Initialize the BFS.
    std::queue<vertex_descriptor> q;
    q.push(vStart);
    std::set<vertex_descriptor> visited;
    visited.insert(vStart);

    // BFS loop.
    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();

        BGL_FORALL_INEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            if(std::binary_search(stopSegments.begin(), stopSegments.end(), e, orderById)){
                reachableStopSegments.push_back(e);
                continue;
            }
            if(std::binary_search(forbiddenSegments.begin(), forbiddenSegments.end(), e, orderById)){
                continue;
            }
            const vertex_descriptor v1 = source(e, assemblyGraph);
            if(not visited.contains(v1)) {
                q.push(v1);
                visited.insert(v1);
            }
        }
    }

    sort(reachableStopSegments.begin(), reachableStopSegments.end(), orderById);
}
