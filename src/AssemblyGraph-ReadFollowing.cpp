#include "AssemblyGraph.hpp"
#include "findDescendants.hpp"
#include "Options.hpp"
#include "ReadFollowing5.hpp"
#include "SegmentGraph.hpp"
#include "Tangle.hpp"
using namespace shasta2;
using namespace ReadFollowing5;



bool AssemblyGraph::readFollowingStrandSymmetric(
    uint64_t tangleId,
    const Tangle& tangle,
    ostream& html)
{
    SHASTA2_ASSERT(not tangle.entrances.empty());
    SHASTA2_ASSERT(not tangle.exits.empty());

    AssemblyGraph& assemblyGraph = tangle.assemblyGraph;
    const bool isSelfComplementary = tangle.isSelfComplementary();

    // For now, only do it for tangles with a number of entrances
    // equal to the number of exits.
    if(tangle.entrances.size() != tangle.exits.size()) {
        if(html) {
            html << "<br>Skipping read following because the number of entrances is not the same as the number of exits.";
        }
        return false;
    }

    // Create the ReadFollowing5::Graph.
    const Graph graph(*this, tangleId, tangle, html);

    // For now, require a one-on-one matching of the entrances
    // and exits.
    for(uint64_t iEntrance=0; iEntrance<tangle.entrances.size(); iEntrance++) {
        if(graph.countReachableExits(iEntrance) != 1) {
            if(html) {
                html << "<br>Skipping read following because of ambiguous or missing reachability.";
            }
            return false;
        }
    }
    for(uint64_t iExit=0; iExit<tangle.exits.size(); iExit++) {
        if(graph.countReachableEntrances(iExit) != 1) {
            if(html) {
                html << "<br>Skipping read following because of ambiguous or missing reachability.";
            }
            return false;
        }
    }

    // If getting here, read following can proceed.
    if(html) {
        html << "<br>Read following can proceed for this tangle.";
    }



    // We need to disconnect (at the beginning or end) the entrances and exits.
    // To do this we use disconnectAtBeginning and disconnectAtEnd, which change
    // edge_descriptors but leave the ids unchanged.
    // Some Segments can at the same time be an entrance and/or exit
    // of this tangle or its reverse complement. To avoid working with invalidated
    // edge_descriptors, we work with Segment ids instead.
    std::map<uint64_t, Segment> segmentMap;
    tangle.createSegmentMap(segmentMap);
    vector<uint64_t> entranceIds;
    for(const Segment e: tangle.entrances) {
        entranceIds.push_back(id(e));
    }
    vector<uint64_t> exitIds;
    for(const Segment e: tangle.exits) {
        exitIds.push_back(id(e));
    }

    // Create disconnected versions of the entrances and exits,
    // while keeping the segmentMap up to date.
    for(const uint64_t entranceId: entranceIds) {
        const Segment eNew = disconnectAtEnd(segmentMap.at(entranceId));
        segmentMap[id(eNew)] = eNew;
        const Segment eNewRc = assemblyGraph[eNew].eRc;
        segmentMap[id(eNewRc)] = eNewRc;
    }
    if(not isSelfComplementary) {
        for(const uint64_t exitId: exitIds) {
            const Segment eNew = disconnectAtBeginning(segmentMap.at(exitId));
            segmentMap[id(eNew)] = eNew;
            const Segment eNewRc = assemblyGraph[eNew].eRc;
            segmentMap[id(eNewRc)] = eNewRc;
        }
    }
    tangle.checkSegmentMap(segmentMap);



    // Now we can connect entrances with exits using the shortest paths stored in the Graph.
    std::set<uint64_t> connectedEntrances;
    for(uint64_t iEntrance=0; iEntrance<tangle.entrances.size(); iEntrance++) {
        const uint64_t entranceId = entranceIds[iEntrance];
        if(isSelfComplementary and (connectedEntrances.contains(entranceId))) {
            continue;
        }
        const Segment entrance = segmentMap.at(entranceId);
        for(uint64_t iExit=0; iExit<tangle.exits.size(); iExit++) {
            const uint64_t exitId = exitIds[iExit];
            const Segment exit = segmentMap.at(exitId);
            const auto& shortestPath = graph.shortestPaths[iEntrance][iExit];
            if(shortestPath.isUnreachable()) {
                continue;
            }

            // Get the Segments internal to the path.
            // We will make disconnected copies of these segments.
            const vector<Segment>& oldAssemblyPath = shortestPath.segments;



            // Create a complete path consisting of the entrance, the
            // newly created copies of the internal segments, and the exit.
            // We will still have to fill in the portions in-between these segments.
            vector<Segment> assemblyPath;
            assemblyPath.push_back(entrance);
            for(const Segment oldSegment: oldAssemblyPath) {
                const Segment newSegment = assemblyGraph.createDisconnectedCopy(oldSegment);
                assemblyPath.push_back(newSegment);
                const vertex_descriptor v0 = source(newSegment, assemblyGraph);
                const vertex_descriptor v1 = target(newSegment, assemblyGraph);
                createReverseComplementVertex(v0);
                createReverseComplementVertex(v1);
                createReverseComplementEdge(newSegment);
            }
            assemblyPath.push_back(exit);

            // Now we have to fill in the portions in-between these segments.
            for(uint64_t i1=1; i1<assemblyPath.size(); i1++) {
                const uint64_t i0 = i1 - 1;
                const Segment segment0 = assemblyPath[i0];
                const Segment segment1 = assemblyPath[i1];
                const Segment newSegment = connect(segment0, segment1);
                createReverseComplementEdge(newSegment);
            }



            if(isSelfComplementary) {
                connectedEntrances.insert(entranceId);
                connectedEntrances.insert(id(assemblyGraph[exit].eRc));
            }
        }
    }



    // Remove all the Segments internal to this tangle
    // and their reverse complements.
    for(const Segment e: tangle.tangleEdges) {
        if(not isSelfComplementary) {
            boost::remove_edge(assemblyGraph[e].eRc, assemblyGraph);
        }
        boost::remove_edge(e, assemblyGraph);
    }

    // If getting here, read following was successful on this tangle.
    return true;
}
