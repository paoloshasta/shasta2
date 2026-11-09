// Shasta.
#include "Tangle.hpp"
#include "Anchor.hpp"
#include "RestrictedAnchorGraph.hpp"
#include "TangleMatrix.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



// Constructor from a set of AssemblyGraph vertices.
Tangle::Tangle(
    AssemblyGraph& assemblyGraph,
    const vector<vertex_descriptor>& tangleVerticesArgument) :
    assemblyGraph(assemblyGraph),
    tangleVertices(tangleVerticesArgument)
{
    // Sort the tangle vertices by id.
    std::ranges::sort(tangleVertices, assemblyGraph.orderById);

    // Find the edges.
    findEntrances();
    findExits();
    findTangleEdges();

    // Compute the tangle matrix.
    ostream html(0);
    tangleMatrixPointer = make_shared<TangleMatrix>(assemblyGraph, entrances, exits, html);
}



// Constructor from a single vertex.
Tangle::Tangle(
    AssemblyGraph& assemblyGraph,
    vertex_descriptor v) :
    Tangle(
        assemblyGraph,
        vector<vertex_descriptor>(1, v))
{}



// Constructor from an edge.
Tangle::Tangle(
    AssemblyGraph& assemblyGraph,
    edge_descriptor e) :
    Tangle(
        assemblyGraph,
        vector<vertex_descriptor>({source(e, assemblyGraph), target(e, assemblyGraph)}))
{}



void Tangle::findEntrances()
{
    entrances.clear();
    for(const vertex_descriptor v1: tangleVertices) {
        BGL_FORALL_INEDGES(v1, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v0 = source(e, assemblyGraph);
            if(not isTangleVertex(v0)) {
                entrances.push_back(e);
            }
        }
    }
    std::ranges::sort(entrances, assemblyGraph.orderById);
}



void Tangle::findExits()
{
    exits.clear();
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if(not isTangleVertex(v1)) {
                exits.push_back(e);
            }
        }
    }
    std::ranges::sort(exits, assemblyGraph.orderById);
}



void Tangle::findTangleEdges()
{
    tangleEdges.clear();
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if(isTangleVertex(v1)) {
                tangleEdges.push_back(e);
            }
        }
    }
    sort(tangleEdges.begin(), tangleEdges.end(), assemblyGraph.orderById);
}



bool Tangle::isTangleVertex(vertex_descriptor v) const
{
    return std::ranges::binary_search(tangleVertices, v, assemblyGraph.orderById);
}



bool Tangle::addConnectPair(uint64_t entranceIndex, uint64_t exitIndex) {
    SHASTA2_ASSERT(entranceIndex < entrances.size());
    SHASTA2_ASSERT(exitIndex < exits.size());
    connectPairs.emplace_back(entranceIndex, exitIndex);
    AssemblyGraphEdge& newEdge = connectPairs.back().newEdge;

    // Store the AssemblyGraphEdge, without adding it to the AssemblyGraph for now.
    const AssemblyGraph::vertex_descriptor v0 = target(entrances[entranceIndex], assemblyGraph);
    const AssemblyGraph::vertex_descriptor v1 = source(exits[exitIndex], assemblyGraph);
    const AnchorId anchorId0 = assemblyGraph[v0].anchorId;
    const AnchorId anchorId1 = assemblyGraph[v1].anchorId;

    if(anchorId0 == anchorId1) {
        // We just generate an empty edge without any steps.
        // If any of these are left after compress, they will have to be
        // removed by collapsing the vertices they join.
        return true;
    } else {

        try {

            // Create the RestrictedAnchorGraph, then:
            // - Remove vertices not accessible from anchorId0 and anchorId1.
            // - Remove cycles.
            // - Find the longest path.
            // - Add one step for each edge of the longest path of the RestrictedAnchorGraph.
            ostream html(0);
            RestrictedAnchorGraph restrictedAnchorGraph(
                assemblyGraph.anchors, assemblyGraph.journeys, tangleMatrix(), entranceIndex, exitIndex, html);
            vector<RestrictedAnchorGraph::edge_descriptor> longestPath;
            // restrictedAnchorGraph.findLongestPath(longestPath);
            restrictedAnchorGraph.findOptimalPath(anchorId0, anchorId1, longestPath);

            for(const RestrictedAnchorGraph::edge_descriptor re: longestPath) {
                const auto& rEdge = restrictedAnchorGraph[re];
                if(rEdge.anchorPair.size() == 0) {
                    newEdge.clear();
                    return false;
                }
                newEdge.push_back(AssemblyGraphEdgeStep(rEdge.anchorPair, rEdge.offset));
            }

        } catch(const std::exception&) {

            // Creation of the RestrictedAnchorGraph failed.
            return false;
        }

        return true;
    }

}



void Tangle::detangle()
{
    // Reroute the entrances to new vertices, so all
    // entrances and exits become temporarily dangling.
    vector<vertex_descriptor> newEntranceVertices;
    rerouteEntrances(newEntranceVertices);

    vector<vertex_descriptor> newExitVertices;
    rerouteExits(newExitVertices);

    // Finally, reconnect the entrance/exit pairs in our ConnectPairs.
    for(ConnectPair& connectPair: connectPairs) {
        const uint64_t entranceIndex = connectPair.entranceIndex;
        const uint64_t exitIndex = connectPair.exitIndex;
        reconnect(
            connectPair,
            newEntranceVertices[entranceIndex],
            newExitVertices[exitIndex]);
    }

    // Now we can remove all the tangle vertices and their edges.
    // This also removes the old entrances and exits.
    removedVertices = tangleVertices;
    for(const vertex_descriptor v: tangleVertices) {
        if(false) {
            cout << "Removing for detangling:";
            BGL_FORALL_INEDGES(v, e, assemblyGraph, AssemblyGraph) {
                cout << " " << assemblyGraph[e].id;
            }
            BGL_FORALL_OUTEDGES(v, e, assemblyGraph, AssemblyGraph) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << endl;
        }
        clear_vertex(v, assemblyGraph);
        remove_vertex(v, assemblyGraph);
    }
}



// Make a copy of each entrance edge, with the target vertex replaced by a new vertex
// with the same AnchorId.
void Tangle::rerouteEntrances(vector<vertex_descriptor>& newEntranceVertices) const
{
    newEntranceVertices.clear();

    for(const edge_descriptor& eOld: entrances) {
        const vertex_descriptor v0Old = source(eOld, assemblyGraph);
        AssemblyGraphEdge& edgeOld = assemblyGraph[eOld];
        const AnchorId lastAnchorId = edgeOld.back().anchorPair.anchorIdB;

        // Create the new vertex.
        const vertex_descriptor v1 = add_vertex(
            AssemblyGraphVertex(lastAnchorId, assemblyGraph.nextVertexId++), assemblyGraph);
        newEntranceVertices.push_back(v1);

        // Create the new edge, with the same steps and id as the old one.
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(v0Old, v1, AssemblyGraphEdge(edgeOld.id), assemblyGraph);
        AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];
        edgeNew.swapSteps(edgeOld);
        edgeNew.wasAssembled = edgeOld.wasAssembled;
    }
}



// Make a copy of each exit edge, with the source vertex replaced by a new vertex
// with the same AnchorId.
void Tangle::rerouteExits(vector<vertex_descriptor>& newExitVertices) const
{
    newExitVertices.clear();
    for(const edge_descriptor& eOld: exits) {
        const vertex_descriptor v1Old = target(eOld, assemblyGraph);
        AssemblyGraphEdge& edgeOld = assemblyGraph[eOld];
        const AnchorId firstAnchorId = edgeOld.front().anchorPair.anchorIdA;

        // Create the new vertex.
        const vertex_descriptor v0 = add_vertex(
            AssemblyGraphVertex(firstAnchorId, assemblyGraph.nextVertexId++), assemblyGraph);
        newExitVertices.push_back(v0);

        // Create the new edge, with the same steps and id as the old one.
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(v0, v1Old, AssemblyGraphEdge(edgeOld.id), assemblyGraph);
        AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];
        edgeNew.swapSteps(edgeOld);
        edgeNew.wasAssembled = edgeOld.wasAssembled;
    }
}



void Tangle::reconnect(
    ConnectPair& connectPair,
    vertex_descriptor v0,
    vertex_descriptor v1
    ) const
{
    const bool debug = false;

    const AnchorId anchorId0 = assemblyGraph[v0].anchorId;
    const AnchorId anchorId1 = assemblyGraph[v1].anchorId;

    if(debug) {
        cout << "Connecting entrance " << assemblyGraph[entrances[connectPair.entranceIndex]].id <<
            " with exit " << assemblyGraph[exits[connectPair.exitIndex]].id << endl;
        cout << "Anchors " << anchorIdToString(anchorId0) << " " <<
            anchorIdToString(anchorId1) << endl;
    }


    // Create a new AssemblyGraphEdge between v0 and v1.
    AssemblyGraph::edge_descriptor e;
    tie(e, ignore) = add_edge(v0, v1, assemblyGraph);
    AssemblyGraphEdge& newEdge = assemblyGraph[e];
    newEdge.swapSteps(connectPair.newEdge);
    newEdge.id = assemblyGraph.nextEdgeId++;
    if(debug) {
        cout << "Created new assembly graph edge " << newEdge.id << endl;
    }
}


// Figure out if the tangle is self-complementary
// (that is, it coincides with its reverse complement).
// This assumes that the AssemblyGraph is strand-symmetric.
bool Tangle::isSelfComplementary() const
{
    const AssemblyGraph::vertex_descriptor v = tangleVertices.front();
    const AssemblyGraph::vertex_descriptor vRc = assemblyGraph[v].vRc;
    return std::ranges::binary_search(tangleVertices, vRc, assemblyGraph.orderById);
}



void Tangle::writeHtml(ostream& html) const
{
    html << "<h2>Tangle entrances</h2>"
        "This tangle has " << entrances.size() << " entrances:<br>";
    for(uint64_t i=0; i<entrances.size(); i++) {
        if(i != 0) {
            html << ",";
        }
        html << assemblyGraph.id(entrances[i]);
    }

    html << "<h2>Tangle exits</h2>"
        "This tangle has " << exits.size() << " exits:<br>";
    for(uint64_t i=0; i<exits.size(); i++) {
        if(i != 0) {
            html << ",";
        }
        html << assemblyGraph.id(exits[i]);
    }

    html << "<h2>Tangle internal segments</h2>"
        "This tangle has " << tangleEdges.size() << " internal segments:<br>"
        "<span style='word-break:break-word'>";
    for(uint64_t i=0; i<tangleEdges.size(); i++) {
        if(i != 0) {
            html << ",";
        }
        html << assemblyGraph.id(tangleEdges[i]);
    }
    html << "</span>";

    html << "<h2>Tangle vertices</h2>"
        "This tangle has " << tangleVertices.size() << " vertices:<br>"
        "<span style='word-break:break-word'>";
    for(uint64_t i=0; i<tangleVertices.size(); i++) {
        if(i != 0) {
            html << ",";
        }
        html << assemblyGraph.id(tangleVertices[i]);
    }
    html << "</span>";

    html << "<h2>Tangle matrix</h2>"
        "Zero values are omitted."
        "<br><br><table><tr><th>";
    html << std::fixed << std::setprecision(1);
    for(uint64_t i=0; i<exits.size(); i++) {
        html << "<th>" << assemblyGraph.id(exits[i]);
    }
    for(uint64_t i=0; i<entrances.size(); i++) {
        html << "<tr><th>" << assemblyGraph.id(entrances[i]);
        for(uint64_t j=0; j<exits.size(); j++) {
            html << "<td class=centered>";
            const double t = tangleMatrix().tangleMatrix[i][j];
            if(t != 0.) {
                html << t;
            }
        }
    }
    html << "</table>";

}



// Create a map that maps ids to edge_descriptors (Segments)
// for all the entrances and exits plus their reverse complements.
// This is used in detangling.
void Tangle::createSegmentMap(std::map<uint64_t, Segment>& segmentMap) const
{
    segmentMap.clear();
    const bool selfComplementary = isSelfComplementary();

    for(const Segment e: entrances) {
        segmentMap.insert({assemblyGraph.id(e), e});
        if(not selfComplementary) {
            const Segment eRc = assemblyGraph[e].eRc;
            segmentMap.insert({assemblyGraph.id(eRc), eRc});
        }
    }

    for(const Segment e: exits) {
        segmentMap.insert({assemblyGraph.id(e), e});
        if(not selfComplementary) {
            const Segment eRc = assemblyGraph[e].eRc;
            segmentMap.insert({assemblyGraph.id(eRc), eRc});
        }
    }
}



void Tangle::checkSegmentMap(const std::map<uint64_t, Segment>& segmentMap) const
{
    for(const auto&[segmentId, e]: segmentMap) {
        SHASTA2_ASSERT(segmentId == assemblyGraph.id(e));
    }

}



bool Tangle::isDegenerate() const
{
    vector<Segment> allSegments;
    std::ranges::copy(entrances, back_inserter(allSegments));    
    std::ranges::copy(exits, back_inserter(allSegments));    
    std::ranges::copy(tangleEdges, back_inserter(allSegments));    
    sort(allSegments.begin(), allSegments.end(), assemblyGraph.orderById);

    // Check for duplicates.
    if(std::ranges::adjacent_find(allSegments) != allSegments.end()) {

        for(uint64_t i=1; i<allSegments.size(); i++) {
            if(allSegments[i-1] == allSegments[i]) {
                cout << "Duplicate " << assemblyGraph.id(allSegments[i]) << endl;
            }
        }
        return true;
    }

#if 0
    // TURN THIS OFF BECAUSE IT CAUSES ALL SELF-COMPLEMENTARY TANGLES
    // TO BE CONSIDERED DEGENERATE.
    // Check for the presence of reverse complement Segments.
    for(const Segment segment: allSegments) {
        const Segment segmentRc = assemblyGraph[segment].eRc;
        if(std::binary_search(allSegments.begin(), allSegments.end(),segmentRc, assemblyGraph.orderById)) {
            return true;
        }
    }
#endif

    return false;
}
