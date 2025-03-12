#include "AssemblyGraph.hpp"
#include "AnchorGraph.hpp"
#include "deduplicate.hpp"
#include "findLinearChains.hpp"
using namespace shasta;



AssemblyGraph::AssemblyGraph(const AnchorGraph& anchorGraph)
{
    AssemblyGraph& assemblyGraph = *this;

    // Find linear chains in the AnchorGraph.
    vector< vector<AnchorGraph::edge_descriptor> > chains;
    findLinearChains(anchorGraph, 0, chains);
    cout << "Found " << chains.size() << " linear chains in the anchor graph." << endl;

    // Get the initial and final AnchorId of all chains.
    // These generate an AssemblyGraph vertex each, after deduplication.
    vector<AnchorId> vertexAnchorIds;
    for(const auto& chain: chains) {
        const AnchorGraph::edge_descriptor eA = chain.front();
        const AnchorGraph::edge_descriptor eB = chain.back();
        const AnchorId anchorIdA = anchorGraph[eA].anchorIdA;
        const AnchorId anchorIdB = anchorGraph[eB].anchorIdB;
        vertexAnchorIds.push_back(anchorIdA);
        vertexAnchorIds.push_back(anchorIdB);
    }
    deduplicate(vertexAnchorIds);

    // Create the vertices.
    std::map<AnchorId, vertex_descriptor> vertexMap;
    for(const AnchorId anchorId: vertexAnchorIds) {
        const vertex_descriptor v = add_vertex(AssemblyGraphVertex(anchorId), assemblyGraph);
        vertexMap.insert(make_pair(anchorId, v));
    }



    // Create the edges.
    for(const auto& chain: chains) {
        const AnchorGraph::edge_descriptor eA = chain.front();
        const AnchorGraph::edge_descriptor eB = chain.back();
        const AnchorId anchorIdA = anchorGraph[eA].anchorIdA;
        const AnchorId anchorIdB = anchorGraph[eB].anchorIdB;
        const vertex_descriptor vA = vertexMap[anchorIdA];
        const vertex_descriptor vB = vertexMap[anchorIdB];

        edge_descriptor e;
        tie(e, ignore) = add_edge(vA, vB, assemblyGraph);
        AssemblyGraphEdge& edge = assemblyGraph[e];
    }

    cout << "The initial assembly graph has " << num_vertices(*this) <<
        " vertices and " << num_edges(*this) << " edges." << endl;
}
