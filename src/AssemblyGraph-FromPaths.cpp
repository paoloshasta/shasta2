#include "AssemblyGraph.hpp"
#include "Anchor.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "performanceLog.hpp"
#include "ReadFollowing.hpp"
#include "RestrictedAnchorGraph.hpp"
#include "TangleMatrix1.hpp"
#include "timestamp.hpp"
using namespace shasta2;

#include <boost/graph/iteration_macros.hpp>



// Note assemblyPaths are not necessarily paths in the AssemblyGraph.
// There may be jumps, which are bridged using local assemblies.
void AssemblyGraph::findAssemblyPaths(vector< vector<edge_descriptor> >& assemblyPaths) const
{
	ReadFollowing::Graph readFollowingGraph(*this);
	readFollowingGraph.findPaths(assemblyPaths);
}



// Note assemlbyPaths are not necessarily paths in the AssemblyGraph.
// There may be jumps, which are bridged using local assemblies.
// The assembly paths must satisfy the following rules:
// - They must consist of at least 2 edge_descriptors.
// - If an edge_descriptor appears at the beginning or end of an assembly path,
//   it cannot appear elsewhere, in other paths or in the same path.
// - Edges that appear inside assembly paths can appear multiple times
//   in the same assembly paths or different assembly paths.
// These rules are satisfied when the assembly paths are computed using
// findAssemblyPaths.
void AssemblyGraph::connectAssemblyPaths(const vector< vector<edge_descriptor> >&  assemblyPaths)
{
	AssemblyGraph& assemblyGraph = *this;

	const bool debug = false;
	if(debug) {
		cout << "Connecting " << assemblyPaths.size() << " assembly paths:" << endl;
		for(const vector<edge_descriptor>& assemblyPath: assemblyPaths) {
			SHASTA2_ASSERT(assemblyPath.size() > 1);
			const edge_descriptor e0 = assemblyPath.front();
			const edge_descriptor e1 = assemblyPath.back();
			cout << "Assembly path with " << assemblyPath.size() << " segments beginning at " <<
				assemblyGraph[e0].id << " and ending at " <<
				assemblyGraph[e1].id << endl;
		}
	}

	// Gather edge_descriptors that appear at the beginning/end of assembly paths.
	std::set<edge_descriptor> pathInitialSegments;
	std::set<edge_descriptor> pathFinalSegments;
	for(const vector<edge_descriptor>& assemblyPath: assemblyPaths) {
		SHASTA2_ASSERT(assemblyPath.size() > 1);
		const edge_descriptor e0 = assemblyPath.front();
		const edge_descriptor e1 = assemblyPath.back();
		SHASTA2_ASSERT(e0 != e1);
		SHASTA2_ASSERT(not pathInitialSegments.contains(e0));
		SHASTA2_ASSERT(not pathFinalSegments.contains(e0));
		SHASTA2_ASSERT(not pathInitialSegments.contains(e1));
		SHASTA2_ASSERT(not pathFinalSegments.contains(e1));
		pathInitialSegments.insert(e0);
		pathFinalSegments.insert(e1);
	}

	// Gather edge_descriptors that appear internally to assembly paths.
	// Count how many times each of them appear.
	// The ones that appear only once will keep their id.
	std::map<edge_descriptor, uint64_t> pathInternalSegments;
	for(const vector<edge_descriptor>& assemblyPath: assemblyPaths) {
		for(uint64_t i=1; i<assemblyPath.size()-1; i++) {
			const edge_descriptor e = assemblyPath[i];
			SHASTA2_ASSERT(not pathInitialSegments.contains(e));
			SHASTA2_ASSERT(not pathFinalSegments.contains(e));
			const auto it = pathInternalSegments.find(e);
			if(it == pathInternalSegments.end()) {
				pathInternalSegments.insert({e, 1});
			} else {
				++it->second;
			}
		}
	}

	if(debug) {
		for(const auto& p: pathInternalSegments) {
			if(p.second > 1) {
				const edge_descriptor e = p.first;
				cout << "Segment " << assemblyGraph[e].id <<
					" appears more than once in assembly paths." << endl;
			}
		}
	}



	// Each assembly path generates a linear sequence of new edges that
	// can be later collapsed into a single edge by calling compress.
	for(const vector<edge_descriptor>& assemblyPath: assemblyPaths) {

		// Generate the new edges.
		vector<edge_descriptor> newAssemblyPath;
		for(uint64_t i=0; i<assemblyPath.size(); i++) {
			const edge_descriptor e = assemblyPath[i];
			const AssemblyGraphEdge& edge = assemblyGraph[e];
			const vertex_descriptor v0 = source(e, assemblyGraph);
			const vertex_descriptor v1 = target(e, assemblyGraph);
			const AnchorId anchorId0 = assemblyGraph[v0].anchorId;
			const AnchorId anchorId1 = assemblyGraph[v1].anchorId;

			// Define the source and target vertices of the new edge, and its id.
			vertex_descriptor v0New = null_vertex();
			vertex_descriptor v1New = null_vertex();
			uint64_t idNew = invalid<uint64_t>;
			if(i == 0) {
				// Initial segment.
				v0New = v0;
				v1New = add_vertex(AssemblyGraphVertex(anchorId1, nextVertexId++), assemblyGraph);
				idNew = edge.id;
			} else if(i == assemblyPath.size() - 1) {
				// Final segment.
				v0New = add_vertex(AssemblyGraphVertex(anchorId0, nextVertexId++), assemblyGraph);
				v1New = v1;
				idNew = edge.id;
			} else {
				// Internal segment.
				v0New = add_vertex(AssemblyGraphVertex(anchorId0, nextVertexId++), assemblyGraph);
				v1New = add_vertex(AssemblyGraphVertex(anchorId1, nextVertexId++), assemblyGraph);
				if(pathInternalSegments[e] == 1) {
					idNew = edge.id;
				} else {
					idNew = nextEdgeId++;
				}
			}

			// Create the new edge.
	        edge_descriptor eNew;
	        tie(eNew, ignore) = add_edge(v0New, v1New, AssemblyGraphEdge(idNew), assemblyGraph);
	        AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];
	        edgeNew = edge;
	        edgeNew.id = idNew;
	        newAssemblyPath.push_back(eNew);
		}



		if(debug) {
		    cout << "Old assembly path " << assemblyGraph[assemblyPath.front()].id << "..." <<
		        assemblyGraph[assemblyPath.back()].id <<
		        " generated new assembly chain " <<
		        assemblyGraph[newAssemblyPath.front()].id << "..." <<
		        assemblyGraph[newAssemblyPath.back()].id << endl;
		}



	    // For each pair of consecutive edges in this path,
	    // generate a new edge in-between to bridge between them.
	    // The code is similar to Tangle1::addConnectPair and Tangle1::detangle,
	    // but simpler.
		for(uint64_t i1=1; i1<newAssemblyPath.size(); i1++) {
			const uint64_t i0 = i1 - 1;
			const edge_descriptor e0 = newAssemblyPath[i0];
			const edge_descriptor e1 = newAssemblyPath[i1];

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
				const TangleMatrix1 tangleMatrix(
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
	}



	// Remove all edges that appear in one or more assembly paths.
	vector<edge_descriptor> edgesToBeRemoved;
	BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
		if(
			pathInitialSegments.contains(e) or
			pathFinalSegments.contains(e) or
			pathInternalSegments.contains(e)) {
			edgesToBeRemoved.push_back(e);
		}
	}
	for(const edge_descriptor e: edgesToBeRemoved) {
		boost::remove_edge(e, assemblyGraph);
	}


	// Write the assembly graph after read following but before compress
	if(debug) {
	    write("ReadFollowing-BeforeCompress");
	}

	// Compress the linear chains we created.
	uint64_t oldCompressDebugLevel = compressDebugLevel;
	if(debug) {
	    compressDebugLevel = 1;
	}
	compress();
    if(debug) {
        compressDebugLevel = oldCompressDebugLevel;
    }
}



void AssemblyGraph::findAndConnectAssemblyPaths()
{
    writePerformanceStatistics("AssemblyGraph::findAndConnectAssemblyPaths begins");

    vector< vector<edge_descriptor> > assemblyPaths;
	findAssemblyPaths(assemblyPaths);
	connectAssemblyPaths(assemblyPaths);

    writePerformanceStatistics("AssemblyGraph::findAndConnectAssemblyPaths ends");
}

