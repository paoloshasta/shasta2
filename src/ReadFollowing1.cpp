// Read following in the AssemblyGraph.

// Shasta.
#include "ReadFollowing1.hpp"
#include "deduplicate.hpp"
#include "findLinearChains.hpp"
#include "Journeys.hpp"
#include "Markers.hpp"
#include "Options.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
#include "transitiveReduction.hpp"
using namespace shasta2;
using namespace ReadFollowing1;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"
#include <iomanip>
#include <random>



Graph::Graph(const AssemblyGraph& assemblyGraph) :
    assemblyGraph(assemblyGraph)
{
    Graph& graph = *this;

    // Initial creation with all possible vertices and edges.
    createVertices();
    createEdges();
    cout << "The initial read following graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;
    write("Initial");


    // Prune short leaves.
    prune();
   	write("Final");
    cout << "After pruning, the read following graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;

    // This fills the outEdges and inEdges fields in each vertex.
    // They are needed for efficient generation of random paths.
    // This must be called after all changes to the graph have been made.
    fillConnectivity();

}



// Create vertices of the ReadFollowing1 graph.
// Each vertex corresponds to a Segment of the AssemblyGraph.
void Graph::createVertices()
{
    // EXPOSE WHEN CODE STABILIZES.
    // const double maxCoverage = 18.;

    Graph& graph = *this;

    // Each Segment generates a Vertex.
    BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {
        const Vertex vertex(assemblyGraph, segment);
        if(true) {
            const vertex_descriptor v = add_vertex(Vertex(assemblyGraph, segment), graph);
            vertexMap.insert(make_pair(segment, v));
        }
    }

}



Vertex::Vertex(
    const AssemblyGraph& assemblyGraph,
    Segment segment) :
    segment(segment)
{
    const AssemblyGraphEdge& edge = assemblyGraph[segment];
    if(edge.wasAssembled) {
        length = edge.sequenceLength();
    } else {
        length = edge.offset();
    }

    coverage = edge.lengthWeightedAverageCoverage();

    // Compute initial/final support.
    const uint32_t representativeRegionStepCount =  uint32_t(assemblyGraph.options.representativeRegionStepCount);
    SegmentStepSupport::getInitialFirst(assemblyGraph, segment, representativeRegionStepCount, initialSupport);
    SegmentStepSupport::getFinalLast   (assemblyGraph, segment, representativeRegionStepCount, finalSupport  );
}



// Create edges of the ReadFollowing1 graph.
// An edge v0->v1 is created if the final support of v0
// shares at least one OrientedReadId with the initial support of v1.
void Graph::createEdges()
{
    performanceLog << timestamp << "ReadFollowing1::Graph::createEdges begins." << endl;

    const bool debug = false;

    Graph& graph = *this;
    const uint64_t orientedReadCount = assemblyGraph.journeys.size();
    const uint64_t minCommonCount = assemblyGraph.options.readFollowingMinCommonCount;
    const double minCorrectedJaccard = assemblyGraph.options.readFollowingMinCorrectedJaccard;

    // For each OrientedReadId, gather the vertices that the OrientedReadId
    // appears in, in the initial/final support.
    vector< vector<vertex_descriptor> > initialSupportVertices(orientedReadCount);
    vector< vector<vertex_descriptor> > finalSupportVertices(orientedReadCount);
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Vertex& vertex = graph[v];

        for(const SegmentStepSupport& s: vertex.initialSupport) {
            initialSupportVertices[s.orientedReadId.getValue()].push_back(v);
        }

        for(const SegmentStepSupport& s: vertex.finalSupport) {
            finalSupportVertices[s.orientedReadId.getValue()].push_back(v);
        }
    }



    // An edge v0->v1 will be created if the final support of v0
    // shares at least minCommonCount OrientedReadId with the initial support of v1.
    vector< pair<vertex_descriptor, vertex_descriptor> > vertexPairs;
    for(uint64_t i=0; i<orientedReadCount; i++) {
        const vector<vertex_descriptor>& initialVertices = initialSupportVertices[i];
        const vector<vertex_descriptor>& finalVertices = finalSupportVertices[i];
        for(const vertex_descriptor v0: finalVertices) {
            for(const vertex_descriptor v1: initialVertices) {
                if(v1 == v0) {
                    continue;
                }

                // This OrientedReadId appears in the final support of v0 and in the
                // initial support of v1, so we will create an edge v0->v1.
                vertexPairs.push_back({v0, v1});
            }
        }
    }

    if(debug) {
        // Debug.
        vector<uint64_t> countCopy;
        vector< pair<vertex_descriptor, vertex_descriptor> > vertexPairsCopy = vertexPairs;
        deduplicateAndCount(vertexPairsCopy, countCopy);
        ofstream csv("ReadFollowing-VertexPairs.csv");
        for(uint64_t i=0; i<vertexPairsCopy.size(); i++) {
            const vertex_descriptor v0 = vertexPairsCopy[i].first;
            const vertex_descriptor v1 = vertexPairsCopy[i].second;
            const Segment segment0 = graph[v0].segment;
            const Segment segment1 = graph[v1].segment;
            csv << assemblyGraph[segment0].id << ",";
            csv << assemblyGraph[segment1].id << ",";
            csv << countCopy[i] << "\n";
        }
    }

    vector<uint64_t> count;
    deduplicateAndCountWithThreshold(vertexPairs, count, minCommonCount);
    SHASTA2_ASSERT(vertexPairs.size() == count.size());



    // Generate an edge for each of these pairs that satisfies our requirements.
    // This should be parallelized because the call to canConnect is expensive.
    cout << "Found " << vertexPairs.size() << " candidate edges for the read following graph." << endl;
    for(uint64_t i=0; i<vertexPairs.size(); i++) {
        const auto& p = vertexPairs[i];

        const vertex_descriptor v0 = p.first;
        const vertex_descriptor v1 = p.second;

        const Segment segment0 = graph[v0].segment;
        const Segment segment1 = graph[v1].segment;

        // Create the candidate edge.
        const Edge edge(assemblyGraph, segment0, segment1);

        // This must be true given the way we constructed the vertex pairs.
        SHASTA2_ASSERT(edge.segmentPairInformation.commonCount >= minCommonCount);
        SHASTA2_ASSERT(edge.segmentPairInformation.commonCount == count[i]);

        if(edge.segmentPairInformation.segmentOffset < 0) {
            if(debug) {
                cout << "Candidate edge " << assemblyGraph[segment0].id << " " << assemblyGraph[segment1].id <<
                    " discarded to negative offset." << endl;
            }
            continue;
        }

        if(edge.segmentPairInformation.correctedJaccard < minCorrectedJaccard) {
            if(debug) {
                cout << "Candidate edge " << assemblyGraph[segment0].id << " " << assemblyGraph[segment1].id <<
                    " discarded to low corrected jaccard " << edge.segmentPairInformation.correctedJaccard << endl;
            }
            continue;
        }

        if(not assemblyGraph.canConnect(segment0, segment1)) {
            if(debug) {
                cout << "Candidate edge " << assemblyGraph[segment0].id << " " << assemblyGraph[segment1].id <<
                    " discarded because they cannot be connected for assembly." << endl;
            }
            continue;
        }

        // Add this edge to the Graph.
        add_edge(v0, v1, edge, graph);
    }
    cout << "Kept " << num_edges(graph) << " edges in the read following graph." << endl;

    performanceLog << timestamp << "ReadFollowing1::Graph::createEdges ends." << endl;
}



Edge::Edge(
    const AssemblyGraph& assemblyGraph,
    Segment segment0,
    Segment segment1)
{
    const uint32_t representativeRegionStepCount =  uint32_t(assemblyGraph.options.representativeRegionStepCount);
    ostream html(0);
    segmentPairInformation = SegmentStepSupport::analyzeSegmentPair(
        html, assemblyGraph, segment0, segment1, representativeRegionStepCount);

    SHASTA2_ASSERT(segmentPairInformation.commonCount > 0);
}



void Graph::write(const string& name) const
{
    cout << "ReadFollowing1-" << name << ": " << num_vertices(*this) <<
        " vertices, " << num_edges(*this) << " edges." << endl;
    writeGraphviz(name);
    writeCsv(name);
}



void Graph::writeGraphviz(const string& name) const
{
    const Graph& graph = *this;

    ofstream dot("ReadFollowing-" + name + ".dot");
    dot << "digraph ReadFollowing {\n";
    dot << std::fixed << std::setprecision(1);

    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Vertex& vertex = graph[v];
        const Segment segment = vertex.segment;
        const AssemblyGraphEdge& assemblyGraphEdge = assemblyGraph[segment];
        dot << assemblyGraphEdge.id;

        // Begin attributes.
        dot << " [";

        // Label.
        dot <<
            "label=\"" << assemblyGraphEdge.id << "\\n" <<
            vertex.length << "\\n" <<
            vertex.coverage << "x" <<
            "\"";

        // Color.
        if(vertex.length >= assemblyGraph.options.readFollowingSegmentLengthThreshold) {
            dot << " style=filled fillcolor=cyan";
        }

        // End attributes.
        dot << "]";

        // End the line for this vertex.
        dot << ";\n";
    }



    dot << std::fixed << std::setprecision(2);
    BGL_FORALL_EDGES(e, graph, Graph) {
        const Edge& edge = graph[e];
        // const int32_t offset = edge.segmentPairInformation.segmentOffset;

        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const Segment segment0 = graph[v0].segment;
        const Segment segment1 = graph[v1].segment;

        dot << assemblyGraph[segment0].id << "->" << assemblyGraph[segment1].id;

        // Begin attributes.
        dot << "[";

        // Label.
#if 0
        dot << "label=\"" <<
            edge.segmentPairInformation.commonCount << "/" <<
            std::fixed << std::setprecision(2) <<
            edge.segmentPairInformation.correctedJaccard << "\\n" <<
            offset << "\"";
#endif

        // Tooltip.
        dot << " tooltip=\"" <<
            edge.segmentPairInformation.commonCount <<  "\"";

        // Thickness is proportional to commonCount.
        dot << " penwidth=" << 0.2 * double(edge.segmentPairInformation.commonCount);

#if 0
        // Color is determined by correctedJaccard for the edge.
        // Green = 1
        // Red = assemblyGraph.options.readFollowingMinCorrectedJaccard.
        double hue;
        if(edge.segmentPairInformation.correctedJaccard >= 1.) {
            hue = 1.;
        } else if(edge.segmentPairInformation.correctedJaccard <= assemblyGraph.options.readFollowingMinCorrectedJaccard) {
            hue = 0.;
        } else {
            hue =
                (edge.segmentPairInformation.correctedJaccard - assemblyGraph.options.readFollowingMinCorrectedJaccard) /
                (1. - assemblyGraph.options.readFollowingMinCorrectedJaccard);
        }
        hue /= 3.;
        dot << std::fixed << std::setprecision(3) << " color=\""  << hue << " 1. 1.\"";
#endif

        // End attributes.
        dot << "]";

        // End the line for this edge.
        dot << ";\n";

    }

    dot << "}\n";
}



void Graph::writeCsv(const string& name) const
{
    writeVerticesCsv(name);
    writeEdgesCsv(name);
}



void Graph::writeVerticesCsv(const string& name) const
{
    const Graph& graph = *this;

    ofstream csv("ReadFollowing-Vertices-" + name + ".csv");
    csv << "Segment,Length,InitialSupport,FinalSupport,\n";

    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Vertex& vertex = graph[v];
        const Segment segment = vertex.segment;
        const AssemblyGraphEdge& assemblyGraphEdge = assemblyGraph[segment];

        csv << assemblyGraphEdge.id << ",";
        csv << vertex.length << ",";
        csv << vertex.initialSupport.size() << ",";
        csv << vertex.finalSupport.size() << ",";
        csv << "\n";
    }
}



void Graph::writeEdgesCsv(const string& name) const
{
    const Graph& graph = *this;

    ofstream csv("ReadFollowing-Edges-" + name + ".csv");
    csv << "Segment0,Segment1,Length0,Length1,FinalSupport0,InitialSupport1,"
        "Common,Missing0,Missing1,MissingTotal,CorrectedJaccard,Offset,\n";
    csv << std::fixed << std::setprecision(2);

    BGL_FORALL_EDGES(e, graph, Graph) {
        const Edge& edge = graph[e];

        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        const Vertex& vertex0 = graph[v0];
        const Vertex& vertex1 = graph[v1];

        const Segment segment0 = vertex0.segment;
        const Segment segment1 = vertex1.segment;

        csv << assemblyGraph[segment0].id << ",";
        csv << assemblyGraph[segment1].id << ",";
        csv << vertex0.length << ",";
        csv << vertex1.length << ",";
        csv << vertex0.finalSupport.size() << ",";
        csv << vertex1.initialSupport.size() << ",";
        csv << edge.segmentPairInformation.commonCount << ",";
        csv << edge.segmentPairInformation.missing0 << ",";
        csv << edge.segmentPairInformation.missing1 << ",";
        csv << edge.segmentPairInformation.missing0 + edge.segmentPairInformation.missing1 << ",";
        csv << edge.segmentPairInformation.correctedJaccard << ",";
        csv << edge.segmentPairInformation.segmentOffset << ",";
        csv << "\n";
    }
}



template<std::uniform_random_bit_generator RandomGenerator> void Graph::findRandomPath(
    vertex_descriptor v, uint64_t direction,
    RandomGenerator& randomGenerator,
    vector<vertex_descriptor>& path,
    const std::set<vertex_descriptor>& stopVertices)
{
    if(direction == 0) {
        findRandomForwardPath(v, randomGenerator, path, stopVertices);
    } else {
        findRandomBackwardPath(v, randomGenerator, path, stopVertices);
    }
}



template<std::uniform_random_bit_generator RandomGenerator> void Graph::findRandomForwardPath(
    vertex_descriptor v,
    RandomGenerator& randomGenerator,
    vector<vertex_descriptor>& path,
    const std::set<vertex_descriptor>& stopVertices)
{
    const Graph& graph = *this;

    // Start with a path consisting of just this vertex.
    path.clear();
    path.push_back(v);

    // At each iteration, add one vertex to the path.
    while(out_degree(v, assemblyGraph) > 0) {

        // Pick a random out-edge.
        const vector<edge_descriptor>& outEdges = graph[v].outEdges;
        const uint64_t n = outEdges.size();
        std::uniform_int_distribution<uint64_t> distribution(0, n - 1);
        const uint64_t i = distribution(randomGenerator);
        const edge_descriptor e = outEdges[i];

        // Add to the path the target of this edge and continue from here.
        v = target(e, graph);
        path.push_back(v);

        if(stopVertices.contains(v)) {
            break;
        }
    }
}



template<std::uniform_random_bit_generator RandomGenerator> void Graph::findRandomBackwardPath(
    vertex_descriptor v,
    RandomGenerator& randomGenerator,
    vector<vertex_descriptor>& path,
    const std::set<vertex_descriptor>& stopVertices)
{
    const Graph& graph = *this;

    // Start with a path consisting of just this vertex.
    path.clear();
    path.push_back(v);

    // At each iteration, add one vertex to the path.
    while(in_degree(v, assemblyGraph) > 0) {

        // Pick a random in-edge.
        const vector<edge_descriptor>& inEdges = graph[v].inEdges;
        const uint64_t n = inEdges.size();
        std::uniform_int_distribution<uint64_t> distribution(0, n - 1);
        const uint64_t i = distribution(randomGenerator);
        const edge_descriptor e = inEdges[i];

        // Add to the path the source of this edge and continue from here.
        v = source(e, graph);
        path.push_back(v);

        if(stopVertices.contains(v)) {
            break;
        }
    }

    // Reverse the path so it goes forward.
    std::ranges::reverse(path);
}



void Graph::writeRandomPath(Segment segment, uint64_t direction)
{
    const Graph& graph = *this;

    const auto it = vertexMap.find(segment);
    SHASTA2_ASSERT(it != vertexMap.end());
    const vertex_descriptor v = it->second;

    vector<vertex_descriptor> path;
    std::set<vertex_descriptor> stopVertices;
    std::random_device randomGenerator;
    findRandomPath(v, direction, randomGenerator, path, stopVertices);

    cout << "Found a path of length " << path.size() << ":" << endl;
    for(const vertex_descriptor v: path) {
        const Segment segment = graph[v].segment;
        cout << assemblyGraph[segment].id << ",";
    }
    cout << endl;
}



void Graph::prune()
{
    while(pruneIteration());
}



bool Graph::pruneIteration()
{
    Graph& graph = *this;

    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if(graph[v].length < assemblyGraph.options.readFollowingPruneLength) {
            const bool isLeaf = (in_degree(v, graph) == 0) or (out_degree(v, graph) == 0);
            if(isLeaf) {
                verticesToBeRemoved.push_back(v);
            }
        }
    }

    for(const vertex_descriptor v: verticesToBeRemoved) {
        boost::clear_vertex(v, graph);
        boost::remove_vertex(v, graph);
    }

    if(not verticesToBeRemoved.empty()) {
        cout << "Pruned " << verticesToBeRemoved.size() << " vertices." << endl;
    }

    return not verticesToBeRemoved.empty();
}



#if 0
// This version results in unbiased paths.
// Each outgoing/incoming edge is chosen with equal probability.
void Graph::fillConnectivity()
{
    Graph& graph = *this;

    BGL_FORALL_VERTICES(v, graph, Graph) {
        Vertex& vertex = graph[v];

        vertex.outEdges.clear();
        BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
            vertex.outEdges.push_back(e);
        }

        vertex.inEdges.clear();
        BGL_FORALL_INEDGES(v, e, graph, Graph) {
            vertex.inEdges.push_back(e);
        }
    }
}
#endif



// This version results in paths weighted by number of common reads.
// Each outgoing/incoming edge is chosen with probability
// proportional to the number of common oriented reads.
void Graph::fillConnectivity()
{
    Graph& graph = *this;

    BGL_FORALL_VERTICES(v, graph, Graph) {
        Vertex& vertex = graph[v];

        vertex.outEdges.clear();
        BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
            const uint64_t commonCount = graph[e].segmentPairInformation.commonCount;
            for(uint64_t i=0; i<commonCount; i++) {
                vertex.outEdges.push_back(e);
            }
        }

        vertex.inEdges.clear();
        BGL_FORALL_INEDGES(v, e, graph, Graph) {
            const uint64_t commonCount = graph[e].segmentPairInformation.commonCount;
            for(uint64_t i=0; i<commonCount; i++) {
                vertex.inEdges.push_back(e);
            }
        }
    }
}



void Graph::findPaths([[maybe_unused]] vector< vector<Segment> >& assemblyPaths)
{
    performanceLog << timestamp << "ReadFollowing1::Graph::findPaths begins." << endl;

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t pathCount = 100;

    const Graph& graph = *this;
    bool debug = false;

    // Random generator used to generate random paths.
    std::mt19937 randomGenerator;


    PathGraph pathGraph(assemblyGraph);
    std::map<Segment, PathGraph::vertex_descriptor> pathGraphVertexMap;

    // Each long Segment generates a PathGraphVertex.
    std::set<vertex_descriptor> longSegments;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Vertex& vertex = graph[v];
        if(vertex.length >= assemblyGraph.options.readFollowingSegmentLengthThreshold) {
            longSegments.insert(v);
            const PathGraph::vertex_descriptor u = boost::add_vertex({vertex.segment}, pathGraph);
            pathGraphVertexMap.insert({vertex.segment, u});
        }
    }
    cout << "The PathGraph has " << pathGraphVertexMap.size() <<
        " vertices, each corresponding to a long segment." << endl;



    // Loop over PathGraph vertices (that is, over long segments).
    vector<vertex_descriptor> path;
    BGL_FORALL_VERTICES(u0, pathGraph, PathGraph) {
        const Segment segment0 = pathGraph[u0].segment;
        const auto it0 = vertexMap.find(segment0);
        SHASTA2_ASSERT(it0 != vertexMap.end());
        const vertex_descriptor v0 = it0->second;

        // Loop over both directions.
        for(uint64_t direction=0; direction<2; direction++) {
            // debug = (assemblyGraph[segment0].id == 83738);
            if(debug) {
                cout << "Working on segment " << assemblyGraph[segment0].id <<
                    " direction " << direction << endl;
            }



            // Generate pathCount random paths starting at v0 and moving in this direction.
            for(uint64_t i=0; i<pathCount; i++) {
                findRandomPath(v0, direction, randomGenerator, path, longSegments);
                const vertex_descriptor v1 = (direction == 0) ? path.back() : path.front();
                const Segment segment1 = graph[v1].segment;

                if(debug) {
                    cout << "Found a path ending at " << assemblyGraph[segment1].id <<
                        " of length " << path.size() << endl;
                }

                // Discard a trivial path.
                SHASTA2_ASSERT(not path.empty());
                if(path.size() == 1) {
                    continue;
                }

                // If this is not a long segment (that is, it does not correspond
                // to a PathGraph vertex), discard this path.
                const auto it1 = pathGraphVertexMap.find(segment1);
                if(it1 == pathGraphVertexMap.end()) {
                    continue;
                }
                const PathGraph::vertex_descriptor u1 = it1->second;

                // Update the PathGraph with this path.

                // Find the PathGraph edge between u0 and u1, creating it if necessary.
                PathGraph::vertex_descriptor uu0 = u0;
                PathGraph::vertex_descriptor uu1 = u1;
                if(direction == 1) {
                    std::swap(uu0, uu1);
                }
                auto[e, edgeExists] = boost::edge(uu0, uu1, pathGraph);
                if(not edgeExists) {
                    tie(e, ignore) = boost::add_edge(uu0, uu1, pathGraph);
                }

                // Update the PathGraphEdge with this path.
                PathGraphEdge& pathGraphEdge = pathGraph[e];
                PathGraphEdge::Info& info = pathGraphEdge.infos[direction];
                if(info.pathCount == 0) {
                    info.pathCount = 1;
                    info.path = path;
                } else {
                    ++info.pathCount;
                    if(path.size() > info.path.size()) {
                        info.path = path;
                    }
                }
            }

        }
    }
    pathGraph.removeNonBestEdges();

    // This will throw if the PathGraph has cycles.
    // We need a more permanent soution.
    transitiveReductionAny(pathGraph);

    pathGraph.writeGraphviz();
    cout << "The PathGraph has " << num_vertices(pathGraph) <<
        " vertices and " << num_edges(pathGraph) << " vertices." << endl;

    // Find linear chains of vertices in the PathGraph.
    vector< vector<PathGraph::vertex_descriptor> > chains;
    findLinearVertexChains(pathGraph, chains);



    // Each linear vertex chain generates a Segment sequence,
    // that is, a sequence of Segments that
    // should be assembled into a single Segment.
    assemblyPaths.clear();
    for(const vector<PathGraph::vertex_descriptor>& chain: chains) {
        if(chain.size() < 2) {
            continue;
        }
        const Segment segment0 = pathGraph[chain.front()].segment;
        const Segment segment1 = pathGraph[chain.back()].segment;

        if(debug) {
            cout << "Found an assembly path that begins at " <<
                assemblyGraph[segment0].id << " and ends at " <<
                assemblyGraph[segment1].id << endl;
        }

        // Create a new assembly path.
        assemblyPaths.emplace_back();
        vector<Segment>& assemblyPath = assemblyPaths.back();
        for(uint64_t i1=1; i1<chain.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const PathGraph::vertex_descriptor u0 = chain[i0];
            const PathGraph::vertex_descriptor u1 = chain[i1];

            // Locate the PathGraph edge.
            PathGraph::edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = boost::edge(u0, u1, pathGraph);
            SHASTA2_ASSERT(edgeExists);
            const PathGraphEdge& pathGraphEdge = pathGraph[e];

            // Get the longest path on this edge.
            const vector<vertex_descriptor>& path = pathGraphEdge.longestPath();

            // Convert it to a sequence of segments.
            vector<Segment> pathSegments;
            for(const vertex_descriptor v: path) {
                pathSegments.push_back(graph[v].segment);
            }

            SHASTA2_ASSERT(pathGraph[u0].segment == pathSegments.front());
            SHASTA2_ASSERT(pathGraph[u1].segment == pathSegments.back());

            // The path for this vertex already contains the segments
            // corresponding to u0 and u1. So, to avoid duplications,
            // for each edge except the last we copy the path
            // without the last segment.
            // For the last edge we copy the entire path.
            auto end = pathSegments.end();
            if(i1 != chain.size() - 1) {
                --end;
            }
            copy(pathSegments.begin(), end, back_inserter(assemblyPath));
        }
    }
    performanceLog << timestamp << "ReadFollowing1::Graph::findPaths ends." << endl;
}



#if 0
// Find assembly paths.
// Note these are paths in the ReadFollowing::Graph but not in the AssemblyGraph.
void Graph::findPaths(vector< vector<Segment> >& assemblyPaths) const
{
    const Graph& graph = *this;
    const bool debug = true;


    // A graph to store the paths we find.
    // Each vertex corresponds to a long segment.
    // An edge u0->u1 contains a path that starts at segment(u0)
    // and ends at segment(u1). We only keep one path u0->u1,
    // even though in many/most cases we will find two,
    // one in each direction.
    class PathGraphVertex {
    public:
        Segment segment;
    };
    class PathGraphEdge {
    public:
        vector<Segment> path;
    };
    using PathGraph = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        PathGraphVertex,
        PathGraphEdge>;
    PathGraph pathGraph;
    std::map<Segment, PathGraph::vertex_descriptor> pathGraphVertexMap;

    // Each long Segment generates a PathGraphVertex.
    std::set<vertex_descriptor> longSegments;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Vertex& vertex = graph[v];
        if(vertex.length >= assemblyGraph.options.readFollowingSegmentLengthThreshold) {
            longSegments.insert(v);
            const PathGraph::vertex_descriptor u = boost::add_vertex({vertex.segment}, pathGraph);
            pathGraphVertexMap.insert({vertex.segment, u});
        }
    }



    // For each PathGraphVertex, compute a path in each direction,
    // always stopping when another long segment is encountered.
    // Each path generates a PathGraphEdge, as long as an edge between
    // the same two PathGraph vertices does not already exist.
    vector<vertex_descriptor> path;
    BGL_FORALL_VERTICES(u, pathGraph, PathGraph) {
        const Segment segment = pathGraph[u].segment;
        const auto it = vertexMap.find(segment);
        SHASTA2_ASSERT(it != vertexMap.end());
        const vertex_descriptor v = it->second;
        for(uint64_t direction=0; direction<2; direction++) {
            findPath(v, direction, path, longSegments);

            if(path.size() > 1) {

                // Find the first/last segment of the path.
                const vertex_descriptor v0 = path.front();
                const vertex_descriptor v1 = path.back();
                const Segment segment0 = graph[v0].segment;
                const Segment segment1 = graph[v1].segment;

                // Locate the corresponding PathGraph vertices.
                const auto it0 = pathGraphVertexMap.find(segment0);
                const auto it1 = pathGraphVertexMap.find(segment1);
                if((it0 != pathGraphVertexMap.end()) and (it1 != pathGraphVertexMap.end())) {
                    const PathGraph::vertex_descriptor u0 = it0->second;
                    const PathGraph::vertex_descriptor u1 = it1->second;

                    // Create an edge, as long as an edge between
                    // the same two PathGraph vertices does not already exist.
                    bool edgeExists = false;
                    tie(ignore, edgeExists) = boost::edge(u0, u1, pathGraph);
                    if(not edgeExists) {

                        // Ok, we are going to add a PathGraph edge.

                        edge_descriptor e;
                        tie(e, ignore) = boost::add_edge(u0, u1, pathGraph);

                        // Fill in the path of this PathGraphEdge.
                        PathGraphEdge& pathGraphEdge = pathGraph[e];
                        for(const vertex_descriptor v: path) {
                            pathGraphEdge.path.push_back(graph[v].segment);
                        }
                    }
                }
            }
        }
    }



    if(debug) {
    	ofstream dot("PathGraph.dot");
        dot << "digraph PathGraph {\n";

        BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        	const Segment segment = pathGraph[v].segment;
        	dot << assemblyGraph[segment].id << ";\n";
        }

        BGL_FORALL_EDGES(e, pathGraph, PathGraph) {
        	const PathGraph::vertex_descriptor v0 = source(e, pathGraph);
        	const PathGraph::vertex_descriptor v1 = target(e, pathGraph);
        	const Segment segment0 = pathGraph[v0].segment;
        	const Segment segment1 = pathGraph[v1].segment;
        	dot << assemblyGraph[segment0].id << "->" << assemblyGraph[segment1].id << ";\n";
        }

    	dot << "}\n";
    }



    // Find linear chains of vertices in the PathGraph.
    vector< vector<PathGraph::vertex_descriptor> > chains;
    findLinearVertexChains(pathGraph, chains);



    // Each linear vertex chain generates a Segment sequence,
    // that is, a sequence of Segments that
    // should be assembled into a single Segment.
    assemblyPaths.clear();
    for(const vector<PathGraph::vertex_descriptor>& chain: chains) {
    	if(chain.size() < 2) {
    		continue;
    	}
        const Segment segment0 = pathGraph[chain.front()].segment;
        const Segment segment1 = pathGraph[chain.back()].segment;

        if(debug) {
			cout << "Found an assembly path that begins at " <<
				assemblyGraph[segment0].id << " and ends at " <<
				assemblyGraph[segment1].id << endl;
        }

        // Create a new assembly path.
        assemblyPaths.emplace_back();
        vector<Segment>& assemblyPath = assemblyPaths.back();
        for(uint64_t i1=1; i1<chain.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const PathGraph::vertex_descriptor u0 = chain[i0];
            const PathGraph::vertex_descriptor u1 = chain[i1];

            // Locate the PathGraph edge.
            PathGraph::edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = boost::edge(u0, u1, pathGraph);
            SHASTA2_ASSERT(edgeExists);
            const PathGraphEdge& pathGraphEdge = pathGraph[e];
            const vector<Segment>& path = pathGraphEdge.path;

            SHASTA2_ASSERT(pathGraph[u0].segment == path.front());
            SHASTA2_ASSERT(pathGraph[u1].segment == path.back());

            // The path for this vertex already contains the segments
            // corresponding to u0 and u1. So, to avoid duplications,
            // for each edge except the last we copy the path
            // without the last segment.
            // For the last edge we copy the entire path.
            auto end = path.end();
            if(i1 != chain.size() - 1) {
                --end;
            }
            copy(path.begin(), end, back_inserter(assemblyPath));
        }
    }
}
#endif



void Graph::writePaths()
{
    vector< vector<Segment> > assemblyPaths;
    findPaths(assemblyPaths);

    ofstream csv("AssemblyPaths.csv");
    cout << "Found " << assemblyPaths.size() << " assembly paths." << endl;
    for(const vector<Segment>& assemblyPath: assemblyPaths) {
        cout << "Assembly path with " << assemblyPath.size() <<
        " segments beginning at " << assemblyGraph[assemblyPath.front()].id <<
        " and ending at " << assemblyGraph[assemblyPath.back()].id << endl;

        for(const Segment& segment: assemblyPath) {
            csv << assemblyGraph[segment].id << ",";
        }
        csv << "\n";
    }
}



PathGraph::PathGraph(const AssemblyGraph& assemblyGraph) :
    assemblyGraph(assemblyGraph)
{}



void PathGraph::writeGraphviz() const
{
    writeGraphviz("PathGraph.dot");
}



void PathGraph::writeGraphviz(const string& fileName) const
{
    ofstream dot(fileName);
    writeGraphviz(dot);
}



void PathGraph::writeGraphviz(ostream& dot) const
{
    const PathGraph& pathGraph = *this;
    dot << std::fixed << std::setprecision(2);

    dot << "digraph PathGraph {\n";



    // Vertices.
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        const Segment segment = pathGraph[v].segment;
        dot << assemblyGraph[segment].id <<
            " ["
            "label=\"" << assemblyGraph[segment].id <<
            // "\\n" << forwardPathCount(v) << "/" << backwardPathCount(v) <<
            "\""
            "]"
            ";\n";
    }



    // Edges.
    BGL_FORALL_EDGES(e, pathGraph, PathGraph) {
        // const PathGraphEdge& edge = pathGraph[e];

        const PathGraph::vertex_descriptor v0 = source(e, pathGraph);
        const PathGraph::vertex_descriptor v1 = target(e, pathGraph);

        const Segment segment0 = pathGraph[v0].segment;
        const Segment segment1 = pathGraph[v1].segment;

        const bool isBestOutEdge = (e == bestOutEdge(v0));
        const bool isBestInEdge  = (e == bestInEdge (v1));

        string color;
        if(isBestOutEdge and isBestInEdge) {
            color = "green";
        } else  if(isBestOutEdge or isBestInEdge) {
            color = "black";
        } else {
            color = "red";
        }

        dot <<
            assemblyGraph[segment0].id << "->" <<
            assemblyGraph[segment1].id <<
            " ["
            "label=\"" <<
            // edge.forwardPathCount() << "/" <<
            // edge.backwardPathCount() << "\\n" <<
            forwardPathFraction(e) << "/" <<
            backwardPathFraction(e) <<
            "\""
            " color=" << color <<
            "];\n";
    }



    dot << "}\n";

}



uint64_t PathGraph::forwardPathCount(vertex_descriptor v) const
{
    const PathGraph& pathGraph = *this;
    uint64_t count = 0;
    BGL_FORALL_OUTEDGES(v, e, pathGraph, PathGraph) {
        count += pathGraph[e].forwardPathCount();
    }
    return count;
}



uint64_t PathGraph::backwardPathCount(vertex_descriptor v) const
{
    const PathGraph& pathGraph = *this;
    uint64_t count = 0;
    BGL_FORALL_INEDGES(v, e, pathGraph, PathGraph) {
        count += pathGraph[e].backwardPathCount();
    }
    return count;
}



double PathGraph::forwardPathFraction(edge_descriptor e) const
{
    const PathGraph& pathGraph = *this;
    const vertex_descriptor v = source(e, pathGraph);
    return double(pathGraph[e].forwardPathCount()) / double(forwardPathCount(v));
}



double PathGraph::backwardPathFraction(edge_descriptor e) const
{
    const PathGraph& pathGraph = *this;
    const vertex_descriptor v = target(e, pathGraph);
    return double(pathGraph[e].backwardPathCount()) / double(backwardPathCount(v));
}



PathGraph::edge_descriptor PathGraph::bestInEdge(vertex_descriptor v) const
{
    const PathGraph& pathGraph = *this;
    SHASTA2_ASSERT(in_degree(v, pathGraph) > 0);

    edge_descriptor eBest = {0, 0, 0};
    uint64_t bestPathCount = 0;
    BGL_FORALL_INEDGES(v, e, pathGraph, PathGraph) {
        const uint64_t pathCount = pathGraph[e].backwardPathCount();
        if(pathCount > bestPathCount) {
            eBest = e;
            bestPathCount = pathCount;
        }
    }
    return eBest;
}



PathGraph::edge_descriptor PathGraph::bestOutEdge(vertex_descriptor v) const
{
    const PathGraph& pathGraph = *this;
    SHASTA2_ASSERT(out_degree(v, pathGraph) > 0);

    edge_descriptor eBest = {0, 0, 0};
    uint64_t bestPathCount = 0;
    BGL_FORALL_OUTEDGES(v, e, pathGraph, PathGraph) {
        const uint64_t pathCount = pathGraph[e].forwardPathCount();
        if(pathCount > bestPathCount) {
            eBest = e;
            bestPathCount = pathCount;
        }
    }
    return eBest;
}



void PathGraph::removeNonBestEdges()
{
    PathGraph& pathGraph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, pathGraph, PathGraph) {
        const PathGraph::vertex_descriptor v0 = source(e, pathGraph);
        const PathGraph::vertex_descriptor v1 = target(e, pathGraph);

        const bool isBestOutEdge = (e == bestOutEdge(v0));
        const bool isBestInEdge  = (e == bestInEdge (v1));

        if(not (isBestOutEdge or isBestInEdge)) {
            edgesToBeRemoved.push_back(e);
        }
    }

    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, pathGraph);
    }

}




