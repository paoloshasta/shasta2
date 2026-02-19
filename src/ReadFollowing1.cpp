// Read following in the AssemblyGraph.

// Shasta.
#include "ReadFollowing1.hpp"
#include "deduplicate.hpp"
// #include "findLinearChains.hpp"
#include "Journeys.hpp"
// #include "Markers.hpp"
#include "Options.hpp"
// #include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
// #include "transitiveReduction.hpp"
using namespace shasta2;
using namespace ReadFollowing1;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"
#include <iomanip>
#include <random>

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<Graph>;



Graph::Graph(const AssemblyGraph& assemblyGraph) :
    MultithreadedObject<Graph>(*this),
    assemblyGraph(assemblyGraph)
{
    Graph& graph = *this;

    createVertices();
    createEdgeCandidates();
    createEdges();
    cout << "The read following graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;
  	write("Final");

    fillConnectivity();

}



// Create vertices of the ReadFollowing1 graph.
// Each vertex corresponds to a Segment of the AssemblyGraph.
void Graph::createVertices()
{
    Graph& graph = *this;

    BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {
        const vertex_descriptor v = add_vertex(Vertex(assemblyGraph, segment), graph);
        vertexMap.insert(make_pair(segment, v));
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



// Create edge candidates of the ReadFollowing1 graph.
// An edge candidate v0->v1 is created if the final support of v0
// shares at least one OrientedReadId with the initial support of v1.
void Graph::createEdgeCandidates()
{
    performanceLog << timestamp << "ReadFollowing1::Graph::createEdgeCandidates begins." << endl;

    Graph& graph = *this;
    const uint64_t orientedReadCount = assemblyGraph.journeys.size();
    const uint64_t minCommonCount = assemblyGraph.options.readFollowingMinCommonCount;

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
                edgeCandidates.push_back(EdgeCandidate({v0, v1}));
            }
        }
    }

    // Only keep the ones that appear at least minCommonCount times.
    vector<uint64_t> count;
    deduplicateAndCountWithThreshold(edgeCandidates, count, minCommonCount);
    SHASTA2_ASSERT(edgeCandidates.size() == count.size());

    cout << "Found " << edgeCandidates.size() << " candidate edges for the read following graph." << endl;
    performanceLog << timestamp << "ReadFollowing1::Graph::createEdgeCandidates ends." << endl;
}



// Generate an edge for each of edge candidate that satisfies our requirements.
void Graph::createEdges()
{
    performanceLog << timestamp << "ReadFollowing1::Graph::createEdges begins." << endl;
    Graph& graph = *this;
    const uint64_t minCommonCount = assemblyGraph.options.readFollowingMinCommonCount;
    const double minCorrectedJaccard = assemblyGraph.options.readFollowingMinCorrectedJaccard;
    const bool debug = false;

    for(uint64_t i=0; i<edgeCandidates.size(); i++) {
        const EdgeCandidate& edgeCandidate = edgeCandidates[i];

        const vertex_descriptor v0 = edgeCandidate.v0;
        const vertex_descriptor v1 = edgeCandidate.v1;

        const Segment segment0 = graph[v0].segment;
        const Segment segment1 = graph[v1].segment;

        // Create the candidate edge.
        const Edge edge(assemblyGraph, segment0, segment1);

        // This must be true given the way we constructed the vertex pairs.
        SHASTA2_ASSERT(edge.segmentPairInformation.commonCount >= minCommonCount);

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
                    " discarded due to low corrected Jaccard " << edge.segmentPairInformation.correctedJaccard << endl;
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
        dot << "label=\"" <<
            edge.segmentPairInformation.commonCount << "/" <<
            std::fixed << std::setprecision(2) <<
            edge.segmentPairInformation.correctedJaccard << "\"";

        // Tooltip.
        dot << " tooltip=\"" <<
            edge.segmentPairInformation.commonCount << "/" <<
            std::fixed << std::setprecision(2) <<
            edge.segmentPairInformation.correctedJaccard << "\"";

        // Thickness is proportional to commonCount.
        dot << " penwidth=" << 0.2 * double(edge.segmentPairInformation.commonCount);

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
    vector<vertex_descriptor>& path)
{
    if(direction == 0) {
        findRandomForwardPath(v, randomGenerator, path);
    } else {
        findRandomBackwardPath(v, randomGenerator, path);
    }
}



template<std::uniform_random_bit_generator RandomGenerator> void Graph::findRandomForwardPath(
    vertex_descriptor v,
    RandomGenerator& randomGenerator,
    vector<vertex_descriptor>& path)
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
    }
}



template<std::uniform_random_bit_generator RandomGenerator> void Graph::findRandomBackwardPath(
    vertex_descriptor v,
    RandomGenerator& randomGenerator,
    vector<vertex_descriptor>& path)
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
    std::random_device randomGenerator;
    findRandomPath(v, direction, randomGenerator, path);

    cout << "Found a path of length " << path.size() << ":" << endl;
    for(const vertex_descriptor v: path) {
        const Segment segment = graph[v].segment;
        cout << assemblyGraph[segment].id << ",";
    }
    cout << endl;
}


void Graph::writePath(Segment segment, uint64_t direction)
{
    const Graph& graph = *this;

    const auto it = vertexMap.find(segment);
    SHASTA2_ASSERT(it != vertexMap.end());
    const vertex_descriptor v = it->second;

    vector<vertex_descriptor> path;
    findPath(v, direction, path);

    cout << "Found a path of length " << path.size() << ":" << endl;
    for(const vertex_descriptor v: path) {
        const Segment segment = graph[v].segment;
        cout << assemblyGraph[segment].id << ",";
    }
    cout << endl;
}



void Graph::findPath(
    vertex_descriptor v,
    uint64_t direction,
    vector<vertex_descriptor>& path) const
{
    if(direction == 0) {
        findForwardPath(v, path);
    } else {
        findBackwardPath(v, path);
    }
}



void Graph::findForwardPath(
    vertex_descriptor v,
    vector<vertex_descriptor>& path) const
{
    const Graph& graph = *this;

    path.clear();
    path.push_back(v);

    vertex_descriptor v0 = v;
    while(true) {

        cout << "Finding next vertex for " << assemblyGraph[graph[v0].segment].id << endl;

        // Loop over all out-edges of v0.
        // Compute the sum of common counts against all vertices already in the path.
        // Keep the one with the greates sum.
        vertex_descriptor vBest = null_vertex();
        uint64_t sumBest = 0;
        BGL_FORALL_OUTEDGES(v0, e, graph, Graph) {
            const vertex_descriptor v1 = target(e, graph);

            // Loop over vertices already in the path.
            uint64_t sum = 0;
            for(const vertex_descriptor vPrevious: path) {
                auto[ePrevious, edgeExists] = boost::edge(vPrevious, v1, graph);
                if(edgeExists) {
                    sum += graph[ePrevious].segmentPairInformation.commonCount;
                }
            }
            cout << "Sum for " << assemblyGraph[graph[v1].segment].id << " is " << sum << endl;
            if(sum > sumBest) {
                sumBest = sum;
                vBest = v1;
            }
        }

        if(vBest== null_vertex()) {
            break;
        }
        cout << "Adding " << assemblyGraph[graph[vBest].segment].id << " to path." << endl;
        path.push_back(vBest);
        v0 = vBest;

    }
}



void Graph::findBackwardPath(
    vertex_descriptor,
    [[maybe_unused]] vector<vertex_descriptor>& path) const
{
    SHASTA2_ASSERT(0);
}



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


