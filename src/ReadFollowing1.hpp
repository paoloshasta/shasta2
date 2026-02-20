#pragma once

// Read following in the AssemblyGraph.

// Shasta.
#include "AssemblyGraph.hpp"
#include "MultithreadedObject.hpp"
#include "SegmentStepSupport.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <map>
#include <random>
#include <set>



namespace shasta2 {

    namespace ReadFollowing1 {

        class Graph;
        class Vertex;
        class Edge;
        using GraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            Vertex,
            Edge>;

        using Path = vector<GraphBaseClass::vertex_descriptor>;

        class PathGraph;
        class PathGraphVertex;
        class PathGraphEdge;
        using PathGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            PathGraphVertex,
            PathGraphEdge>;

        // A Segment is an edge of the AssemblyGraph.
        using Segment = AssemblyGraph::edge_descriptor;
    }
}



class shasta2::ReadFollowing1::Vertex {
public:
    // A Segment is an edge of the AssemblyGraph.
    Segment segment;

    // The sequence length or estimated offset of this AssemblyGraph edge.
    uint64_t length = invalid<uint64_t>;

    double coverage = 0.;

    Vertex(const AssemblyGraph&, Segment);

    vector<SegmentStepSupport> initialSupport;
    vector<SegmentStepSupport> finalSupport;

    // Vectors used for the generation of random paths.
    // Each outgoing/incoming edge of this vertex appears
    // a number of tinmes equal to its commonCount.
    // This way, the random paths are biased by commonCount.
    // These are filled in by fillConnectivity.
    vector<GraphBaseClass::edge_descriptor> outEdges;
    vector<GraphBaseClass::edge_descriptor> inEdges;

};



class shasta2::ReadFollowing1::Edge {
public:
    Edge(const AssemblyGraph&, Segment, Segment);
    SegmentPairInformation segmentPairInformation;
};



class shasta2::ReadFollowing1::Graph :
    public GraphBaseClass,
    public MultithreadedObject<Graph> {
public:
    Graph(const AssemblyGraph&);

    // EXPOSE WHEN CODE STABILIZES.
    static const uint64_t pathCount = 100;
    static constexpr double minPathCountFraction = 0.4;

private:
    const AssemblyGraph& assemblyGraph;

    // Vertex creation.
    std::map<Segment, vertex_descriptor> vertexMap;
    void createVertices();
    std::set<vertex_descriptor> longVertices;

    // Create edge candidates.
    class EdgeCandidate {
    public:
        vertex_descriptor v0;
        vertex_descriptor v1;
        bool operator==(const EdgeCandidate& that) const
        {
            return tie(v0, v1) == tie(that.v0, that.v1);
        }
        bool operator<(const EdgeCandidate& that) const
        {
            return tie(v0, v1) < tie(that.v0, that.v1);
        }
    };
    vector<EdgeCandidate> edgeCandidates;
    void createEdgeCandidates();

    // Edge creation.
    void createEdges();
    void createEdgesMultithreaded();
    void createEdgesThreadFunction(uint64_t);

    // Prune removes all vertices that are not accessible from "long"
    // vertices in both directions.
    void prune();

    // This fills in the outEdges and inEdges vectors
    // of all vertices, which are then used to generate
    // random paths.
    void fillConnectivity();

    // Output.
    void write(const string& name) const;
    void writeCsv(const string& name) const;
    void writeVerticesCsv(const string& name) const;
    void writeEdgesCsv(const string& name) const;
    void writeGraphviz(const string& name) const;

public:



    // Random paths.
    // These functions find a random path starting at the given vertex.
    // Direction is 0 for forward and 1 backward.
    // Note these are paths in the ReadFollowing1::Graph
    // but not in the AssemblyGraph.
    template<std::uniform_random_bit_generator RandomGenerator> void findRandomPath(
        vertex_descriptor, uint64_t direction,
        RandomGenerator&,
        vector<vertex_descriptor>& path,
        const std::set<vertex_descriptor>& stopVertices);
    template<std::uniform_random_bit_generator RandomGenerator> void findRandomForwardPath(
        vertex_descriptor,
        RandomGenerator&,
        vector<vertex_descriptor>& path,
        const std::set<vertex_descriptor>& stopVertices);
    template<std::uniform_random_bit_generator RandomGenerator> void findRandomBackwardPath(
        vertex_descriptor,
        RandomGenerator&,
        vector<vertex_descriptor>& path,
        const std::set<vertex_descriptor>& stopVertices);
    // Python callable.
    void writeRandomPath(Segment, uint64_t direction);

    // Non-random paths.
    void findPath(
        vertex_descriptor,
        uint64_t direction,
        vector<vertex_descriptor>& path) const;
    void findForwardPath(
        vertex_descriptor,
        vector<vertex_descriptor>& path) const;
    void findBackwardPath(
        vertex_descriptor,
        vector<vertex_descriptor>& path) const;
    // Python callable.
    void writePath(Segment, uint64_t direction);

    // Assembly paths.
    shared_ptr<PathGraph> createPathGraph();
    void findPaths(vector< vector<Segment> >& assemblyPaths);
    void writePaths(const vector< vector<Segment> >& assemblyPaths);
    void findAndWritePaths();

};



// Each PathGraph vertex corresponds to a long segment.
class shasta2::ReadFollowing1::PathGraphVertex {
public:
    Segment segment;

    // The number of paths starting at this segment, in each direction.
    array<uint64_t, 2> pathCount = {0, 0};
};



class shasta2::ReadFollowing1::PathGraphEdge {
public:

    // The Paths between v0 and v1 of this edge found in each direction.
    array< vector<Path>, 2> paths;

};



// Class used to store paths between long segments.
class shasta2::ReadFollowing1::PathGraph : public PathGraphBaseClass {
public:
    PathGraph(const AssemblyGraph&);
    void removeLowPathCountFractionEdges();
    void removeNonBestEdges();
    void cleanupBubbles();

    // This returns all the non-trivial connected components
    // of the PathGraph (that is, the ones with at least two vertices).
    vector< shared_ptr<PathGraph> > findConnectedComponents();

    // This computes an assembly path, assuming it is working
    // on a PathGraph with a single connected component.
    void findAssemblyPath(vector<Segment>&);

    // Graphviz output.
    void writeGraphviz(const string& name) const;
    void writeGraphviz(ostream&) const;

private:
    const AssemblyGraph& assemblyGraph;

    double pathCountFraction(edge_descriptor, uint64_t direction) const;
    double forwardPathCountFraction(edge_descriptor) const;
    double backwardPathCountFraction(edge_descriptor) const;

    edge_descriptor bestInEdge(vertex_descriptor) const;
    edge_descriptor bestOutEdge(vertex_descriptor) const;

    using Bubble = vector< vector<edge_descriptor> >;
    void findBubbles(vector<Bubble>&) const;
};
