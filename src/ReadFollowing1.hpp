#pragma once

// Read following in the AssemblyGraph.

// Shasta.
#include "AssemblyGraph.hpp"
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

    // These are needed for efficient creation of random paths.
    // We cannot use vecS for the adjacency_list because we
    // need to remove vertices and edges (for pruning).
    // These are filled by fillConnectivity after
    // pruning is done and the graph will no longer change.
    vector<GraphBaseClass::edge_descriptor> outEdges;
    vector<GraphBaseClass::edge_descriptor> inEdges;

};



class shasta2::ReadFollowing1::Edge {
public:
    Edge(const AssemblyGraph&, Segment, Segment);
    SegmentPairInformation segmentPairInformation;
};



class shasta2::ReadFollowing1::Graph : public GraphBaseClass {
public:
    Graph(const AssemblyGraph&);

private:
    const AssemblyGraph& assemblyGraph;

    // Initial creation.
    std::map<Segment, vertex_descriptor> vertexMap;
    void createVertices();
    void createEdges();

    // Prune short leaves.
    void prune();
    bool pruneIteration();

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
    // The path ends when one of the stop vertices is encountered,
    // but can end sooner.
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



    // Find assembly paths.
    void findPaths(vector< vector<Segment> >& assemblyPaths);

    // Python callable.
    void writeRandomPath(Segment, uint64_t direction);
    void writePaths();
};



// Each PathGraph vertex corresponds to a long segment.
class shasta2::ReadFollowing1::PathGraphVertex {
public:
    Segment segment;
};



class shasta2::ReadFollowing1::PathGraphEdge {
public:

    // Store information for each direction.
    class Info {
    public:

        // The number of paths between these two vertices found in each direction.
        uint64_t pathCount = 0;

        // The longest of the paths.
        // These are paths in the ReadFollowing1::Graph
        // but not in the AssemblyGraph.
        vector<Graph::vertex_descriptor> path;
    };

    array<Info, 2> infos;

    uint64_t forwardPathCount() const
    {
        return infos[0].pathCount;
    }
    uint64_t backwardPathCount() const
    {
        return infos[1].pathCount;
    }

    const vector<Graph::vertex_descriptor>& longestPath() const
    {
        return
            (infos[0].path.size() >= infos[1].path.size()) ?
            infos[0].path :
            infos[1].path;

    }
};



// Class used to store paths between long segments.
class shasta2::ReadFollowing1::PathGraph : public PathGraphBaseClass {
public:
    PathGraph(const AssemblyGraph&);
    void removeNonBestEdges();

    // Graphviz output.
    void writeGraphviz() const;
    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;

private:
    const AssemblyGraph& assemblyGraph;

    uint64_t forwardPathCount(vertex_descriptor) const;
    uint64_t backwardPathCount(vertex_descriptor) const;

    double forwardPathFraction(edge_descriptor) const;
    double backwardPathFraction(edge_descriptor) const;

    edge_descriptor bestInEdge(vertex_descriptor) const;
    edge_descriptor bestOutEdge(vertex_descriptor) const;

};
