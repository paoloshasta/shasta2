#pragma once

// Local read following, one Tangle at a time.

// Shasta2.
#include "AssemblyGraphBaseClass.hpp"
#include "invalid.hpp"

// Standard library.
#include "fstream.hpp"
#include "iosfwd.hpp"
#include "string.hpp"
#include "vector.hpp"

namespace shasta2 {
    namespace ReadFollowing {

    class Graph;
    class Vertex;
    class Edge;

    using GraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        Vertex,
        Edge>;

    class CondensedGraphVertex;
    class CondensedGraphEdge;
    using CondensedGraph =
        boost::adjacency_list<boost::listS, boost::vecS, boost::bidirectionalS, CondensedGraphVertex, CondensedGraphEdge>;

    }

    class OrientedReadId;
    class Tangle;
}


// There is a vertex for each Segment internal to the Tangle,
// plus a vertex for each entrance and exit of the tangle.
// If a Segment is both an entrance and an exit
// there will be two vertices associated with it.
class shasta2::ReadFollowing::Vertex {
public:
    Segment segment;
    bool isEntrance;
    bool isExit;

    // The reverse complement of this vertex.
    // This is only filled in for a self-complementary tangle.
    GraphBaseClass::vertex_descriptor vRc = GraphBaseClass::null_vertex();

    vector<OrientedReadId> initialSupport;
    vector<OrientedReadId> finalSupport;

    Vertex(
        const AssemblyGraph&,
        Segment segment,
        bool isEntrance,
        bool isExit);
};




// An edge can be created for one or both of two reasons:
// - The two Segments have sufficiently similar read compositions,
//   as determined based on commonCount and logP.
// - The two Segments are adjacent in the AssemblyGraph.
class shasta2::ReadFollowing::Edge {
public:

    bool areAdjacent = false;
    bool haveSimilarReadCompositions = false;

    uint64_t commonCount = invalid<uint64_t>;
    double logP = invalid<double>;
    double weight = invalid<double>;

    Edge(uint64_t commonCount, double logP) :
        haveSimilarReadCompositions(true),
        commonCount(commonCount),
        logP(logP),
        weight(std::pow(10., -0.1 * logP))
    {}

    Edge() : areAdjacent(true) {}

};



class shasta2::ReadFollowing::Graph: public GraphBaseClass {
public:
    Graph(
        AssemblyGraph&,
        uint64_t tangleId,
        const Tangle&,
        ostream& html,
        const string& debugOutputBaseName);

    bool updateAssemblyGraph1();
    bool updateAssemblyGraph2();

private:

    // Constructor arguments.
    AssemblyGraph& assemblyGraph;
    uint64_t tangleId;
    const Tangle& tangle;
    ostream& html;
    string debugOutputBaseName;

    bool isSelfComplementaryTangle;

    vector<vertex_descriptor> entranceVertices;
    vector<vertex_descriptor> exitVertices;

    void createVertices();
    std::map<Segment, vertex_descriptor> vertexMap;
    void createReadSimilarityEdges();
    void createAdjacencyEdges();

    void partialTransitiveReduction();
    void disconnectBranchVertices();

    // This returns the number of adjacency edges that were removed.
    uint64_t cleanupEdges();

    // Remove inaccessible vertices.
    // This returns false if any entrances or exits
    // are inaccessible.
    bool removeInaccessibleVertices();

    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;
    void writeGraphToHtml(
        const string& title,
        const string& fileNameSuffix) const;

    // Shortest paths on the graph between each entrance/exit pair.
    class ShortestPath {
    public:
        Segment entrance;
        Segment exit;

        // The distance on the Graph between this entrance and this exit
        // (sum of weights). If the exit is not reachable from this entrance,
        // this is left to its initial invalid value.
        double distance = invalid<double>;
        bool isUnreachable() const
        {
            return distance == invalid<double>;
        }
        bool isReachable() const
        {
            return not isUnreachable();
        }

        // The path Segments, excluding the entrance and the exit.
        vector<Segment> segments;
    };
    vector< vector<ShortestPath> > shortestPaths;
    void findShortestPaths();
    void writeShortestPaths();

    // Return the number of exits reachable from an entrance.
    uint64_t countReachableExits(uint64_t iEntrance) const;

    // Return the number of entrances from which an exit is reachable.
    uint64_t countReachableEntrances(uint64_t iExit) const;

};



class shasta2::ReadFollowing::CondensedGraphVertex {
public:
    vector<Graph::vertex_descriptor> vertices;
};



class shasta2::ReadFollowing::CondensedGraphEdge {
public:
    Graph::edge_descriptor e = {0, 0, 0};
};
