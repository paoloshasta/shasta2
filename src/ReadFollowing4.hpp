#pragma once

/*****************************************************************

Read following in the AssemblyGraph.

The SearchGraph has one vertex for each
Segment, regardless of length. They are used to find assembly paths
that start and end at long Segments and can use zero or more
short Segments in-between.

The ConnectGraph has one vertex for each long segment.
Its edges store connection information between long segments.
See below for mode details.

*****************************************************************/

// Shasta.
#include "AssemblyGraphBaseClass.hpp"
#include "invalid.hpp"
#include "MultithreadedObject.hpp"
#include "ReadId.hpp"

// Standard library.
#include "array.hpp"
#include "string.hpp"
#include "utility.hpp"



namespace shasta2 {

    namespace ReadFollowing4 {

        class ReadFollower;

        class SearchGraph;
        class SearchGraphVertex;
        class SearchGraphEdge;
        using SearchGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            SearchGraphVertex,
            SearchGraphEdge>;

        class ConnectGraph;
        class ConnectGraphVertex;
        class ConnectGraphEdge;
        using ConnectGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            ConnectGraphVertex,
            ConnectGraphEdge>;

        // A Segment is an edge of the AssemblyGraph.
        using Segment = AssemblyGraphBaseClass::edge_descriptor;

        // EXPOSE WHEN CODE STABILIZES.
        const double logPThreshold = 10.;    // dB
        const double a = 3.;                // dB
        const double b = 10.;               // dB
        const uint64_t transitiveReductionMaxDistance = 4;
    }

}



class shasta2::ReadFollowing4::SearchGraphVertex {
public:
    // A Segment is an edge of the AssemblyGraph.
    Segment segment;

    // The sequence length or estimated offset of this AssemblyGraph edge.
    uint64_t length = invalid<uint64_t>;

    // This is set for long vertices (length >= readFollowingSegmentLengthThreshold).
    bool isLong = false;

    SearchGraphVertex(Segment, uint64_t length, bool isLong);
};



class shasta2::ReadFollowing4::SearchGraphEdge {
public:
    uint64_t commonCount;
    uint64_t missingCount0;
    uint64_t missingCount1;

    SearchGraphEdge(
        uint64_t commonCount,
        uint64_t missingCount0,
        uint64_t missingCount1);

    // These are computed by the constructor from the three above.
    // logs are in dB.
    double logP;
    double weight;

    bool isReverseComplement(const SearchGraphEdge&) const;
};



class shasta2::ReadFollowing4::SearchGraph : public SearchGraphBaseClass {
public:

    // A map that gives the vertex_descriptor corresponding to each Segment.
    std::map<Segment, vertex_descriptor> vertexMap;

    // Create a vertex and update the vertexMap.
    void createVertex(Segment, uint64_t length, bool isLong);

    // Prune removes all vertices that are not accessible from long
    // vertices in both directions.
    void prune();

    // Checks that the SearchGraph is strand-symmetric.
    void check(const AssemblyGraph&) const;

    // The vertex index map is needed to compute shortest paths.
    // It must be created when no more changes will be made to the graph.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    void createVertexIndexMap();

    void findShortestPath(Segment, vector<Segment>&) const;

    void writeGraphviz(
        const AssemblyGraph&,
        const string& name) const;

};



class shasta2::ReadFollowing4::ConnectGraphVertex {
public:
    // A Segment is an edge of the AssemblyGraph.
    Segment segment;

    // The sequence length or estimated offset of this AssemblyGraph edge.
    uint64_t length = invalid<uint64_t>;

    ConnectGraphVertex(Segment, uint64_t length);
};




// A ConnectGraph edge between two Segments segments0 and segment1
// can be created for one or both of two reasons:
// - The SegmentPairInformation between the two segment
//   satisfies our criteria for minCommonCount
//   or one or more of the logP values.
// - We found a shortest path between segment0 and segment1
//   in the SearchGraph.
class shasta2::ReadFollowing4::ConnectGraphEdge {
public:

    // Direct connection between the source and target Segments
    // of this edge, without any intervening short Segments.
    class DirectConnectInformation {
    public:
        uint64_t commonCount = 0;
        uint64_t missingCount0 = 0;
        uint64_t missingCount1 = 0;

        // These are computed by the constructor from the three above.
        // logs are in dB.
        double logP = std::numeric_limits<double>::min();
        double logPForward = std::numeric_limits<double>::min();
        double logPBackward = std::numeric_limits<double>::min();
        double maxLogP() const;

        DirectConnectInformation(
            uint64_t commonCount,
            uint64_t missingCount0,
            uint64_t missingCount1);
        DirectConnectInformation() {}
        bool isReverseComplement(const DirectConnectInformation&) const;
    };
    DirectConnectInformation directConnectInformation;

    bool hasDirectConnection() const
    {
        return directConnectInformation.commonCount > 0;
    }

    enum class DirectConnectionType {
        None,
        Bidirectional,
        Forward,
        Backward,
        Ambiguous
    };
    DirectConnectionType directConnectionType() const;

    // Constructor that initializes the DirectConnectInformation.
    ConnectGraphEdge(
        uint64_t commonCount,
        uint64_t missingCount0,
        uint64_t missingCount1);

    ConnectGraphEdge() {}

    // The assembly paths between segment0 and segment1 found
    // for each direction.
    array<vector<Segment>, 2> assemblyPaths;

    bool isReverseComplement(const ConnectGraphEdge&) const;
};



class shasta2::ReadFollowing4::ConnectGraph : public ConnectGraphBaseClass {
public:

    // A map that gives the vertex_descriptor corresponding to each Segment.
    std::map<Segment, vertex_descriptor> vertexMap;

    // Create a vertex and update the vertexMap.
    void createVertex(Segment, uint64_t length);

    void writeGraphviz(
        const AssemblyGraph&,
        const string& name) const;

    // This removes edges without a direct connection
    // and that don't have paths in both directions.
    void removeWeakEdges();

    // Transitive reduction.
    void transitiveReduction();
    bool transitiveReductionCanRemove(edge_descriptor) const;

    vector<Segment> getAssemblyPath(edge_descriptor) const;

    // Check that it is strand-symmetric.
    void check(const AssemblyGraph&) const;
};



class shasta2::ReadFollowing4::ReadFollower :
    public MultithreadedObject<ReadFollower> {

public:
    ReadFollower(const AssemblyGraph&);

    // Use the SearchGraph to find a shortest path starting at segment0
    // and ending at a long Segment, with path length defined by SearchGraphEdge::weight.
    void findShortestPath(
        Segment segment0,
        uint64_t direction,     // 0 = forward, 1 = backward
        vector<Segment>& path
        ) const;
    void findShortestPathForward(
        Segment segment0,
        vector<Segment>& path
        ) const;
    void findShortestPathBackward(
        Segment segment0,
        vector<Segment>& path
        ) const;
    void findAndWriteShortestPath(Segment, uint64_t direction) const; // Python callable

    // Initial and final support for each Segment.
    std::map<Segment, vector<OrientedReadId> > initialSupportMap;
    std::map<Segment, vector<OrientedReadId> > finalSupportMap;
    void fillSupportMaps();

    // Segment pairs (segment0, segment1) such that the final support
    // of segment0 shares at least readFollowingMinCommonCount
    // OrientedReadIds with the initial support of segment1.
    // Each of them is a candidate edge for the SearchGraph.
    vector< pair<Segment, Segment> > segmentPairs;
    void findSegmentPairs();

    // The search graph used for shortest paths.
    SearchGraph searchGraph;

    // Also store a ConnectGraph, which contains only vertices corresponding to long Segments.
    // Information on intervening short segments is stored in the edges.
    ConnectGraph graph;



    // Use the SearchGraph to find shortest paths between long segments
    // and store them in the ConnectGraph.
    void findShortestPathsMultithreaded(uint64_t threadCount);
    void findShortestPathsThreadFunction(uint64_t threadId);
    class FindShortestPathsData {
    public:
        vector<ConnectGraph::vertex_descriptor> graphVertices;
    };
    FindShortestPathsData findShortestPathsData;



    void createVertices();
    void createEdges();
    void createEdgesThreadFunction(uint64_t threadId);

    const AssemblyGraph& assemblyGraph;

    // Use the ConnectGraph to update the AssemblyGraph.
    void updateAssemblyGraph(AssemblyGraph& writableAssemblyGraph) const;

private:
    bool isLong(Segment) const;

    // Make a disconnected copy of a Segment. The copy keeps the same id,
    // so the original Segment will have to be removed for the AssemblyGraph
    // to remain valid.
    Segment createDisconnectedSegmentCopy(
        AssemblyGraph& writableAssemblyGraph,
        Segment oldSegment) const;
};
