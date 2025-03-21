#pragma once


// Shasta.
#include "AnchorPair.hpp"
#include "AssemblerOptions.hpp"
#include "MappedMemoryOwner.hpp"
#include "shastaTypes.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/serialization/vector.hpp>

// Standard library.
#include <tuple.hpp>

namespace shasta {

    class AssemblyGraph;
    class AssemblyGraphVertex;
    class AssemblyGraphEdge;
    class AssemblyGraphVertexOrderByAnchorId;
    class AssemblyGraphVertexCompareEqualByAnchorId;
    class AssemblyGraphEdgeOrderById;
    class AssemblyGraphEdgeCompareEqualById;
    class AssemblyGraphStep;

    using AssemblyGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        AssemblyGraphVertex,
        AssemblyGraphEdge>;

    class Anchor;
    class AnchorGraph;
    class Detangler;

}



class shasta::AssemblyGraphVertex {
public:
    AnchorId anchorId;

    AssemblyGraphVertex(AnchorId anchorId = invalid<AnchorId>) : anchorId(anchorId) {}

    // This is used for the BFS in AssemblyGraph::transitiveReduction.
    uint64_t bfsDistance = invalid<uint64_t>;

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorId;
    }
};



// For now this is just an AnchorPair.
class shasta::AssemblyGraphStep : public AnchorPair {
public:

    AssemblyGraphStep(const AnchorPair& anchorPair) :
        AnchorPair(anchorPair) {}

    AssemblyGraphStep(
        const AnchorPair& anchorPair,
        const vector<OrientedReadId>& excludedOrientedReadIds) :
        AnchorPair(anchorPair, excludedOrientedReadIds) {}

    // It would be good to remove the default constructor but it is currently
    // used in detangling.
    AssemblyGraphStep() = default;

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object<AnchorPair>(*this);
    }
};



class shasta::AssemblyGraphEdge : public vector<AssemblyGraphStep> {
public:
    uint64_t id = invalid<uint64_t>;

    // Compute the sum of offsets of all the AssemblyGraphSteps.
    uint32_t getAverageOffset(const Anchors&) const;
    void getOffsets(
        const Anchors&,
        uint32_t& averageOffset,
        uint32_t& minOffset,
        uint32_t& maxOffset) const;

    // The length of an AssemblyGraphEdge is the estimated length of its sequence,
    // equal to averageOffset, which is the sum of the averageOffsets
    // of all the AssemblyGraphSteps of this edge.
    uint64_t length(const Anchors& anchors) const
    {
        uint32_t averageOffset;
        uint32_t minOffset;
        uint32_t maxOffset;
        getOffsets(anchors, averageOffset, minOffset, maxOffset);
        return averageOffset;
    }

    AnchorId anchorIdA() const
    {
        return (*this).front().anchorIdA;
    }
    AnchorId anchorIdB() const
    {
        return (*this).back().anchorIdB;
    }

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object< vector<AssemblyGraphStep> >(*this);
        ar & id;
    }
};



class shasta::AssemblyGraph:
    public AssemblyGraphBaseClass,
    public MappedMemoryOwner {
public:

    // Initial construction from the AnchorGraph.
    AssemblyGraph(
        const Anchors&,
        const AnchorGraph&,
        const AssemblerOptions::AssemblyGraphOptions&);

    // Deserialize.
    AssemblyGraph(
        const Anchors&,
        const string& assemblyStage);


    edge_descriptor addEdge(vertex_descriptor v0, vertex_descriptor v1)
    {
        edge_descriptor e;
        tie(e, ignore) = boost::add_edge(v0, v1, *this);
        (*this)[e].id = nextEdgeId++;
        return e;
    }

    const Anchors& anchors;

    void detangleVertices(Detangler&);
    void detangleEdges(Detangler&);

private:
    uint64_t nextEdgeId = 0;

    void transitiveReduction(
        uint64_t threshold,
        uint64_t a,
        uint64_t b);
    void compress();

    void write(const string& name) const;
    void writeGfa(const string& fileName) const;
    void writeGraphviz(const string& fileName) const;
    void writeSegments(const string& fileName) const;
    void writeSegmentDetails(const string& fileName) const;



    // Serialization.
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object<AssemblyGraphBaseClass>(*this);
        ar & nextEdgeId;
    }
    void save(ostream&) const;
    void load(istream&);

    // These do save/load to/from mapped memory.
    // The file name is AssemblyGraph-Stage.
    void save(const string& stage) const;
    void load(const string& stage);
};



// A function object class that orders AssemblyGraph::vertex_descriptors
// considering only their AnchorId.
class shasta::AssemblyGraphVertexOrderByAnchorId {
public:
    AssemblyGraphVertexOrderByAnchorId(const AssemblyGraph& assemblyGraph) : assemblyGraph(assemblyGraph) {}
    using vertex_descriptor = AssemblyGraph::vertex_descriptor;
    bool operator()(const vertex_descriptor& x, const vertex_descriptor& y) const
    {
        const uint64_t xAnchorId = assemblyGraph[x].anchorId;
        const uint64_t yAnchorId = assemblyGraph[y].anchorId;
        return xAnchorId < yAnchorId;
    }
private:
    const AssemblyGraph& assemblyGraph;
};



// A function object class that orders AssemblyGraph::vertex_descriptors
// considering only their AnchorId.
class shasta::AssemblyGraphVertexCompareEqualByAnchorId {
public:
    AssemblyGraphVertexCompareEqualByAnchorId(const AssemblyGraph& assemblyGraph) : assemblyGraph(assemblyGraph) {}
    using vertex_descriptor = AssemblyGraph::vertex_descriptor;
    bool operator()(const vertex_descriptor& x, const vertex_descriptor& y) const
    {
        const uint64_t xAnchorId = assemblyGraph[x].anchorId;
        const uint64_t yAnchorId = assemblyGraph[y].anchorId;
        return xAnchorId == yAnchorId;
    }
private:
    const AssemblyGraph& assemblyGraph;
};



// A function object class that orders AssemblyGraph::edge_descriptors
// considering only their id.
class shasta::AssemblyGraphEdgeOrderById {
public:
    AssemblyGraphEdgeOrderById(const AssemblyGraph& assemblyGraph) : assemblyGraph(assemblyGraph) {}
    using edge_descriptor = AssemblyGraph::edge_descriptor;
    bool operator()(const edge_descriptor& x, const edge_descriptor& y) const
    {
        const uint64_t xId = assemblyGraph[x].id;
        const uint64_t yId = assemblyGraph[y].id;
        return xId < yId;
    }
private:
    const AssemblyGraph& assemblyGraph;
};




class shasta::AssemblyGraphEdgeCompareEqualById {
public:
    AssemblyGraphEdgeCompareEqualById(const AssemblyGraph& assemblyGraph) : assemblyGraph(assemblyGraph) {}
    using edge_descriptor = AssemblyGraph::edge_descriptor;
    bool operator()(const edge_descriptor& x, const edge_descriptor& y) const
    {
        const uint64_t xId = assemblyGraph[x].id;
        const uint64_t yId = assemblyGraph[y].id;
        return xId == yId;
    }
private:
    const AssemblyGraph& assemblyGraph;
};
