#pragma once


// Shasta.
#include "AnchorPair.hpp"
#include "AssemblerOptions.hpp"
#include "Base.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"
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



class shasta::AssemblyGraphEdge : public vector<AnchorId> {
public:
    uint64_t id = invalid<uint64_t>;

    // The length of an AssemblyGraphEdge is the estimated length of its sequence,
    // equal to the sum of the base offsets for adjacent pairs of AnchorIds in this edge.
    uint64_t length(const Anchors& anchors) const;

    AnchorId second() const
    {
        SHASTA_ASSERT(size() > 1);
        return (*this)[1];
    }

    AnchorId secondToLast() const
    {
        SHASTA_ASSERT(size() > 1);
        return (*this)[size() - 2];
    }

    // Assembled sequence.
    // We store the assembled sequencse between each of the size()-1 pairs
    // of adjacent AnchorIds of this AssemblyGraphEdge.
    // Assembled sequence for this edge is the concatenation of these sequences.
    bool wasAssembled = false;  // Set if all steps of this edge have been assembled.
    vector< vector<Base> > sequences;
    uint64_t sequenceLength() const;
    void getSequence(vector<Base>&) const;

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object< vector<AnchorId> >(*this);
        ar & id;
        ar & wasAssembled;
        ar & sequences;
    }
};



class shasta::AssemblyGraph:
    public AssemblyGraphBaseClass,
    public MappedMemoryOwner,
    public MultithreadedObject<AssemblyGraph> {
public:

    // Initial construction from the AnchorGraph.
    AssemblyGraph(
        const AssemblerOptions&,
        const Anchors&,
        const AnchorGraph&,
        uint64_t threadCount);

    // Deserialize.
    AssemblyGraph(
        const AssemblerOptions&,
        const Anchors&,
        const string& assemblyStage);


    edge_descriptor addEdge(vertex_descriptor v0, vertex_descriptor v1)
    {
        edge_descriptor e;
        tie(e, ignore) = boost::add_edge(v0, v1, *this);
        (*this)[e].id = nextEdgeId++;
        return e;
    }

    const AssemblerOptions& assemblerOptions;
    const Anchors& anchors;

    void detangleVertices(Detangler&);
    void detangleEdges(Detangler&);

    void transitiveReduction(
        uint64_t threshold,
        uint64_t a,
        uint64_t b);

    // For each set of parallel edges with identical AnchorId sequences,
    // keep only one.
    void cleanupTrivialBubbles();

private:
    uint64_t nextEdgeId = 0;

    void compress();

    void write(const string& name) const;
    void writeGfa(const string& fileName) const;
    void writeFasta(const string& fileName) const;
    void writeGraphviz(const string& fileName) const;
    void writeSegments(const string& fileName) const;
    void writeSegmentDetails(const string& fileName) const;

    // For a given edge, this returns the minimum common count
    // for pairs of adjacent anchors in the edge.
    uint64_t minCommonCountOnEdge(edge_descriptor) const;

    // Same, but only counting journey offsets equal to 1.
    uint64_t minCommonCountOnEdgeAdjacent(edge_descriptor) const;


    using TangleTemplate = boost::adjacency_list<boost::listS, boost::vecS, boost::bidirectionalS>;
    static void writeGraphviz(ostream&, const TangleTemplate&);
    static TangleTemplate reverse(const TangleTemplate&);
    vector<TangleTemplate> tangleTemplates;
    void createTangleTemplates();
    void detangle(const TangleTemplate&, Detangler&);



    // Sequence assembly.

    // Assemble sequence for all edges.
    void assembleAll(uint64_t threadCount);

    // Assemble sequence for the specified edge.
    void assemble(edge_descriptor, uint64_t threadCount);

    // Assemble sequence for step i of the specified edge.
    // The sequences vector for the edge must have already been sized to the correct length.
    void assembleStep(edge_descriptor, uint64_t i);

    // Assemble sequence for all edges in the edgesToBeAssembled vector.
    // This fills in the edgeStepsToBeAssembled with all steps of those edges,
    // then assembles each of the steps in parallel.
    void assemble(uint64_t threadCount);
    void assembleThreadFunction(uint64_t threadId);
    vector<edge_descriptor> edgesToBeAssembled;
    vector< pair<edge_descriptor, uint64_t> > edgeStepsToBeAssembled;


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
