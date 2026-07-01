#pragma once

// In the AnchorGraph, each vertex corresponds to an AnchorId
// and each edge corresponds to an AnchorPair.
// It uses boost::vecS as its second template argument,
// and as a result its vertex descriptors are AnchorIds.

// Shasta.
#include "Anchor.hpp"
#include "AnchorPair.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/serialization/base_object.hpp>

// Standard library.
#include "utility.hpp"
#include "vector.hpp"



namespace shasta2 {

    class AnchorGraph;
    class AnchorGraphEdge;
    using AnchorGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::bidirectionalS,
        boost::no_property,
        AnchorGraphEdge>;

    class Anchors;
    class Journeys;
    class ReadLengthDistribution;
}



class shasta2::AnchorGraphEdge {
public:
    AnchorId anchorIdA = invalid<AnchorId>;
    AnchorId anchorIdB = invalid<AnchorId>;

    // Begin/end indexes in the AnchorGraph::orientedReadIds vector
    // for the OrientedReadIds that belong to this edge.
    uint64_t orientedReadIdsBegin = invalid<uint64_t>;
    uint64_t orientedReadIdsEnd = invalid<uint64_t>;

    uint64_t id = invalid<uint64_t>;
    bool useForAssembly = false;

    AnchorGraphEdge(
        AnchorId anchorIdA,
        AnchorId anchorIdB,
        uint64_t orientedReadIdsBegin,
        uint64_t orientedReadIdsEnd,
        uint64_t id,
        bool useForAssembly) :
        anchorIdA(anchorIdA),
        anchorIdB(anchorIdB),
        orientedReadIdsBegin(orientedReadIdsBegin),
        orientedReadIdsEnd(orientedReadIdsEnd),
        id(id),
        useForAssembly(useForAssembly)
    {}

    AnchorGraphEdge() {}

    uint64_t coverage() const
    {
        return orientedReadIdsEnd - orientedReadIdsBegin;
    }

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorIdA;
        ar & anchorIdB;
        ar & orientedReadIdsBegin;
        ar & orientedReadIdsEnd;
        ar & id;
        ar & useForAssembly;
    }
};



class shasta2::AnchorGraph :
    public AnchorGraphBaseClass,
    public MappedMemoryOwner,
    public MultithreadedObject<AnchorGraph> {
public:

    // Construct the AnchorGraph from the Journeys.
    // Only include edges with at least the specified minCoverage.
    AnchorGraph(const Anchors&, const Journeys&, uint64_t minEdgeCoverage);

    // Constructor from binary data.
    AnchorGraph(const MappedMemoryOwner&, const string& name);

    AnchorPair getAnchorPair(edge_descriptor e) const
    {
        const AnchorGraphEdge& edge = (*this)[e];
        return AnchorPair(
            edge.anchorIdA,
            edge.anchorIdB,
            span<const OrientedReadId>(
                orientedReadIds.begin() + edge.orientedReadIdsBegin,
                orientedReadIds.begin() + edge.orientedReadIdsEnd));
    }

    // To reduce memory fragmentation, the OrientedReadIds of all edges
    // are stored together in this vector.
    // AnchorGraphEdge::orientedReadIdsBegin and AnchorGraphEdge::orientedReadIdsEnd
    // are indexes that identify the OrientedReadIds that belong to each edge.
    vector<OrientedReadId> orientedReadIds;

    uint64_t nextEdgeId = 0;
    edge_descriptor addEdge(
        AnchorId anchorIdA,
        AnchorId anchorIdB,
        const vector<OrientedReadId>& edgeOrientedReadIds,
        bool useForAssembly
        )
    {
        const uint64_t orientedReadIdsBegin = orientedReadIds.size();
        std::ranges::copy(edgeOrientedReadIds, back_inserter(orientedReadIds));
        const uint64_t orientedReadIdsEnd = orientedReadIds.size();
        auto[e, ignore] = add_edge(
            anchorIdA, anchorIdB,
            AnchorGraphEdge(anchorIdA, anchorIdB, orientedReadIdsBegin, orientedReadIdsEnd, nextEdgeId++, useForAssembly),
            *this);
        return e;
    }

    void transitiveReduction(
        uint64_t transitiveReductionMaxEdgeCoverage,
        uint64_t maxDistance);

    // Return the reverse complement of a vertex or edge.
    vertex_descriptor reverseComplement(vertex_descriptor) const;
    edge_descriptor reverseComplement(edge_descriptor) const;

private:
    bool transitiveReductionCanRemove(edge_descriptor, uint64_t transitiveReductionMaxDistance) const;
public:

    // Serialization.
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object<AnchorGraphBaseClass>(*this);
        ar & orientedReadIds;
    }
    void save(ostream&) const;
    void load(istream&);

    // These do save/load to/from mapped memory.
    void save(const string& name) const;
    void load(const string& name);
};

