#pragma once

// Shasta2.
#include "AssemblyGraph.hpp"
#include "invalid.hpp"
#include "ReadId.hpp"

// Standard library.
#include "array.hpp"
#include "fstream.hpp"
#include <map>
#include "string.hpp"
#include "tuple.hpp"
#include "vector.hpp"

namespace shasta2 {
    namespace StrandSeparation1 {
        class ReadOccurrence;
        class SegmentInfo;
        class StrandContact;
        class SegmentPair;

        class BipartiteGraphVertex;
        class BipartiteGraphEdge;
        using BipartiteGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::vecS,
            boost::undirectedS,
            BipartiteGraphVertex,
            BipartiteGraphEdge>;
        class BipartiteGraph;
        class BipartiteGraphVertexStatistics;
        class BipartiteGraphSplit;

        class SegmentPairGraphVertex;
        class SegmentPairGraphEdge;
        using SegmentPairGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::vecS,
            boost::undirectedS,
            SegmentPairGraphVertex,
            SegmentPairGraphEdge>;
        class SegmentPairGraph;
        class SegmentGraphSplit;
 
        class SegmentGraphVertex;
        class SegmentGraphEdge;
        using SegmentGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::undirectedS,
            SegmentGraphVertex,
            SegmentGraphEdge>;
        class SegmentGraph;
 
        class ConnectGraphVertex;
        class ConnectGraphEdge;
        class ConnectGraph;
        using ConnectGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::undirectedS,
            ConnectGraphVertex,
            ConnectGraphEdge>;
   }
}



// Information about a single Segment in a StrandContact.
class shasta2::StrandSeparation1::SegmentInfo {
public:
    Segment segment;
    uint64_t id;
    bool isEntrance = false;
    bool isExit = false;
    BipartiteGraphBaseClass::vertex_descriptor v = BipartiteGraphBaseClass::null_vertex();
    uint64_t componentId = invalid<uint64_t>;
    ConnectGraphBaseClass::vertex_descriptor cv = ConnectGraphBaseClass::null_vertex();
};



// A pair of reverse complemented Segments in a StrandContact.
// The Segment with the lower id is stored first.
// Ordered by the id of the first segment in the pair.
class shasta2::StrandSeparation1::SegmentPair {
public:
    array<SegmentInfo, 2> segmentInfos;
    uint64_t length;
    double coverage;
    double crossStrandEdgeFrequencyRatio = invalid<double>;

    uint64_t id0() const
    {
        return segmentInfos[0].id;
    }
};



// Class used to count occurrences of reads in Segments.
class shasta2::StrandSeparation1::ReadOccurrence {
public:
    uint64_t segmentPairId = invalid<uint64_t>;
    Strand strand = invalid<Strand>;
    uint64_t frequency = 0;
    ReadOccurrence() {}
    ReadOccurrence(uint64_t segmentPairId, Strand strand) :
        segmentPairId(segmentPairId), strand(strand) {}

    bool operator==(const ReadOccurrence& that) const
    {
        return tie(segmentPairId, strand) == tie(that.segmentPairId, that.strand);
    }
    bool operator<(const ReadOccurrence& that) const
    {
        return tie(segmentPairId, strand) < tie(that.segmentPairId, that.strand);
    }
};



// A vertex of the BipartiteGraph can represent an OrientedReadId
// or a Segment.
class shasta2::StrandSeparation1::BipartiteGraphVertex {
public:
    bool isSegment;

    OrientedReadId orientedReadId = OrientedReadId(invalid<ReadId>, 0);

    uint64_t segmentPairId = invalid<uint64_t>;
    uint64_t segmentIndexInPair = invalid<uint64_t>;    // 0 or 1

    BipartiteGraphVertex(uint64_t segmentPairId, uint64_t segmentIndexInPair) :
        isSegment(true),
        segmentPairId(segmentPairId),
        segmentIndexInPair(segmentIndexInPair)
    {}

    BipartiteGraphVertex(OrientedReadId orientedReadId) :
        isSegment(false),
        orientedReadId(orientedReadId)
    {}

    BipartiteGraphVertex()
    {}
};



class shasta2::StrandSeparation1::BipartiteGraphEdge {
public:

    // The number of times (steps) this OrientedReadId appears
    // in this Segment.
    uint64_t frequency;
};



class shasta2::StrandSeparation1::BipartiteGraphVertexStatistics {
public:
    uint64_t totalEdgeCount = 0;
    uint64_t crossStrandEdgeCount = 0;
    uint64_t totalEdgeFrequency = 0;
    uint64_t crossStrandEdgeFrequency = 0;

    double crossStrandEdgeRatio() const;
    double crossStrandEdgeFrequencyRatio() const;
    uint64_t sameStrandEdgeCount() const;
    uint64_t sameStrandEdgeFrequency() const;
    int64_t crossStrandEdgeFrequencyExcess() const;
};



class shasta2::StrandSeparation1::BipartiteGraph : public BipartiteGraphBaseClass {
public:

    // Add a vertex representing a Segment.
    vertex_descriptor addVertex(uint64_t segmentPairId, uint64_t segmentIndexInPair);

    // Add a vertex representing an OrientedReadId.
    vertex_descriptor addVertex(OrientedReadId);

    // Return the reverse complement of a vertex.
    // Because vertices are added in reverse complemented pairs,
    // pairs of reverse complemented vertices have consecutive vertex_descriptors.
    vertex_descriptor reverseComplement(vertex_descriptor) const;

    // The pairs of reverse complemented edges.
    // Sorted by decreasing frequency.
    class EdgePair {
    public:
        edge_descriptor e;
        edge_descriptor eRc;
        uint64_t frequency;
    };
    vector<EdgePair> edgePairs;

    // This creates a BipartiteGraphSplit obtained by processing the EdgePairs
    // in the order described by the edgePairsIndexes.
    void computeSplit(
        const vector<uint64_t>& edgePairsIndexes,
        BipartiteGraphSplit&) const;

    BipartiteGraphVertexStatistics getVertexStatistics(
        const BipartiteGraphSplit&,
        vertex_descriptor) const;

    void writeGraphviz(
        const string& dotFileName,
        const vector<SegmentPair>&,
        const AssemblyGraph&,
        const BipartiteGraphSplit&) const;
};



// In the SegmentPairGraph, each vertex corresponds to a SegmentPair.
// The edges keep track of reads that support a pair of SegmentPairs
// in the same or opposite strands.
// The vertex_descriptor of the SegmentPairGraph is the segmentPairId
// (index in the segmentPairs vector).
class shasta2::StrandSeparation1::SegmentPairGraphVertex {
public:
};
class shasta2::StrandSeparation1::SegmentPairGraphEdge {
public:
    class Data {
    public:
        ReadId readId;
        bool isSameStrand;
        uint64_t frequencyProduct;
    };
    vector<Data> data;
    uint64_t coverage() const;
    uint64_t sameStrandCoverage() const;
    uint64_t oppositeStrandCoverage() const;

    // Orientation is +1 is all Data::isSameStrand are true,
    // -1 is all Data::isSameStrand are false,
    double orientation() const;
};
class shasta2::StrandSeparation1::SegmentPairGraph : public SegmentPairGraphBaseClass {
public:
    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;
    void writeCsv(const string& fileName, const AssemblyGraph&, const vector<SegmentPair>&) const;
    void writeCsv(ostream&, const AssemblyGraph&, const vector<SegmentPair>&) const;
};



// Each vertex of the SegmentPairGraph generates two vertices of the SegmentGraph.
class shasta2::StrandSeparation1::SegmentGraphVertex {
public:
    uint64_t segmentPairId;
    uint64_t segmentIndexInPair;
};
class shasta2::StrandSeparation1::SegmentGraphEdge {
public:
    uint64_t coverage;
};
class shasta2::StrandSeparation1::SegmentGraph : public SegmentGraphBaseClass {
public:

    // The vertexTable is indexed by [segmentPairId][segmentIndexInPair]
    vector< array<vertex_descriptor, 2> > vertexTable;

    // Return the vertex index.
    uint64_t vertexIndex(vertex_descriptor v) const
    {
        const auto& vertex = (*this)[v];
        return 2 * vertex.segmentPairId + vertex.segmentIndexInPair;
    }

    vertex_descriptor reverseComplement(vertex_descriptor v) const
    {
        const auto& vertex = (*this)[v];
        return vertexTable[vertex.segmentPairId][vertex.segmentIndexInPair ^ 1];
    }

    // The pairs of reverse complemented edges.
    // Sorted by decreasing frequency.
    class EdgePair {
    public:
        edge_descriptor e;
        edge_descriptor eRc;
        uint64_t coverage;
    };
    vector<EdgePair> edgePairs;

    void computeSplit(
        const vector<uint64_t>& edgePairsIndexes,
        SegmentGraphSplit&) const;

    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;
};



class shasta2::StrandSeparation1::SegmentGraphSplit {
public:
    vector< vector<SegmentGraph::vertex_descriptor> > components;

    // This gives the componentId (index in the components vector)
    // that each vertex belongs to.
    std::map<SegmentGraphBaseClass::vertex_descriptor, uint64_t> vertexComponent;

    // The cross-strand edge pairs that generated this SegmentGraphSplit.
    // These are indexes in the SegmentGraph::edgePairs vector.
    // They are sorted so we can do binary searches in them.
    vector<uint64_t> crossStrandEdgePairIndexes;
    bool isCrossStrandEdgePair(uint64_t) const;

    // The sum of the coverages of the crossStrandEdgePairs.
    // An optimal SegmentGraphSplit minimizes this.
    uint64_t crossStrandCoverage = 0;

    void clear();
};




// A possible way to separate the BipartiteGraph in two strands.
// Here, the connected components are computed without using
// the cross-strand edges. The Segment and each OrientedReadId
// are guaranteed not to be in the same component.
// If the BipartiteGraph has n components, a BipartiteGraphSplit will have 2*n.
// If the BipartiteGraph is connected, a BipartiteGraphSplit will have two components.
// The components are ordered by decreasing side, and with
// reverse complemented components consecutively numbered.
// Each component corresponds to a strand.
// The vertex_descriptors in each component are sorted.
// The BipartiteGraphSplit of a BipartiteGraph is not unique.
// An optimal BipartiteGraphSplit minimizes the sum of the frequencies
// of the cross-strand edges.
class shasta2::StrandSeparation1::BipartiteGraphSplit {
public:
    vector< vector<BipartiteGraph::vertex_descriptor> > components;

    // This gives the componentId (index in the components vector)
    // that each vertex belongs to. It is indexes by the
    // vertex_descriptor, which for the BipartiteGraph is simply uint64_t.
    vector<uint64_t> vertexComponent;

    // The cross-strand edge pairs that generated this BipartiteGraphSplit.
    // These are indexes in the BipartiteGraph::edgePairs vector.
    // They are sorted so we can do binary searches in them.
    vector<uint64_t> crossStrandEdgePairIndexes;
    bool isCrossStrandEdgePair(uint64_t) const;

    // The sum of the frequencies of the crossStrandEdgePairs.
    // An optimal BipartiteGraphSplit minimizes this.
    uint64_t crossStrandFrequency = 0;

    void clear();
};



class shasta2::StrandSeparation1::ConnectGraphVertex {
public:
    uint64_t segmentPairId = invalid<uint64_t>;
    uint64_t segmentIndexInPair = invalid<uint64_t>;
    uint64_t componentId = invalid<uint64_t>;
    bool isAmbiguous = false;

    ConnectGraphVertex() {}
    ConnectGraphVertex(
        uint64_t segmentPairId,
        uint64_t segmentIndexInPair,
        uint64_t componentId, 
        bool isAmbiguous) :
        segmentPairId(segmentPairId),
        segmentIndexInPair(segmentIndexInPair),
        componentId(componentId),
        isAmbiguous(isAmbiguous)
    {}
};



class shasta2::StrandSeparation1::ConnectGraphEdge {
public:
};



class shasta2::StrandSeparation1::ConnectGraph : public shasta2::StrandSeparation1::ConnectGraphBaseClass {
public:
    std::map<Segment, vertex_descriptor> vertexMap;
    void keepStrand0();
    void writeGraphviz(
        const string& dotFileName, 
        const vector<SegmentPair>&, 
        const AssemblyGraph&) const;
};



class shasta2::StrandSeparation1::StrandContact {
public:

    // EXPOSE WHEN CODE STABILIZES.
    const double maxCrossStrandFrequencyRatio = 0.2;
    const double segmentGraphOrientationThreshold = 0.8;
    const uint64_t segmentGraphCoverageThreshold = 10;

    // StrandContact constructor.
    // The strandContactVerticesmust be sorted by id.
    // The debugOutputBaseName and strandContactId are only used for debug output.
    StrandContact(
        AssemblyGraph&,
        const vector<AssemblyGraphBaseClass::vertex_descriptor>& strandContactVertices,
        const string& debugOutputBaseName,
        uint64_t strandContactId
        );

private:

    // Constructor arguments.
    AssemblyGraph& assemblyGraph;
    const vector<AssemblyGraphBaseClass::vertex_descriptor>& strandContactVertices;
    const string& debugOutputBaseName;
    uint64_t strandContactId;

    // The html is defined mutable to allow more functions to be const.
    mutable ofstream html;

    // Pairs of reverse complemented Segments in the StrandContact.
    // This includes Segments internal to the StrandContact plus
    // entrances and exits.
    // Ordered by the id of the first segment in the pair.
    // An index into this vector is a segmentPairId.
    vector<SegmentPair> segmentPairs;
    void gatherSegmentPairs();
    void writeSegmentPairs();

    uint64_t id(AssemblyGraphBaseClass::vertex_descriptor) const;
    uint64_t id(Segment) const;

    // Count occurrences of reads in the first Segment of SegmentPair.
    std::map<ReadId, vector<ReadOccurrence> > readOccurrenceMap;
    void countReadOccurrences();
    void writeReadOccurrences();

    // The BipartiteGraph has a vertex for Segment
    // and a vertex for each OrientedReadId that appears
    // in more than one Segment.
    // An edge (Segment--OrientedReadId) contains the number of times (steps)
    // that OrientedReadId occures in that Segment.
    // The vertex corresponding to a Segment is stored in the SegmentInfo.
    // The vertex corresponding to an OrientedReadId is stored in the
    // orientedReadIdVertexMap.
    BipartiteGraph bipartiteGraph;
    std::map<OrientedReadId, BipartiteGraph::vertex_descriptor> orientedReadIdVertexMap;
    void createBipartiteGraph();
    void writeBipartiteGraphGraphviz(const BipartiteGraphSplit&);
    void writeBipartiteGraphCustom(const BipartiteGraphSplit&);
    void computeSplit(BipartiteGraphSplit&) const;
    void writeSplitSummary(const BipartiteGraphSplit&) const;
    void writeSplitDetails(const BipartiteGraphSplit&) const;

    // Store component information in the SegmentPairs.
    void storeSegmentInformation(const BipartiteGraphSplit&);

    // SegmentPairGraph.
    void createSegmentPairGraph();
    SegmentPairGraph segmentPairGraph;

    // The SegmentGraph is created using the SegmentPairGraph.
    void createSegmentGraph();
    SegmentGraph segmentGraph;
    void computeSplit(SegmentGraphSplit&) const;

    // The ConnectGraph describes trhe connections between Segments that
    // exist in the AssemblyGraph and/or should be removed
    // for strand separation.
    ConnectGraph connectGraph;
    void createConnectGraph();
    void writeConnectGraph(const string& stage);
 
    // Initial version.
    void updateAssemblyGraph();

    // Version that uses the ConnectGraph.
    void updateAssemblyGraph1();

public:
    static string componentColor(uint64_t componentId);
    static string ambiguousColor();
};
