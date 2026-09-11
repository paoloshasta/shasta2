#pragma once

// Shasta2.
#include "AssemblyGraphBaseClass.hpp"
#include "Tangle.hpp"

// Standard library.
#include "fstream.hpp"
#include "string.hpp"
#include "vector.hpp"

namespace shasta2 {
    class StrandSplitter;

    class AssemblyGraph;
}



// This takes as input a self-complementary tangle
// in an AssemblyGraph and attempts to split the strands.
class shasta2::StrandSplitter {
public:


    // The last two arguments are only used for debug output.
    StrandSplitter(
        AssemblyGraph&,
        const vector<AssemblyGraphBaseClass::vertex_descriptor>& tangleVertices,
        uint64_t tangleId,
        const string& debugOutputBaseName);

    ~StrandSplitter();

private:

    // EXPOSE WHEN CODE STABILIZES.
    const double maxCoverage = 16.;

    AssemblyGraph& assemblyGraph;
    uint64_t id(Segment) const;

    // If debug is set to true, debug output is written to html.
    bool debug = true;
    uint64_t tangleId;
    string debugOutputBaseName;
    ofstream html;
    void writeInitialDebugOutput();

    Tangle tangle;



    // The Tangle Segments.
    void gatherSegments();

    // All Tangle Segments (entrances, exits, and internal segments).
    // Sorted by id.
    vector<Segment> allTangleSegments;

    // The segmentPairs and segments vector only include Tangle Segments
    // with coverage no greater than maxCoverage. These are the ones
    // that are considered reliably single-copy and are used
    // for strand separation.

    // Pairs of reverse complemented segments
    // with coverage no greater than maxCoverage.
    vector< pair<Segment, Segment> > segmentPairs;
    void writeSegmentPairs();

    // Segments with coverage no greater than maxCoverage,
    // ordered by their appearance in segmentPairs.
    vector<Segment> segments;



    // Gather occurrences of reads in the first Segment of each pair.
    void findReadOccurrences();
    class ReadOccurrence {
    public:
        uint64_t segmentPairIndex = invalid<uint64_t>;
        Strand strand = invalid<Strand>;
        uint64_t frequency = 0;
        ReadOccurrence() {}
        ReadOccurrence(uint64_t segmentPairIndex, Strand strand) :
            segmentPairIndex(segmentPairIndex), strand(strand) {}
        bool operator==(const ReadOccurrence& that) const
        {
            return tie(segmentPairIndex, strand) == tie(that.segmentPairIndex, that.strand);
        }
        bool operator<(const ReadOccurrence& that) const
        {
            return tie(segmentPairIndex, strand) < tie(that.segmentPairIndex, that.strand);
        }
    };
    std::map<ReadId, vector<ReadOccurrence> > readOccurrenceMap;



    // For strand separation we construct an undirected graph
    // with a vertex for each of the segments in the segmentsVector.
    // The vertex_descriptor of this graph is the index of the Segment
    // in the segments vector.
    // An edge s0-s1 is created if there are reads that appear in the
    // same strand on s0 and s1.
    class Vertex {
    public:
        Segment segment;
        uint64_t component = invalid<uint64_t>;
        Vertex(Segment segment = assemblyGraphNullEdge) : segment(segment) {}
    };
    class Edge {
    public:
        uint64_t frequency;
        bool isCrossStrandEdge = false;
        Edge(uint64_t frequency) : frequency(frequency) {}
    };
    using GraphBaseClass = boost::adjacency_list<
        boost::setS,
        boost::vecS,
        boost::undirectedS,
        Vertex,
        Edge>;
    class Graph: public GraphBaseClass {
    public:
        void addToEdge(
            uint64_t segmentIndex0,
            uint64_t segmentIndex1,
            uint64_t frequency);

        // Pairs of reverse complemented edges in the Graph,
        // sorted by decreasing frequency.
        class EdgePair {
        public:
            Graph::edge_descriptor e;
            Graph::edge_descriptor eRc;
            uint64_t frequency;
            bool operator<(const EdgePair& that) const
            {
                return frequency > that.frequency;
            }
        };
        vector<EdgePair> edgePairs;
        void findEdgePairs();
    };
    Graph graph;
    void createGraph();

    void separateStrands();
};
