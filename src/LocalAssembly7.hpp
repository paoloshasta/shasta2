#pragma once

// Shasta.
#include "Base.hpp"
#include "invalid.hpp"
#include "orderPairs.hpp"
#include "ReadId.hpp"
#include "shastaTypes.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "iosfwd.hpp"
#include "memory.hpp"
#include "vector.hpp"

namespace shasta2 {
    class LocalAssembly7;

    class Anchors;
    class HomopolymerModel;
}



class shasta2::LocalAssembly7 {
public:

    enum class Method {
        Adaptive,
        Abpoa,
        TheseusOnly,    // Only use oriented reads that are on both anchors
        TheseusAll,     // Also use oriented reads that are om just one anchor.
        DeBruijn,
        Invalid
    };


    class Options {
    public:

        // For reads fixed on one side only, we use a sequence length
        // equal to aExtend * offset + bExtend.
        double aExtend = 0.1;
        double bExtend = 30.;

        // The Method chosen for the local assembly.
        Method method = Method::Adaptive;
        void setMethod(const string&);

        // Whether to repair the bad homopolymer regions of the alignment
        // computed by whichever Method ran, using msaRepair (see msaRepair.hpp).
        // Independent of Method, so any Method that produces a row alignment
        // (Adaptive, Abpoa, TheseusOnly, TheseusAll) can be repaired
        // or not. Mirrors the global Options::useMsaRepair - see OptionsDefine.hpp.
        bool useMsaRepair = false;

        // If the number of oriented reads on both anchors is at least
        // equal to commonThreshold, the adaptive method uses one of:
        // - Fast path.
        // - Abpoa.
        // Otherwise it uses Theseus.
        uint64_t commonCoverageThreshold = 6;

        // The Fast path is used if all of the following is true:
        // - The number of oriented reads on both anchors is at least
        //   equal to commonThreshold.
        // - allowFastPath is true.
        // - Of the oriented reads on both anchors, at least a fraction
        //   equal to fastPathFractionThreshold have the same sequence.
        bool allowFastPath = true;
        double fastPathFractionThreshold = 0.6;

        // If the number of oriented reads on both anchors is at least
        // equal to commonThreshold and the fast path cannot be used,
        // we use abpoa if the maximum length of a sequence
        // fixed on both sides is less than maxAbpoaLength, and
        // theseus otherwise.
        uint64_t maxAbpoaLength = 5000;

        Options() {}
        Options(
            double aExtend,
            double bExtend,
            const string& methodString,
            uint64_t commonCoverageThreshold,
            bool allowFastPath,
            double fastPathFractionThreshold,
            uint64_t maxAbpoaLength);


        // THE REMAINING OPTIONS ARE ONLY USED WHEN USING THE DE BRUIJN GRAPH.
        // THIS IS EXPERIMENTAL AND SO THESE VALUES ARE HARDWIRED.
        // EXPOSE WHEN CODE STABILIZES.

        // The k-mer length for De Bruijn graphs.
        const uint64_t k = 32;

        // Coefficient to compute edge weights for the DeBruijn graph.
        // logP = logPCoefficient * coverage, with logPCoefficient in dB.
        // weight = pow(10, -0.1 * logP).
        // Optimal paths are computed using this weight.
        const double logPCoefficient = 10.;

        // These are used to decide if a k-mer should be split
        // into multiple vertices.
        const double aDrift = 0.02;
        const double bDrift = 50.;

    };



    // This assembles sequence between two anchors using an input vector
    // of oriented reads, which must be sorted.
    // * Of the oriented reads given on input, only the ones that appear
    //   in at least one of the two anchors can be used for assembly.
    // * At least one of the input oriented reads must appear on both anchors.
    // When this completes successfully, the assembled sequence is stored
    // in the sequence vector below. If an error occurs,
    // this throws a std::runtime_error.
    // If homopolymerModelPointer is not null, msaRepair uses it to choose the length of
    // long homopolymer runs (see HomopolymerModel.hpp). Only used if
    // Options::useMsa1 is set.
    LocalAssembly7(
        const Options&,
        const Anchors&,
        shared_ptr<const HomopolymerModel> homopolymerModelPointer,
        AnchorId anchorIdA,
        AnchorId anchorIdB,
        ostream& html,
        const vector<OrientedReadId>& orientedReadIds);

    // Assembled sequence and its coverage.
    bool success = false;
    vector<Base> sequence;

private:

    // Parameters filled in by the constructor.
    Options options;
    const Anchors& anchors;
    shared_ptr<const HomopolymerModel> homopolymerModelPointer;
    AnchorId anchorIdA;     // Left Anchor.
    AnchorId anchorIdB;     // Right Anchor.
    ostream& html;          // Pass ostream(0) to suppress html output.



    // The oriented reads used in this local assembly.
    class OrientedRead {
    public:
        OrientedReadId orientedReadId;

        // The base positions of this OrientedReadId in the two Anchors, if any.
        // These are positions in the oriented read sequence
        // of the mid point of the marker.
        uint32_t positionA = invalid<uint32_t>;
        uint32_t positionB = invalid<uint32_t>;
        uint32_t positionOffsetAB() const
        {
            SHASTA2_ASSERT(isOnBothAnchors());
            SHASTA2_ASSERT(positionB > positionA);
            return positionB - positionA;
        }

        // Whether this OrientedReadId appears in the two Anchors.
        bool isOnAnchorA() const {return isValid(positionA);}
        bool isOnAnchorB() const {return isValid(positionB);}
        bool isOnBothAnchors() const
        {
            return isOnAnchorA() and isOnAnchorB();
        }

        // The above portion is filled in by gatherOrientedReads().
        // The rest is filled in by gatherSequences().

        // The region of this oriented read that will be used in this local assembly.
        // The positions are positions of the marker midpoints.
        uint32_t positionBegin = invalid<uint32_t>;
        uint32_t positionEnd = invalid<uint32_t>;

        uint32_t positionOffsetForAssembly() const
        {
            return positionEnd - positionBegin;
        }

        uint64_t sequenceId = invalid<uint64_t>;

    };
    vector<OrientedRead> orientedReads;
    void gatherOrientedReads(const vector<OrientedReadId>&);
    void removeOutliers();
    static bool checkOffsets(uint64_t, uint64_t);

    // Use the reads fixed on both sides to estimate the offset
    // between anchorIdA and anchorIdB.
    uint32_t offset;
    void estimateOffset();



    // The distinct sequences to be used for assembly.
    using Kmer = vector<Base>;
    class SequenceInfo {
    public:
        bool isOnAnchorA = false;
        bool isOnAnchorB = false;
        vector<OrientedReadId> orientedReadIds;
        uint64_t coverage() const
        {
            return orientedReadIds.size();
        }
        vector<Base> sequence;
        SequenceInfo(
            bool isOnAnchorA,
            bool isOnAnchorB,
            OrientedReadId orientedReadId,
            const vector<Base>& sequence) :
            isOnAnchorA(isOnAnchorA),
            isOnAnchorB(isOnAnchorB),
            orientedReadIds(1, orientedReadId),
            sequence(sequence)
        {
        }

        // The sequence used for the DeBruijn graph is the same as the sequence, but:
        // - If isOnAnchorA, k copies of Base::fromInteger(10) are added at the beginning.
        // - If isOnAnchorB, k copies of Base::fromInteger(20) are added at the end.
        vector<Base> deBruijnSequence;
        void constructDeBruijnSequence(uint64_t k);

        // The k-mers of the deBruijnSequence.
        vector<Kmer> kmers;
        void constructKmers(uint64_t k);

    };
    vector<SequenceInfo> sequences;
    void gatherSequences();

    void writeOrientedReads() const;
    void writeSequences() const;



    // Get  sequenceIds for the SequenceInfos on both anchors,
    // sorted by decreasing coverage.
    // These are used for assembly with abpoa or theseus.
    // SequenceId is the index in the sequences vector above.
    void getSequencesOnBothAnchors(vector<uint64_t>&) const;

    // Same, but for SequennceInfos on one anchor only.
    // These are used for local assembly with theseus.
    void getSequencesOnAnchorA(vector<uint64_t>&) const;
    void getSequencesOnAnchorB(vector<uint64_t>&) const;



    // This returns the total coverage in sequences that are on both anchors.
    uint64_t getTotalCommonCoverage() const;

    // This returns the maximul length of sequences that are on both anchors.
    uint64_t getMaxLengthCommon() const;

    void run();
    void runFastPath();
    void runAdaptive();
    void runAbpoa();
    void runTheseus(bool useAll);


    // Functions and data to find the consensus using a De Bruijn graph
    void runDeBruijn();

    // An occurrence of a k-mer in one of our sequences.
    class KmerOccurrence {
    public:
        // The id of the sequence where this k-mer occurs.
        uint64_t sequenceId;
        // The position of the first base of the k-mer in the deBruijnSequence.
        uint64_t position;
    };

    // The k-mers of all the sequences and their occurrences.
    vector< pair<Kmer, vector<KmerOccurrence> > > kmers;
    void gatherKmers();

    // A vertex of the De Bruijn graph.
    class Vertex {
    public:
        uint64_t vertexId;  // Only for display purposes.
        uint64_t kmerId;
        vector<KmerOccurrence> occurrences;
        uint64_t coverage = invalid<uint64_t>;
        bool isAVertex = false;
        bool isBVertex = false;
        bool isOnAssemblyPath = false;
        Vertex() {}
        Vertex(
            uint64_t vertexId,
            uint64_t kmerId,
            const vector<KmerOccurrence>& occurrences,
            uint64_t coverage) :
            vertexId(vertexId),
            kmerId(kmerId),
            occurrences(occurrences),
            coverage(coverage)
        {}
    };



    // An edge of the De Bruijn graph.
    class Edge {
    public:
        uint64_t coverage = 0;
        double weight = 0.;
    };



    // The De Bruijn graph.
    using GraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        Vertex,
        Edge>;
    class Graph : public GraphBaseClass {
    public:
        using Base = shasta2::Base;
        vertex_descriptor vA;
        vertex_descriptor vB;
        uint64_t nextVertexId = 0;
        void removeUnreachableVertices();

        void writeGraphviz(const string& fileName, const vector< pair<Kmer, vector<KmerOccurrence> > >& kmers) const;
        void writeGraphviz(ostream&, const vector< pair<Kmer, vector<KmerOccurrence> > >& kmers) const;
        void writeVertices(
            const vector< pair<Kmer, vector<KmerOccurrence> > >& kmers,
            const string& fileName) const;
        void writeVertices(
            const vector< pair<Kmer, vector<KmerOccurrence> > >& kmers,
            ostream&) const;

        // The assembly path is a minimum weight path between vA and vB.
        vector<vertex_descriptor> assemblyPath;
        void computeAssemblyPath();
    };
    using vertex_descriptor = Graph::vertex_descriptor;
    using edge_descriptor = Graph::edge_descriptor;
    void createGraph(uint64_t k, Graph&);
    void writeGraph(uint64_t k, const Graph&);

    void writeKmerOccurrences(const Graph&, const string& fileName) const;
    void writeKmerOccurrences(const Graph&, ostream&) const;

    // Use a De Bruijn graph to assemble sequence.
    void assemble(uint64_t k, Graph&);

    static char getCoverageCharacter(uint64_t coverage);
    void writeConsensus(const vector< pair<Base, uint64_t> >&) const;
    void writeAlignment(
        const vector< vector<AlignedBase> >& alignment,
        const vector<AlignedBase>& alignedConsensus,
        const vector< pair<Base, uint64_t> >& consensus,
        const vector< pair<uint64_t, uint64_t> > sequenceIdsWithWeight
        ) const;
    void writeSequence() const;
    void writeAssemblyPath(const Graph&) const;
};
