#pragma once

// Shasta.
#include "Kmer.hpp"
#include "invalid.hpp"
#include "MappedMemoryOwner.hpp"
#include "MarkerInfo.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "ReadId.hpp"

// Standard library.
#include "cstdint.hpp"
#include "memory.hpp"
#include "span.hpp"



namespace shasta {

    class Base;
    class Marker;
    class MarkerKmers;
    class Markers;
    class MarkerInfo;
    class Reads;

    using AnchorId = uint64_t;
    class Anchor;
    class AnchorMarkerInfo;
    class Anchors;
    class AnchorInfo;
    class AnchorPairInfo;

    using AnchorBaseClass = span<const AnchorMarkerInfo>;

    class Journeys;

    string anchorIdToString(AnchorId);
    AnchorId anchorIdFromString(const string&);
}



// An Anchor is a set of AnchorMarkerInfos.
class shasta::AnchorMarkerInfo : public MarkerInfo {
public:
    uint32_t positionInJourney = invalid<uint32_t>;

    // Default constructor.
    AnchorMarkerInfo() {}

    // Constructor from a MarkerInfo.
    AnchorMarkerInfo(
        const MarkerInfo& markerInfo) :
        MarkerInfo(markerInfo)
    {}

    // Constructor from OrientedReadId and ordinal.
    AnchorMarkerInfo(
        OrientedReadId orientedReadId,
        uint32_t ordinal) :
        MarkerInfo(orientedReadId, ordinal)
    {}
};



class shasta::AnchorInfo {
public:
    // The k-mer index in the MarkerKmers for the k-mer
    // that generated this anchor and its reverse complement.
    uint64_t kmerIndex = invalid<uint64_t>;

    AnchorInfo(uint64_t kmerIndex = invalid<uint64_t>) : kmerIndex(kmerIndex) {}
};



// An Anchor is a set of AnchorMarkerInfos.
class shasta::Anchor : public AnchorBaseClass {
public:

    Anchor(const AnchorBaseClass& s) : AnchorBaseClass(s) {}

    void check() const;

    uint64_t coverage() const
    {
        return size();
    }

    // Return the number of common oriented reads with another Anchor.
    uint64_t countCommon(const Anchor& that, bool ignoreNegativeOffsets = false) const;
};



class shasta::Anchors :
    public MultithreadedObject<Anchors>,
    public MappedMemoryOwner {
public:

    Anchors(
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const Markers& markers,
        shared_ptr<MarkerKmers>,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage,
        uint64_t threadCount);

    // This constructor accesses existing Anchors.
    Anchors(
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const Markers& markers,
        bool writeAccess);

    Anchor operator[](AnchorId) const;
    uint64_t size() const;

    // This returns the sequence of the marker k-mer
    // that this anchor was created from.
    vector<Base> anchorKmerSequence(AnchorId) const;

    // Return the number of common oriented reads between two Anchors.
    uint64_t countCommon(AnchorId, AnchorId, bool ignoreNegativeOffsets = false) const;

    // Analyze the oriented read composition of two anchors.
    void analyzeAnchorPair(AnchorId, AnchorId, AnchorPairInfo&) const;
    void writeHtml(AnchorId, AnchorId, AnchorPairInfo&, ostream&) const;

    void writeCoverageHistogram() const;

    MemoryMapped::VectorOfVectors<AnchorMarkerInfo, uint64_t> anchorMarkerInfos;

private:
    // A MemoryMapped::Vector that gives, for each k-mer in the Marker K-mers,
    // the AnchordId of the first of the two anchors generated by that k-mer,
    // or invalid<AnchorId> if that k-mer did not generate any anchors.
    MemoryMapped::Vector<AnchorId> kmerToAnchorTable;

public:

    // Get the ordinal for the AnchorMarkerInfo corresponding to a
    // given AnchorId and OrientedReadId.
    // This asserts if the given AnchorId does not contain an AnchorMarkerInfo
    // for the requested OrientedReadId.
    uint32_t getOrdinal(AnchorId, OrientedReadId) const;

    const Reads& reads;
    uint64_t k;
    uint64_t kHalf;
    const Markers& markers;

private:

    void check() const;

public:

    // For a given AnchorId, follow the read journeys forward/backward by one step.
    // Return a vector of the AnchorIds reached in this way.
    // The count vector is the number of oriented reads each of the AnchorIds.
    void findChildren(
        const Journeys&,
        AnchorId,
        vector<AnchorId>&,
        vector<uint64_t>& count,
        uint64_t minCoverage = 0) const;
    void findParents(
        const Journeys&,
        AnchorId,
        vector<AnchorId>&,
        vector<uint64_t>& count,
        uint64_t minCoverage = 0) const;


    // In addition to the marker intervals, we also store an AnchorInfo for each Anchor.
    MemoryMapped::Vector<AnchorInfo> anchorInfos;
    void storeAnchorInfo(AnchorId)
    {
        // For now there is nothing to store.
        // AnchorInfo& anchorInfo = anchorInfos[anchorId];
    }

    // Read following.
    void followOrientedReads(
        const Journeys& journeys,
        AnchorId,
        uint64_t direction,                         // 0 = forward, 1 = backward
        uint64_t minCommonCount,
        double minJaccard,
        double minCorrectedJaccard,
        vector< pair<AnchorId, AnchorPairInfo> >&
        ) const;

private:

    // Data and functions used when constructing the Anchors.
    class ConstructData {
    public:
        uint64_t minAnchorCoverage;
        uint64_t maxAnchorCoverage;

        shared_ptr<MarkerKmers> markerKmers;

        // During multithreaded pass 1 we loop over all marker k-mers
        // and for each one we find out if it can be used to generate
        // a pair of anchors or not. If it can be used,
        // we also fill in the coverage - that is,
        // the number of usable MarkerInfos that will go in each of the
        // two anchors.
        MemoryMapped::Vector<uint64_t> coverage;
    };
    ConstructData constructData;
    void constructThreadFunctionPass1(uint64_t threadId);
    void constructThreadFunctionPass2(uint64_t threadId);

};



// Information about the read composition similarity of two anchors A and B.
class shasta::AnchorPairInfo {
public:

    // The total number of OrientedReadIds in each of the anchors A and B.
    uint64_t totalA = 0;
    uint64_t totalB = 0;

    // The number of common oriented reads.
    uint64_t common = 0;

    // The number of oriented reads present in A but not in B.
    uint64_t onlyA = 0;

    // The number of oriented reads present in B but not in A.
    uint64_t onlyB = 0;

    // The rest of the statistics are only valid if the number
    // of common oriented reads is not 0.

    // The estimated offset between the two Anchors.
    // The estimate is done using the common oriented reads.
    int64_t offsetInMarkers = invalid<int64_t>;
    int64_t offsetInBases = invalid<int64_t>;

    // The number of onlyA reads which are too short to be on edge B,
    // based on the above estimated offset.
    uint64_t onlyAShort = invalid<uint64_t>;

    // The number of onlyB reads which are too short to be on edge A,
    // based on the above estimated offset.
    uint64_t onlyBShort = invalid<uint64_t>;

    uint64_t intersectionCount() const
    {
        return common;
    }
    uint64_t unionCount() const {
        return totalA + totalB - common;
    }
    uint64_t correctedUnionCount() const
    {
        return unionCount() - onlyAShort - onlyBShort;
    }
    double jaccard() const
    {
        return double(intersectionCount()) / double(unionCount());
    }
    double correctedJaccard() const
    {
        return double(intersectionCount()) / double(correctedUnionCount());
    }

    void reverse()
    {
        swap(totalA, totalB);
        swap(onlyA, onlyB);
        swap(onlyAShort, onlyBShort);
        offsetInMarkers = - offsetInMarkers;
        offsetInBases = - offsetInBases;
    }

};
