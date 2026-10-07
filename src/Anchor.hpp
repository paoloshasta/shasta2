#pragma once

// Shasta.
#include "AnchorPair.hpp"
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



namespace shasta2 {

    class Base;
    class MarkerKmers;
    class MarkerInfo;
    class Reads;

    using AnchorId = uint64_t;
    class Anchor;
    class AnchorMarkerInfo;
    class Anchors;
    class AnchorInfo;
    class AnchorData;
    class AnchorPairInfo;

    using AnchorBaseClass = span<const AnchorMarkerInfo>;

    class Journeys;

    string anchorIdToString(AnchorId);
    AnchorId anchorIdFromString(const string&);

    inline AnchorId reverseComplementAnchorId(AnchorId anchorId)
    {
        return anchorId ^ 1UL;

    }
}



// An Anchor is a set of AnchorMarkerInfos.
class shasta2::AnchorMarkerInfo : public MarkerInfo {
public:
    uint32_t positionInJourney = invalid<uint32_t>;

    // Default constructor.
    AnchorMarkerInfo() {}

    // Constructor from a MarkerInfo.
    AnchorMarkerInfo(
        const MarkerInfo& markerInfo) :
        MarkerInfo(markerInfo)
    {}

    // Constructor from OrientedReadId and position.
    AnchorMarkerInfo(
        OrientedReadId orientedReadId,
        uint32_t position) :
        MarkerInfo(orientedReadId, position)
    {}

    bool operator<(const AnchorMarkerInfo& that) const
    {
        return orientedReadId < that.orientedReadId;
    }
};



class shasta2::AnchorInfo {
public:
    // The k-mer index in the MarkerKmers for the k-mer
    // that generated this anchor and its reverse complement.
    // This is only used in constructThreadFunctionPass2.
    // When using ExternalAnchors, it is not filled in.
    uint64_t kmerIndex = invalid<uint64_t>;

    AnchorInfo(uint64_t kmerIndex = invalid<uint64_t>) : kmerIndex(kmerIndex) {}
};



class shasta2::AnchorData {
public:
    bool isBad = false;
};



// An Anchor is a set of AnchorMarkerInfos.
class shasta2::Anchor : public AnchorBaseClass {
public:

    Anchor(const AnchorBaseClass& s) : AnchorBaseClass(s) {}

    void check() const;

    uint64_t coverage() const
    {
        return size();
    }
};



class shasta2::Anchors :
    public MultithreadedObject<Anchors>,
    public MappedMemoryOwner {
public:

    // Constructor to create Anchors from MarkerKmers.
    Anchors(
        const string& baseName,
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const MarkerKmers&,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage,
        const vector<uint64_t>& maxAnchorRepeatLength,
        const vector<uint64_t>& minAnchorDistinctSubkmerCount,
        uint64_t threadCount);

    // Constructor to read Anchors from ExternalAnchors.
    Anchors(
        const string& baseName,
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const string& externalAnchorsName);

    // Constructor to accesses existing Anchors from binary data.
    Anchors(
        const string& baseName,
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        bool writeAccess = false);

    void remove();

    Anchor operator[](AnchorId) const;
    uint64_t size() const;

    // This returns the sequence of the marker k-mer
    // that this anchor was created from.
    Kmer anchorKmer(AnchorId) const;

    // Return the number of common oriented reads between two Anchors,
    // counting only oriented reads that have a greater ordinal on anchorId1
    // than they have on anchorId0.
    uint64_t countCommon(AnchorId anchorId0, AnchorId anchorId1) const;

    // Same as above, but also compute the average offset in bases.
    uint64_t countCommon(AnchorId anchorId0, AnchorId anchorId1, uint64_t& baseOffset) const;

    // Analyze the oriented read composition of two anchors.
    void analyzeAnchorPair(AnchorId, AnchorId, AnchorPairInfo&) const;
    void writeHtml(AnchorId, AnchorId, AnchorPairInfo&, const Journeys&
        , ostream&) const;

    void writeCoverageHistogram() const;

    MemoryMapped::VectorOfVectors<AnchorMarkerInfo, uint64_t> anchorMarkerInfos;

    // Get the position for the AnchorMarkerInfo corresponding to a
    // given AnchorId and OrientedReadId.
    // This asserts if the given AnchorId does not contain an AnchorMarkerInfo
    // for the requested OrientedReadId.
    uint32_t getPosition(AnchorId, OrientedReadId) const;

    // Get the positioInJourney for the AnchorMarkerInfo corresponding to a
    // given AnchorId and OrientedReadId.
    // This asserts if the given AnchorId does not contain an AnchorMarkerInfo
    // for the requested OrientedReadId.
    uint32_t getPositionInJourney(AnchorId, OrientedReadId) const;

    // Get the AnchorMarkerInfo corresponding to a given AnchorId and OrientedReadId.
    // This asserts if the given AnchorId does not contain an AnchorMarkerInfo
    // for the requested OrientedReadId.
    const AnchorMarkerInfo& getAnchorMarkerInfo(AnchorId, OrientedReadId) const;


    // Find out if the given AnchorId contains the specified OrientedReadId.
    bool anchorContains(AnchorId, OrientedReadId) const;

    void flagBadAnchors(const Journeys&);

    const string baseName;
    const Reads& reads;
    const uint64_t k;
    const uint64_t kHalf;
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


    MemoryMapped::Vector<AnchorInfo> anchorInfos;
    MemoryMapped::Vector<AnchorData> anchorData;
private:


    // Data and functions used when constructing the Anchors.
    class ConstructData {
    public:
        const MarkerKmers* markerKmersPointer = 0;
        uint64_t minAnchorCoverage;
        uint64_t maxAnchorCoverage;
        vector<uint64_t> maxAnchorRepeatLength;
        vector<uint64_t> minAnchorDistinctSubkmerCount;

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
class shasta2::AnchorPairInfo {
public:

    // The total number of OrientedReadIds in each of the anchors A and B.
    uint64_t totalA = 0;
    uint64_t totalB = 0;

    // The number of common oriented reads with positive offset
    // (anchor B occurs after anchor A in the oriented read).
    // Count separately the ones with journey offset equal to 1
    // ("adjacent") and greater than 1 ("non-adjacent").
    uint64_t commonForwardAdjacent = 0;
    uint64_t commonForwardNonAdjacent = 0;
    uint64_t commonForward() const
    {
        return commonForwardAdjacent + commonForwardNonAdjacent;
    }

    double adjacentFraction() const
    {
        return double(commonForwardAdjacent) / double(commonForward());
    }

    // The number of common oriented reads with negative offset
    // (anchor B occurs before anchor A in the oriented read).
    // Zero offset is not possible if the two anchors are distinct.
    uint64_t commonBackward = 0;

    // The number of oriented reads present in A but not in B.
    uint64_t onlyA = 0;

    // The number of oriented reads present in B but not in A.
    uint64_t onlyB = 0;

    // The rest of the statistics are only valid if the number
    // of common oriented reads is not 0.

    // The estimated offset between the two Anchors.
    // The estimate is done using the common oriented reads
    // with positive offset.
    uint64_t offsetInBases = invalid<int64_t>;
    uint64_t minOffsetInBases  = invalid<int64_t>;
    uint64_t maxOffsetInBases  = invalid<int64_t>;

    // The number of onlyA reads which are too short to be on anchor B,
    // based on the above estimated offset.
    uint64_t onlyAShort = invalid<uint64_t>;

    // The number of onlyB reads which are too short to be on anchor A,
    // based on the above estimated offset.
    uint64_t onlyBShort = invalid<uint64_t>;

    uint64_t intersectionCountPositiveOffset() const
    {
        return commonForward();
    }
    uint64_t unionCount() const {
        return totalA + totalB - commonForward() - commonBackward;
    }
    uint64_t correctedUnionCount() const
    {
        return unionCount() - onlyAShort - onlyBShort;
    }
    double correctedJaccard() const
    {
        return double(intersectionCountPositiveOffset()) / double(correctedUnionCount());
    }

    uint64_t missingCount() const
    {
        return onlyA + onlyB - onlyAShort - onlyBShort + commonBackward;
    }

    uint64_t missingA() const
    {
        return onlyA - onlyAShort - commonBackward;
    }
    uint64_t missingB() const
    {
        return onlyB - onlyBShort - commonBackward;
    }

    double missingAFraction() const
    {
        return double(missingA()) / double(totalA);
    }
    double missingBFraction() const
    {
        return double(missingB()) / double(totalB);
    }
    double minMissingFraction() const
    {
        return min(missingAFraction(), missingBFraction());
    }

};
