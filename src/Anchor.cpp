// Shasta.
#include "Anchor.hpp"
#include "color.hpp"
#include "deduplicate.hpp"
#include "ExternalAnchors.hpp"
#include "html.hpp"
#include "invalid.hpp"
#include "Journeys.hpp"
#include "MarkerInfo.hpp"
#include "MarkerKmers.hpp"
#include "orderPairs.hpp"
#include "orderVectors.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "runCommandWithTimeout.hpp"
#include "timestamp.hpp"
#include "tmpDirectory.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/dynamic_bitset.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include <cmath>
#include <queue>
#include <sstream>

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<Anchors>;



// Constructor to access existing Anchors.
Anchors::Anchors(
    const string& baseName,
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    bool writeAccess) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    baseName(baseName),
    reads(reads),
    k(k),
    kHalf(k/2)
{
    anchorMarkerInfos.accessExisting(largeDataName(baseName + "-AnchorMarkerInfos"), writeAccess);
    anchorData.accessExisting(largeDataName(baseName + "-AnchorData"), writeAccess);
}



Anchor Anchors::operator[](AnchorId anchorId) const
{
    return anchorMarkerInfos[anchorId];
}



Kmer Anchors::anchorKmer(AnchorId anchorId) const
{
    // Get the first AnchorMarkerInterval for this Anchor.
    const Anchor anchor = (*this)[anchorId];
    const AnchorMarkerInfo& firstMarkerInfo = anchor.front();

    return firstMarkerInfo.getKmer(k, reads);
}



uint64_t Anchors::size() const
{
    return anchorMarkerInfos.size();
}



void Anchors::check() const
{
    const Anchors& anchors = *this;

    for(AnchorId anchorId=0; anchorId<size(); anchorId++) {
        const Anchor& anchor = anchors[anchorId];
        anchor.check();
    }
}



void Anchor::check() const
{
    const Anchor& anchor = *this;

    // Check that the ReadIds are in strictly increasing order.
    for(uint64_t i=1; i<size(); i++) {
        SHASTA2_ASSERT(anchor[i-1].orientedReadId.getReadId() < anchor[i].orientedReadId.getReadId());
    }
}



// Return the number of common oriented reads between two Anchors,
// counting only oriented reads that have a greater ordinal on anchorId1
// than they have on anchorId0.
uint64_t Anchors::countCommon(
    AnchorId anchorId0,
    AnchorId anchorId1) const
{
    const Anchors& anchors = *this;
    const Anchor anchor0 = anchors[anchorId0];
    const Anchor anchor1 = anchors[anchorId1];

    auto it0 = anchor0.begin();
    auto it1 = anchor1.begin();

    const auto end0 = anchor0.end();
    const auto end1 = anchor1.end();

    uint64_t count = 0;
    while((it0 != end0) and (it1 != end1)) {
        const OrientedReadId orientedReadId0 = it0->orientedReadId;
        const OrientedReadId orientedReadId1 = it1->orientedReadId;
        if(orientedReadId0 < orientedReadId1) {
            ++it0;
        } else if(orientedReadId1 < orientedReadId0) {
            ++it1;
        } else {
            if(it0->position < it1->position) {
                ++count;
            }
            ++it0;
            ++it1;
        }
    }

    return count;
}



// Same as above, but also compute the average offset in bases.
uint64_t Anchors::countCommon(
    AnchorId anchorId0,
    AnchorId anchorId1,
    uint64_t& baseOffset) const
{
    const Anchors& anchors = *this;
    const Anchor anchor0 = anchors[anchorId0];
    const Anchor anchor1 = anchors[anchorId1];

    auto it0 = anchor0.begin();
    auto it1 = anchor1.begin();

    const auto end0 = anchor0.end();
    const auto end1 = anchor1.end();

    uint64_t count = 0;
    uint64_t sumBaseOffsets = 0;
    while((it0 != end0) and (it1 != end1)) {
        const OrientedReadId orientedReadId0 = it0->orientedReadId;
        const OrientedReadId orientedReadId1 = it1->orientedReadId;
        if(orientedReadId0 < orientedReadId1) {
            ++it0;
        } else if(orientedReadId1 < orientedReadId0) {
            ++it1;
        } else {

            // We found a common oriented read.
            const uint32_t position0 = it0->position;
            const uint32_t position1 = it1->position;
            if(position0 < position1) {
                ++count;
                sumBaseOffsets += position1 - position0;
            }

            ++it0;
            ++it1;
        }
    }

    baseOffset = uint64_t(std::round(double(sumBaseOffsets) / double(count)));
    return count;
}



void Anchors::analyzeAnchorPair(
    AnchorId anchorIdA,
    AnchorId anchorIdB,
    AnchorPairInfo& info
    ) const
{
    const Anchors& anchors = *this;

    // Prepare for the joint loop over OrientedReadIds of the two Anchors.
    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];
    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    // Store the total number of OrientedReadIds on the two edges.
    info.totalA = endA - beginA;
    info.totalB = endB - beginB;


    // Joint loop over the MarkerIntervals of the two Anchors,
    // to count the common oriented reads with positive offset
    // and compute average offsets.
    info.commonForwardAdjacent = 0;
    info.commonForwardNonAdjacent = 0;
    info.commonBackward = 0;
    info.minOffsetInBases = std::numeric_limits<uint64_t>::max();
    info.maxOffsetInBases = 0;
    int64_t sumBaseOffsets = 0;
    auto itA = beginA;
    auto itB = beginB;
    while(itA != endA and itB != endB) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common OrientedReadId.

        // Compute the offset in bases.
        const uint64_t positionA = itA->position;
        const uint64_t positionB = itB->position;

        // Update.
        if(positionA < positionB) {
            const uint64_t journeyOffset = itB->positionInJourney - itA->positionInJourney;
            if(journeyOffset == 1) {
                ++info.commonForwardAdjacent;
            } else {
                ++info.commonForwardNonAdjacent;
            }
            const uint64_t offsetInBases = positionB - positionA;
            sumBaseOffsets += offsetInBases;
            info.minOffsetInBases = min(info.minOffsetInBases, offsetInBases);
            info.maxOffsetInBases = max(info.maxOffsetInBases, offsetInBases);
        } else {
            ++info.commonBackward;
        }

        // Continue the joint loop.
        ++itA;
        ++itB;

    }
    info.onlyA = info.totalA - info.commonForward() - info.commonBackward;
    info.onlyB = info.totalB - info.commonForward() - info.commonBackward;

    // If there are no common reads with positive offset, this is all we can do.
    if(info.commonForward() == 0) {
        info.offsetInBases = invalid<uint64_t>;
        info.onlyAShort = invalid<uint64_t>;
        info.onlyBShort = invalid<uint64_t>;
        info.minOffsetInBases = invalid<uint64_t>;
        info.maxOffsetInBases = invalid<uint64_t>;
        return;
    }

    // Compute the estimated offsets.
    info.offsetInBases = uint64_t(std::round(double(sumBaseOffsets) / double(info.commonForward())));



    // Now do the joint loop again, and count the onlyA and onlyB oriented reads
    // that are too short to appear in the other edge.
    itA = beginA;
    itB = beginB;
    uint64_t onlyACheck = 0;
    uint64_t onlyBCheck = 0;
    info.onlyAShort = 0;
    info.onlyBShort = 0;
    while(true) {
        if(itA == endA and itB == endB) {
            break;
        }

        else if(itB == endB or ((itA!=endA) and (itA->orientedReadId < itB->orientedReadId))) {
            // This oriented read only appears in Anchor A.
            ++onlyACheck;
            const OrientedReadId orientedReadId = itA->orientedReadId;
            const int64_t lengthInBases = int64_t(reads.getReadSequenceLength(orientedReadId.getReadId()));

            const int64_t positionA = itA->position;

            // Find the hypothetical positions of anchor B, assuming the estimated base offset.
            const int64_t positionB = positionA + info.offsetInBases;

            // If this ends up outside the read, this counts as onlyAShort.
            if(positionB < 0 or positionB >= lengthInBases) {
                ++info.onlyAShort;
            }

            ++itA;
            continue;
        }

        else if(itA == endA or ((itB!=endB) and (itB->orientedReadId < itA->orientedReadId))) {
            // This oriented read only appears in Anchor B.
            ++onlyBCheck;
            const OrientedReadId orientedReadId = itB->orientedReadId;
            const int64_t lengthInBases = int64_t(reads.getReadSequenceLength(orientedReadId.getReadId()));

            // Get the positions of edge B in this oriented read.
            const int64_t positionB = itB->position;

            // Find the hypothetical positions of anchor A, assuming the estimated base offset.
            const int64_t positionA = positionB - info.offsetInBases;

            // If this ends up outside the read, this counts as onlyBShort.
            if(positionA < 0 or positionA >= lengthInBases) {
                ++info.onlyBShort;
            }

            ++itB;
            continue;
        }

        else {
            // This oriented read appears in both anchors. In this loop, we
            // don't need to do anything.
            ++itA;
            ++itB;
        }
    }
    SHASTA2_ASSERT(onlyACheck == info.onlyA);
    SHASTA2_ASSERT(onlyBCheck == info.onlyB);
}



void Anchors::writeHtml(
    AnchorId anchorIdA,
    AnchorId anchorIdB,
    AnchorPairInfo& info,
    const Journeys& journeys,
    ostream& html) const
{
    const Anchors& anchors = *this;

    // Begin the summary table.
    html <<
        "<table>"
        "<tr><th><th>On<br>anchor A<th>On<br>anchor B";

    // Total.
    html <<
        "<tr><th class=left>Total<td class=centered>" << info.totalA << "<td class=centered>" << info.totalB;

    // Common.
    html << "<tr><th class=left>Common, forward<td class=centered colspan=2>" <<
        info.commonForward();
    html << "<tr><th class=left>Common, backward<td class=centered colspan=2>" <<
        info.commonBackward;

    // Only.
    html <<
        "<tr><th class=left>Only ";
    writeInformationIcon(html, "The number of oriented reads that appear in one anchor but not the other.");
    html <<
        "<td class=centered>" << info.onlyA << "<td class=centered>" << info.onlyB;

    // The rest of the summary table can only be written if there are common reads with positive offset.
    if(info.commonForward() > 0) {

        // Only, short.
        html <<
            "<tr><th class=left>Only, short<td class=centered>" <<
            info.onlyAShort << "<td class=centered>" << info.onlyBShort;

            // Only, missing.
            html <<
                "<tr><th class=left>Only, missing<td class=centered>" <<
                info.onlyA - info.onlyAShort << "<td class=centered>" << info.onlyB - info.onlyBShort;
    }

    // End the summary table.
    html << "</table>";



    // Only write out the rest if there are common reads with positive offset.
    if(info.commonForward() == 0) {
        return;
    }

    // Write the table with Jaccard similarity and estimated offsets.
    using std::fixed;
    using std::setprecision;
    html <<
        "<br><table>"
        "<tr><th class=left>Corrected Jaccard similarity<td class=centered>" <<
        fixed << setprecision(2) << info.correctedJaccard() <<
        "<tr><th class=left>Estimated offset in bases<td class=centered>" << info.offsetInBases <<
        "<tr><th class=left>Minimum offset in bases<td class=centered>" << info.minOffsetInBases <<
        "<tr><th class=left>Maximum offset in bases<td class=centered>" << info.maxOffsetInBases <<
        "</table>";



    // Write the details table.
    html <<
        "<br>In the following table, positions in red are hypothetical, based on the above "
        "estimated base offset."
        "<p><table>";

    // Header row.
    html <<
        "<tr>"
        "<th class=centered rowspan=2>Oriented<br>read id"
        "<th class=centered colspan=2>Length"
        "<th colspan=2>Anchor A"
        "<th colspan=2>Anchor B"
        "<th colspan=2>Offset"
        "<th rowspan=2>Classification"
        "<tr>"
        "<th>Bases"
        "<th>Anchors"
        "<th>Base<br>Position"
        "<th>Position<br>in journey"
        "<th>Base<br>Position"
        "<th>Position<br>in journey"
        "<th>Base<br>Position"
        "<th>Position<br>in journey";

    // Prepare for the joint loop over OrientedReadIds of the two anchors.
    const auto markerIntervalsA = anchors[anchorIdA];
    const auto markerIntervalsB = anchors[anchorIdB];
    const auto beginA = markerIntervalsA.begin();
    const auto beginB = markerIntervalsB.begin();
    const auto endA = markerIntervalsA.end();
    const auto endB = markerIntervalsB.end();

    // Joint loop over the AnchorMarkerIntervals of the two Anchors.
    auto itA = beginA;
    auto itB = beginB;
    while(true) {
        if(itA == endA and itB == endB) {
            break;
        }

        else if(itB == endB or ((itA!=endA) and (itA->orientedReadId < itB->orientedReadId))) {
            // This oriented read only appears in Anchor A.
            const OrientedReadId orientedReadId = itA->orientedReadId;
            const int64_t lengthInBases = int64_t(reads.getReadSequenceLength(orientedReadId.getReadId()));
            const auto journey = journeys[orientedReadId];

            // Get the positions of Anchor A in this oriented read.
            const int64_t positionA = itA->position;

            // Find the hypothetical positions of Anchor B, assuming the estimated base offset.
            const int64_t positionB = positionA + info.offsetInBases;
            const bool isShort = positionB<0 or positionB >= lengthInBases;

            html <<
                "<tr><td class=centered>"
                "<a href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&strand=" << orientedReadId.getStrand() << "'>" << orientedReadId << "</a>"
                "<td class=centered>" << lengthInBases <<
                "<td class=centered>" << journey.size() <<
                "<td class=centered>" << positionA <<
                "<td class=centered>" << itA->positionInJourney <<
                "<td class=centered style='color:Red'>" << positionB <<
                "<td class=centered style='color:Red'>" << "<td><td>"
                "<td class=centered>OnlyA, " << (isShort ? "short" : "missing");

            ++itA;
            continue;
        }

        else if(itA == endA or ((itB!=endB) and (itB->orientedReadId < itA->orientedReadId))) {
            // This oriented read only appears in Anchor B.
            const OrientedReadId orientedReadId = itB->orientedReadId;
            const int64_t lengthInBases = int64_t(reads.getReadSequenceLength(orientedReadId.getReadId()));
            const auto journey = journeys[orientedReadId];

            // Get the positions of Anchor B in this oriented read.
            const int64_t positionB = itB->position;

            // Find the hypothetical positions of edge A, assuming the estimated base offset.
            const int64_t positionA = positionB - info.offsetInBases;
            const bool isShort = positionA<0 or positionA >= lengthInBases;

            html <<
                "<tr><td class=centered>"
                "<a href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&strand=" << orientedReadId.getStrand() << "'>" << orientedReadId << "</a>"
                "<td class=centered>" << lengthInBases <<
                "<td class=centered>" << journey.size() <<
                "<td class=centered style='color:Red'>" << positionA <<
                "<td>"
                "<td class=centered>" << positionB <<
                "<td class=centered>" << itB->positionInJourney <<
                "<td class=centered>" << "<td>"
                "<td class=centered>OnlyB, " << (isShort ? "short" : "missing");

            ++itB;
            continue;
        }

        else {
            // This oriented read appears in both Anchors.
            const OrientedReadId orientedReadId = itA->orientedReadId;
            const int64_t lengthInBases = int64_t(reads.getReadSequenceLength(orientedReadId.getReadId()));
            const auto journey = journeys[orientedReadId];

            // Get the positions of Anchor A in this oriented read.
            const int64_t positionA = itA->position;

            // Get the positions of Anchor B in this oriented read.
            const int64_t positionB = itB->position;

            // Compute estimated offsets.
            const int64_t baseOffset = positionB - positionA;

            html <<
                "<tr><td class=centered>"
                "<a href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&strand=" << orientedReadId.getStrand() << "'>" << orientedReadId << "</a>"
                "<td class=centered>" << lengthInBases <<
                "<td class=centered>" << journey.size() <<
                "<td class=centered>" << positionA <<
                "<td class=centered>" << itA->positionInJourney <<
                "<td class=centered>" << positionB <<
                "<td class=centered>" << itB->positionInJourney <<
                "<td class=centered>" << baseOffset <<
                "<td class=centered>" << int64_t(itB->positionInJourney) - int64_t(itA->positionInJourney) <<
                "<td class=centered>Common";

            ++itA;
            ++itB;
        }
    }

    // Finish the details table.
    html << "</table>";

}



// Anchors are numbered such that each pair of reverse complemented
// AnchorIds are numbered (n, n+1), where n is even, n = 2*m.
// We represent an AnchorId as a string as follows:
// - AnchorId n is represented as m+
// - AnchorId n+1 is represented as m-
// For example, the reverse complemented pair (150, 151) is represented as (75+, 75-).

string shasta2::anchorIdToString(AnchorId n)
{
    std::ostringstream s;

    const AnchorId m = (n >> 1);
    s << m;

    if(n & 1) {
        s << "-";
    } else {
        s << "+";
    }

    return s.str();
}



AnchorId shasta2::anchorIdFromString(const string& s)
{

    if(s.size() < 2) {
        return invalid<AnchorId>;
    }

    const char cLast = s.back();

    uint64_t lastBit;
    if(cLast == '+') {
        lastBit = 0;
    } else if(cLast == '-') {
        lastBit = 1;
    } else {
        return invalid<AnchorId>;
    }

    const uint64_t m = std::stoul(s.substr(0, s.size()-1));

    return 2* m + lastBit;
}



// For a given AnchorId, follow the read journeys forward by one step.
// Return a vector of the AnchorIds reached in this way.
// The count vector is the number of oriented reads each of the AnchorIds.
void Anchors::findChildren(
    const Journeys& journeys,
    AnchorId anchorId,
    vector<AnchorId>& children,
    vector<uint64_t>& count,
    uint64_t minCoverage) const
{
    children.clear();
    for(const auto& markerInfo: anchorMarkerInfos[anchorId]) {
        const OrientedReadId orientedReadId = markerInfo.orientedReadId;
        const auto journey = journeys[orientedReadId];
        const uint64_t position = markerInfo.positionInJourney;
        const uint64_t nextPosition = position + 1;
        if(nextPosition < journey.size()) {
            const AnchorId nextAnchorId = journey[nextPosition];
            children.push_back(nextAnchorId);
        }
    }

    deduplicateAndCountWithThreshold(children, count, minCoverage);
}



// For a given AnchorId, follow the read journeys backward by one step.
// Return a vector of the AnchorIds reached in this way.
// The count vector is the number of oriented reads each of the AnchorIds.
void Anchors::findParents(
    const Journeys& journeys,
    AnchorId anchorId,
    vector<AnchorId>& parents,
    vector<uint64_t>& count,
    uint64_t minCoverage) const
{
    parents.clear();
    for(const auto& markerInfo: anchorMarkerInfos[anchorId]) {
        const OrientedReadId orientedReadId = markerInfo.orientedReadId;
        const auto journey = journeys[orientedReadId];
        const uint64_t position = markerInfo.positionInJourney;
        if(position > 0) {
            const uint64_t previousPosition = position - 1;
            const AnchorId previousAnchorId = journey[previousPosition];
            parents.push_back(previousAnchorId);
        }
    }

    deduplicateAndCountWithThreshold(parents, count, minCoverage);
}



// Get the position for the AnchorMarkerInfo corresponding to a
// given AnchorId and OrientedReadId.
// This asserts if the given AnchorId does not contain an AnchorMarkerInfo
// for the requested OrientedReadId.
uint32_t Anchors::getPosition(AnchorId anchorId, OrientedReadId orientedReadId) const
{
    for(const auto& markerInfo: anchorMarkerInfos[anchorId]) {
        if(markerInfo.orientedReadId == orientedReadId) {
            return markerInfo.position;
        }
    }

    SHASTA2_ASSERT(0);
}



// Get the positionInJourney for the AnchorMarkerInfo corresponding to a
// given AnchorId and OrientedReadId.
// This asserts if the given AnchorId does not contain an AnchorMarkerInfo
// for the requested OrientedReadId.
uint32_t Anchors::getPositionInJourney(AnchorId anchorId, OrientedReadId orientedReadId) const
{
    for(const auto& markerInfo: anchorMarkerInfos[anchorId]) {
        if(markerInfo.orientedReadId == orientedReadId) {
            return markerInfo.positionInJourney;
        }
    }

    SHASTA2_ASSERT(0);
}


// Get the AnchorMarkerInfo corresponding to a given AnchorId and OrientedReadId.
// This asserts if the given AnchorId does not contain an AnchorMarkerInfo
// for the requested OrientedReadId.
const AnchorMarkerInfo& Anchors::getAnchorMarkerInfo(AnchorId anchorId, OrientedReadId orientedReadId) const
{
    for(const auto& markerInfo: anchorMarkerInfos[anchorId]) {
        if(markerInfo.orientedReadId == orientedReadId) {
            return markerInfo;
        }
    }

    SHASTA2_ASSERT(0);
}


// Find out if the given AnchorId contains the specified OrientedReadId.
bool Anchors::anchorContains(AnchorId anchorId, OrientedReadId orientedReadId) const
{
    for(const auto& markerInfo: anchorMarkerInfos[anchorId]) {
        if(markerInfo.orientedReadId == orientedReadId) {
            return true;
        }
    }

    return false;
}



void Anchors::writeCoverageHistogram() const
{
    vector<uint64_t> histogram;
    for(AnchorId anchorId=0; anchorId<anchorMarkerInfos.size(); anchorId++) {
        const uint64_t coverage = anchorMarkerInfos.size(anchorId);
        if(coverage >= histogram.size()) {
            histogram.resize(coverage + 1, 0);
        }
        ++histogram[coverage];
    }

    ofstream csv("AnchorCoverageHistogram.csv");
    csv << "Coverage,Frequency\n";
    for(uint64_t coverage=0; coverage<histogram.size(); coverage++) {
        csv << coverage << "," << histogram[coverage] << "\n";
    }
}



// Constructor to create Anchors from MarkerKmers.
Anchors::Anchors(
    const string& baseName,
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MarkerKmers& markerKmers,
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage,
    const vector<uint64_t>& maxAnchorRepeatLength,
    const vector<uint64_t>& minAnchorDistinctSubkmerCount,
    uint64_t threadCount) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    baseName(baseName),
    reads(reads),
    k(k),
    kHalf(k/2)
{

    performanceLog << timestamp << "Anchor creation begins." << endl;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store arguments so all threads can see them.
    ConstructData& data = constructData;
    data.markerKmersPointer = &markerKmers;
    data.minAnchorCoverage = minAnchorCoverage;
    data.maxAnchorCoverage = maxAnchorCoverage;
    data.maxAnchorRepeatLength = maxAnchorRepeatLength;
    data.minAnchorDistinctSubkmerCount = minAnchorDistinctSubkmerCount;

    // During multithreaded pass 1 we loop over all marker k-mers
    // and for each one we find out if it can be used to generate
    // a pair of anchors or not. If it can be used,
    // we also fill in the coverage - that is,
    // the number of usable MarkerInfos that will go in each of the
    // two anchors.
    const uint64_t markerKmerCount = markerKmers.size();
    data.coverage.createNew(largeDataName("tmp-kmerToAnchorInfos"), largeDataPageSize);
    data.coverage.resize(markerKmerCount);
    const uint64_t batchSize = 1000;
    setupLoadBalancing(markerKmerCount, batchSize);
    runThreads(&Anchors::constructThreadFunctionPass1, threadCount);



    // Assign AnchorIds to marker k-mers and allocate space for
    // each anchor.
    anchorMarkerInfos.createNew(
            largeDataName(baseName + "-AnchorMarkerInfos"),
            largeDataPageSize);
    anchorInfos.createNew(largeDataName(baseName + "-AnchorInfos"), largeDataPageSize);
    AnchorId anchorId = 0;
    for(uint64_t kmerIndex=0; kmerIndex<markerKmerCount; kmerIndex++) {
        const uint64_t coverage = data.coverage[kmerIndex];
        if(coverage == 0) {
            // This k-mer does not generate any anchors.
        } else {
            // This k-mer generates two anchors.
            anchorMarkerInfos.appendVector(coverage);
            anchorMarkerInfos.appendVector(coverage);
            anchorInfos.push_back(AnchorInfo(kmerIndex));
            anchorInfos.push_back(AnchorInfo(kmerIndex));
            anchorId += 2;
        }
    }
    data.coverage.remove();
    const uint64_t anchorCount = anchorMarkerInfos.size();
    SHASTA2_ASSERT(anchorId == anchorCount);
    SHASTA2_ASSERT(anchorInfos.size() == anchorCount);

    anchorMarkerInfos.unreserve();
    anchorInfos.unreserve();



    // In pass 2 we fill in the AnchorMarkerInfos for each anchor.
    SHASTA2_ASSERT((batchSize % 2) == 0);
    setupLoadBalancing(anchorCount, batchSize);
    runThreads(&Anchors::constructThreadFunctionPass2, threadCount);

    anchorInfos.remove();

    // Initialize the AnchorData.
    anchorData.createNew(largeDataName(baseName + "-AnchorData"), largeDataPageSize);
    anchorData.resize(anchorCount);
    std::ranges::fill(anchorData, AnchorData());    // Probably superfluous.

    cout << "Number of anchors per strand: " << anchorCount / 2 << endl;
    performanceLog << timestamp << "Anchor creation from marker kmers ends." << endl;

    // check();

}




void Anchors::constructThreadFunctionPass1(uint64_t /* threadId */)
{

    ConstructData& data = constructData;
    const MarkerKmers& markerKmers = *(data.markerKmersPointer);
    const uint64_t minAnchorCoverage = data.minAnchorCoverage;
    const uint64_t maxAnchorCoverage = data.maxAnchorCoverage;
    const vector<uint64_t> maxAnchorRepeatLength = data.maxAnchorRepeatLength;
    const vector<uint64_t> minAnchorDistinctSubkmerCount = data.minAnchorDistinctSubkmerCount;

    // Loop over batches of marker Kmers assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over marker k-mers assigned to this batch.
        for(uint64_t kmerIndex=begin; kmerIndex!=end; kmerIndex++) {

            SHASTA2_ASSERT(data.coverage[kmerIndex] == 0);

            // Get the MarkerInfos for this marker Kmer.
            const span<const MarkerInfo> markerInfos = markerKmers[kmerIndex];

            // If coverage is too low, don't generate an anchor.
            if(markerInfos.size() < minAnchorCoverage) {
                continue;
            }

            // Check for high coverage using all of the marker infos.
            if(markerInfos.size() > maxAnchorCoverage) {
                continue;
            }

            // Count the usable MarkerInfos.
            // These are the ones for which the ReadId is different from the ReadId
            // of the previous and next MarkerInfo.
            uint64_t  usableMarkerInfosCount = 0;
            for(uint64_t i=0; i<markerInfos.size(); i++) {
                const MarkerInfo& markerInfo = markerInfos[i];
                bool isUsable = true;

                // Check if same ReadId of previous MarkerInfo.
                if(i != 0) {
                    isUsable =
                        isUsable and
                        (markerInfo.orientedReadId.getReadId() != markerInfos[i-1].orientedReadId.getReadId());
                }

                // Check if same ReadId of next MarkerInfo.
                if(i != markerInfos.size() - 1) {
                    isUsable =
                        isUsable and
                        (markerInfo.orientedReadId.getReadId() != markerInfos[i+1].orientedReadId.getReadId());
                }

                if(isUsable) {
                    ++usableMarkerInfosCount;
                }
            }

            if(markerInfos.size() - usableMarkerInfosCount > 0) {
                continue;
            }

            // If coverage is too low, don't generate an anchor.
            if(usableMarkerInfosCount < minAnchorCoverage) {
                continue;
            }

            // Check for repeats.
            bool skipDueToRepeats = false;
            const Kmer kmer = markerInfos.front().getKmer(k, reads);
            for(uint64_t i=0; i<maxAnchorRepeatLength.size(); i++) {
                const uint64_t period = i + 1;
                const uint64_t maxAllowedCopyNumber = maxAnchorRepeatLength[i];
                if(kmer.countExactRepeatCopies(period, k) > maxAllowedCopyNumber) {
                    skipDueToRepeats = true;
                    break;
                }
            }
            if(skipDueToRepeats) {
                continue;
            }

#if 1
            // Check for low complexity by counting distinct sub-k-mers.
            bool skipDueToLowComplexitySequence = false;
            for(uint64_t i=0; i<minAnchorDistinctSubkmerCount.size(); i++) {
                const uint64_t subKmerLength = i + 1;
                const uint64_t minAllowedCount = minAnchorDistinctSubkmerCount[i];
                if(kmer.count(subKmerLength, k) < minAllowedCount) {
                    skipDueToLowComplexitySequence = true;
                    // kmer.write(cout, k);
                    // cout << " skipped due to low complexity sequence at length " << subKmerLength << endl;
                    break;
                }
            }
            if(skipDueToLowComplexitySequence) {
                continue;
            }
#endif

            // If getting here, we will generate a pair of Anchors corresponding to this Kmer.
            data.coverage[kmerIndex] = usableMarkerInfosCount;
        }
    }
}



void Anchors::constructThreadFunctionPass2(uint64_t /* threadId */)
{

    ConstructData& data = constructData;
    const MarkerKmers& markerKmers = *(data.markerKmersPointer);
    const uint64_t minAnchorCoverage = data.minAnchorCoverage;
    const uint64_t maxAnchorCoverage = data.maxAnchorCoverage;


    // A vector used below and defined here to reduce memory allocation activity.
    // It will contain the MarkerInfos for a marker Kmer, excluding
    // the ones for which the same ReadId appears more than once in the same Kmer.
    // There are the ones that will be used to generate anchors.
    vector<MarkerInfo> usableMarkerInfos;

    // Loop over batches of AnchorIds assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over marker k-mers assigned to this batch.
        for(AnchorId anchorId=begin; anchorId!=end; anchorId+=2) {
            const uint64_t kmerIndex = anchorInfos[anchorId].kmerIndex;

            // Get the MarkerInfos for this marker Kmer.
            const span<const MarkerInfo> markerInfos = markerKmers[kmerIndex];

            // We already checked for high coverage during pass 1.
            SHASTA2_ASSERT(markerInfos.size() <= maxAnchorCoverage);

            // Gather the usable MarkerInfos.
            // These are the ones for which the ReadId is different from the ReadId
            // of the previous and next MarkerInfo.
            usableMarkerInfos.clear();
            for(uint64_t i=0; i<markerInfos.size(); i++) {
                const MarkerInfo& markerInfo = markerInfos[i];
                bool isUsable = true;

                // Check if same ReadId of previous MarkerInfo.
                if(i != 0) {
                    isUsable =
                        isUsable and
                        (markerInfo.orientedReadId.getReadId() != markerInfos[i-1].orientedReadId.getReadId());
                }

                // Check if same ReadId of next MarkerInfo.
                if(i != markerInfos.size() - 1) {
                    isUsable =
                        isUsable and
                        (markerInfo.orientedReadId.getReadId() != markerInfos[i+1].orientedReadId.getReadId());
                }

                if(isUsable) {
                    usableMarkerInfos.push_back(markerInfo);
                }
            }

            if(markerInfos.size() - usableMarkerInfos.size() > 0) {
                continue;
            }

            // We already checked for low coverage durign pass1.
            SHASTA2_ASSERT(usableMarkerInfos.size() >= minAnchorCoverage);

            // Fill in the AnchorMarkerInfos for this anchor.
            const auto& anchorMarkerInfos0 = anchorMarkerInfos[anchorId];
            SHASTA2_ASSERT(anchorMarkerInfos0.size() == usableMarkerInfos.size());
            copy(usableMarkerInfos.begin(), usableMarkerInfos.end(), anchorMarkerInfos0.begin());

            // Reverse complement the usableMarkerInfos, then
            // generate the second anchor in the pair.
            for(MarkerInfo& markerInfo: usableMarkerInfos) {
                markerInfo = markerInfo.reverseComplement(reads);
            }
            const auto& anchorMarkerInfos1 = anchorMarkerInfos[anchorId + 1];
            SHASTA2_ASSERT(anchorMarkerInfos1.size() == usableMarkerInfos.size());
            copy(usableMarkerInfos.begin(), usableMarkerInfos.end(), anchorMarkerInfos1.begin());
        }
    }
}



// Constructor to read Anchors from ExternalAnchors.
Anchors::Anchors(
    const string& baseName,
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const string& externalAnchorsName) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    baseName(baseName),
    reads(reads),
    k(k),
    kHalf(k/2)
{

    // Access the ExternalAnchors.
    SHASTA2_ASSERT(not externalAnchorsName.empty());
    if(externalAnchorsName[0] != '/') {
        throw runtime_error("--external-anchors-name must specify an absolute path.");
    }
    cout << "Reading external anchors " << externalAnchorsName << endl;
    const ExternalAnchors externalAnchors(externalAnchorsName, ExternalAnchors::AccessExisting());
    cout << "Found " << externalAnchors.data.size() <<
        " external anchors with average coverage " <<
        externalAnchors.data.totalSize() / externalAnchors.data.size() << endl;

    // Initialize the binary data owned by Anchors.
    anchorMarkerInfos.createNew(
        largeDataName(baseName + "-AnchorMarkerInfos"),
        largeDataPageSize);



    // Loop over external anchors.
    // Each external anchor generates a pair of Anchors.
    vector<AnchorMarkerInfo> markerInfos;
    for(uint64_t i=0; i<externalAnchors.data.size(); i++) {
        const span<const ExternalAnchors::OrientedRead> externalAnchor = externalAnchors.data[i];

        // Create the MarkerInfos for this external anchor.
        // Also check that the Kmers for all oriented reads are identical.
        Kmer kmer;
        markerInfos.clear();
        for(const ExternalAnchors::OrientedRead& orientedRead: externalAnchor) {
            const OrientedReadId orientedReadId = orientedRead.orientedReadId;
            const uint32_t position = orientedRead.position;

            // Check the Kmer.
            const Kmer orientedReadKmer = reads.getKmer(k, orientedReadId, position);
            if(markerInfos.empty()) {
                kmer = orientedReadKmer;
            } else {
                if(orientedReadKmer != kmer) {
                    std::ostringstream message;
                    message << "Inconsistent kmer at oriented read " << orientedReadId <<
                        " position " << position << endl;
                    cout << message.str() << endl;
                    cout << "Offending external anchor:" << endl;
                    externalAnchors.write(cout, i, k, reads);
                    throw runtime_error(message.str());
                }
            }

            // Store this MarkerInfo.
            markerInfos.emplace_back(orientedReadId, position + kHalf);
        }

        // Sort the MarkerInfos by OrientedReadId.
        sort(markerInfos.begin(), markerInfos.end());

        // Check that the ReadIds are all distinct.
        for(uint64_t i1=1; i1<markerInfos.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            if(markerInfos[i0].orientedReadId.getReadId() >= markerInfos[i1].orientedReadId.getReadId()) {
                std::ostringstream message;
                message << "Duplicate ReadId " << markerInfos[i0].orientedReadId.getReadId() << endl;
                cout << "Offending external anchor:" << endl;
                externalAnchors.write(cout, i, k, reads);
                throw runtime_error(message.str());
            }
        }

        // Generate the Anchor corresponding to this external anchor,
        // without reverse complementing.
        // We are not filling AnchorInfo::kmerIndex.
        anchorMarkerInfos.appendVector(markerInfos);

        // Reverse complement the MarkerInfos, then generate
        // the reverse complemented Anchor.
        // We are not filling AnchorInfo::kmerIndex.
        for(AnchorMarkerInfo& markerInfo: markerInfos) {
            markerInfo = markerInfo.reverseComplement(reads);
        }
        anchorMarkerInfos.appendVector(markerInfos);
    }

    // Initialize the AnchorData.
    anchorData.createNew(largeDataName(baseName + "-AnchorData"), largeDataPageSize);
    anchorData.resize(anchorMarkerInfos.size());
    std::ranges::fill(anchorData, AnchorData());    // Probably superfluous.

    cout << "Generated " << anchorMarkerInfos.size() << " anchors from " <<
        externalAnchors.data.size() << " external anchors." << endl;
}



void Anchors::remove()
{
    anchorMarkerInfos.remove();
    anchorInfos.remove();
}



#if 0
// This is the old version that requires the Journeys to be updated
// at each iteration by removing the bad anchors from the Journeys.
uint64_t Anchors::flagBadAnchors(
    const Journeys& journeys,
    uint64_t coverageThreshold)
{

    Anchors& anchors = *this;
    vector<AnchorId> nextOrPrevious;
    vector<uint64_t> count;

    ofstream csv("BadAnchors.csv");

    // Loop over positive (even) anchors.
    uint64_t flaggedCount = 0;
    for(AnchorId anchorId=0; anchorId<anchors.size(); anchorId+=2) {
        if(anchorData[anchorId].isBad) {
            continue;
        }
        const Anchor anchor = anchors[anchorId];

        // Find the next anchors in Journeys.
        nextOrPrevious.clear();
        for(const auto& markerInfo: anchor) {
            const OrientedReadId orientedReadId = markerInfo.orientedReadId;
            const auto journey = journeys[orientedReadId];
            const uint64_t position = markerInfo.positionInJourney;
            const uint64_t nextPosition = position + 1;
            if(nextPosition < journey.size()) {
                const AnchorId nextAnchorId = journey[nextPosition];
                nextOrPrevious.push_back(nextAnchorId);
            }
        }
        // Count how many times each of them appears.
        deduplicateAndCount(nextOrPrevious, count);
        const uint64_t maxForwardCoverage = (count.empty() ? 0 : std::ranges::max(count));

        // Do the same, backward.
        nextOrPrevious.clear();
        for(const auto& markerInfo: anchor) {
            const OrientedReadId orientedReadId = markerInfo.orientedReadId;
            const auto journey = journeys[orientedReadId];
            const uint64_t position = markerInfo.positionInJourney;
            if(position > 0) {
                const uint64_t previousPosition = position - 1;
                const AnchorId previousAnchorId = journey[previousPosition];
                nextOrPrevious.push_back(previousAnchorId);
            }
        }
        // Count how many times each of them appears.
        deduplicateAndCount(nextOrPrevious, count);
        const uint64_t maxBackwardCoverage = (count.empty() ? 0 : std::ranges::max(count));

        // If this is a terminal Anchor, don't flag it it as bad.
        if(maxForwardCoverage == 0) {
            continue;
        }
        if(maxBackwardCoverage == 0) {
            continue;
        }

        // The minimum of maxForwardCoverage and maxBackwardCoverage
        // must must be at least equal to coverageThreshold.
        // If that is not the case, the anchor is flagged as bad.
        const uint64_t n = min(maxForwardCoverage, maxBackwardCoverage);
        if((n < coverageThreshold)) {
            anchorData[anchorId].isBad = true;
            anchorData[anchorId + 1].isBad = true;
            flaggedCount += 2;
            csv << anchorIdToString(anchorId) << "," << anchorIdToString(anchorId + 1) << "\n";
        }
    }

    cout << "Flagged " << flaggedCount << " anchors as bad out of " <<
        anchors.size() << " total." << endl;

    return flaggedCount;
}
#endif



// This is the new version that works on the initial Jurneys
// created using all of the Anchors.
uint64_t Anchors::flagBadAnchors(
    const Journeys& journeys,
    uint64_t coverageThreshold)
{

    Anchors& anchors = *this;
    vector<AnchorId> nextOrPrevious;
    vector<uint64_t> count;

    // Loop over positive (even) anchors.
    uint64_t flaggedCount = 0;
    for(AnchorId anchorId=0; anchorId<anchors.size(); anchorId+=2) {
        if(anchorData[anchorId].isBad) {
            continue;
        }
        const Anchor anchor = anchors[anchorId];

        // Find the next anchors in Journeys, excluding bad anchors.
        nextOrPrevious.clear();
        for(const auto& markerInfo: anchor) {
            const OrientedReadId orientedReadId = markerInfo.orientedReadId;
            const auto journey = journeys[orientedReadId];
            const uint64_t position = markerInfo.positionInJourney;
            for(uint64_t nextPosition=position+1; nextPosition<journey.size(); nextPosition++) {
                const AnchorId nextAnchorId = journey[nextPosition];
                if(not anchorData[nextAnchorId].isBad) {
                    nextOrPrevious.push_back(nextAnchorId);
                    break;
                }
            }
        }
        // Count how many times each of them appears.
        deduplicateAndCount(nextOrPrevious, count);
        const uint64_t maxForwardCoverage = (count.empty() ? 0 : std::ranges::max(count));

        // Do the same, backward.
        nextOrPrevious.clear();
        for(const auto& markerInfo: anchor) {
            const OrientedReadId orientedReadId = markerInfo.orientedReadId;
            const auto journey = journeys[orientedReadId];
            const uint64_t position = markerInfo.positionInJourney;
            if(position > 0) {
                for(uint64_t previousPosition=position-1; /* Check later */ ; previousPosition--) {
                    const AnchorId previousAnchorId = journey[previousPosition];
                    if(not anchorData[previousAnchorId].isBad) {
                        nextOrPrevious.push_back(previousAnchorId);
                        break;
                    }
                    if(previousPosition == 0) {
                        break;
                    }
                }
            }
        }
        // Count how many times each of them appears.
        deduplicateAndCount(nextOrPrevious, count);
        const uint64_t maxBackwardCoverage = (count.empty() ? 0 : std::ranges::max(count));

        // If this is a terminal Anchor, don't flag it it as bad.
        if(maxForwardCoverage == 0) {
            continue;
        }
        if(maxBackwardCoverage == 0) {
            continue;
        }

        // The minimum of maxForwardCoverage and maxBackwardCoverage
        // must must be at least equal to coverageThreshold.
        // If that is not the case, the anchor is flagged as bad.
        const uint64_t n = min(maxForwardCoverage, maxBackwardCoverage);
        if((n < coverageThreshold)) {
            anchorData[anchorId].isBad = true;
            anchorData[anchorId + 1].isBad = true;
            flaggedCount += 2;
        }
    }

    return flaggedCount;
}

