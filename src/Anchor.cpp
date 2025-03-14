#include "Anchor.hpp"
#include "deduplicate.hpp"
#include "html.hpp"
#include "invalid.hpp"
#include "Journeys.hpp"
#include "MarkerInfo.hpp"
#include "MarkerKmers.hpp"
#include "Markers.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace shasta;

#include <cmath>

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<Anchors>;



// This constructor accesses existing Anchors.
Anchors::Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const Markers& markers,
    bool writeAccess) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    k(k),
    markers(markers)
{
    SHASTA_ASSERT((k %2) == 0);
    kHalf = k / 2;

    anchorMarkerInfos.accessExisting(largeDataName("AnchorMarkerInfos"), writeAccess);
    anchorInfos.accessExistingReadOnly(largeDataName("AnchorInfos"));
    kmerToAnchorTable.accessExistingReadOnly(largeDataName("KmerToAnchorTable"));
}



Anchor Anchors::operator[](AnchorId anchorId) const
{
    return anchorMarkerInfos[anchorId];
}



// This returns the sequence of the marker k-mer
// that this anchor was created from.
vector<Base> Anchors::anchorKmerSequence(AnchorId anchorId) const
{
    // Get the first AnchorMarkerInterval for this Anchor.
    const Anchor anchor = (*this)[anchorId];
    const AnchorMarkerInfo& firstMarkerInfo = anchor.front();

    // Get the OrientedReadId and the ordinals.
    const OrientedReadId orientedReadId = firstMarkerInfo.orientedReadId;
    const uint32_t ordinal = firstMarkerInfo.ordinal;

    // Access the markers of this OrientedReadId.
    const auto orientedReadMarkers = markers[orientedReadId.getValue()];

    const Marker& marker = orientedReadMarkers[ordinal];

    const uint32_t begin = marker.position;
    const uint32_t end = begin + uint32_t(k);

    vector<Base> sequence;
    for(uint32_t position=begin; position!=end; position++) {
        sequence.push_back(reads.getOrientedReadBase(orientedReadId, position));
    }

    return sequence;
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

    for(uint64_t i=1; i<size(); i++) {
        SHASTA_ASSERT(anchor[i-1].orientedReadId.getReadId() < anchor[i].orientedReadId.getReadId());
    }
}



// Return the number of common oriented reads between two Anchors.
uint64_t Anchors::countCommon(
    AnchorId anchorId0,
    AnchorId anchorId1,
    bool ignoreNegativeOffsets) const
{
    const Anchors& anchors = *this;
    const Anchor anchor0 = anchors[anchorId0];
    const Anchor anchor1 = anchors[anchorId1];

    return anchor0.countCommon(anchor1, ignoreNegativeOffsets);
}



// Return the number of common oriented reads with another Anchor.
// Oriented reads in each Anchor are sorted and not duplicated.
uint64_t Anchor::countCommon(const Anchor& that, bool ignoreNegativeOffsets) const
{
    const Anchor& anchor0 = *this;
    const Anchor& anchor1 = that;

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
            if(ignoreNegativeOffsets){
                if(it0->ordinal < it1->ordinal) {
                    ++count;
                }
            } else {
                ++count;
            }
            ++it0;
            ++it1;
        }
    }

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
    // to count the common oreiented reads and compute average offsets.
    info.common = 0;
    int64_t sumMarkerOffsets = 0;
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
        ++info.common;
        const OrientedReadId orientedReadId = itA->orientedReadId;
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];

        // Compute the offset in markers.
        const uint32_t ordinalA = itA->ordinal;
        const uint32_t ordinalB = itB->ordinal;
        sumMarkerOffsets += int64_t(ordinalB) - int64_t(ordinalA);

        // Compute the offset in bases.
        const int64_t positionA = int64_t(orientedReadMarkers[ordinalA].position);
        const int64_t positionB = int64_t(orientedReadMarkers[ordinalB].position);
        sumBaseOffsets += positionB - positionA;

        // Continue the joint loop.
        ++itA;
        ++itB;

    }
    info.onlyA = info.totalA - info.common;
    info.onlyB = info.totalB - info.common;

    // If there are no common reads, this is all we can do.
    if(info.common == 0) {
        info.offsetInMarkers = invalid<int64_t>;
        info.offsetInBases = invalid<int64_t>;
        info.onlyAShort = invalid<uint64_t>;
        info.onlyBShort = invalid<uint64_t>;
        return;
    }

    // Compute the estimated offsets.
    info.offsetInMarkers = int64_t(std::round(double(sumMarkerOffsets) / double(info.common)));
    info.offsetInBases = int64_t(std::round(double(sumBaseOffsets) / double(info.common)));



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
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(reads.getReadSequenceLength(orientedReadId.getReadId()));

            const uint32_t ordinalA = itA->ordinal;
            const int64_t positionA = int64_t(orientedReadMarkers[ordinalA].position);

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
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(reads.getReadSequenceLength(orientedReadId.getReadId()));

            // Get the positions of edge B in this oriented read.
            const uint32_t ordinalB = itB->ordinal;
            const int64_t positionB = int64_t(orientedReadMarkers[ordinalB].position);

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
            // This oriented read appears in both edges. In this loop, we
            // don't need to do anything.
            ++itA;
            ++itB;
        }
    }
    SHASTA_ASSERT(onlyACheck == info.onlyA);
    SHASTA_ASSERT(onlyBCheck == info.onlyB);
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
        "<tr><th class=left>Total ";
    writeInformationIcon(html, "The total number of oriented reads on each of the two anchors.");
    html << "<td class=centered>" << info.totalA << "<td class=centered>" << info.totalB;

    // Common.
    html << "<tr><th class=left>Common ";
    writeInformationIcon(html, "The number of common oriented reads between the two anchors.");
    html <<
        "<td class=centered colspan = 2>" << info.common;

    // Only.
    html <<
        "<tr><th class=left>Only ";
    writeInformationIcon(html, "The number of oriented reads that appear in one anchor but not the other.");
    html <<
        "<td class=centered>" << info.onlyA << "<td class=centered>" << info.onlyB;

    // The rest of the summary table can only be written if there are common reads.
    if(info.common > 0) {

        // Only, short.
        html <<
            "<tr><th class=left>Only, short ";
        writeInformationIcon(html, "The number of oriented reads that appear in one anchor only "
            " and are too short to appear on the other anchor, based on the estimated base offset.");
        html <<
            "<td class=centered>" << info.onlyAShort << "<td class=centered>" << info.onlyBShort;

        // Only, missing.
        html <<
            "<tr><th class=left>Only, missing ";
        writeInformationIcon(html, "The number of oriented reads that appear in one anchor only "
            " and are not too short to appear on the other anchor, based on the estimated base offset.");
        html <<
            "<td class=centered>" << info.onlyA - info.onlyAShort << "<td class=centered>" << info.onlyB - info.onlyBShort;
    }

    // End the summary table.
    html << "</table>";

    // Only write out the rest if there are common reads.
    if(info.common == 0) {
        return;
    }

    // Write the table with Jaccard similarities and estimated offsets.
    using std::fixed;
    using std::setprecision;
    html <<
        "<br><table>"
        "<tr><th class=left>Jaccard similarity<td class=centered>" <<
        fixed << setprecision(2) << info.jaccard() <<
        "<tr><th class=left>Corrected Jaccard similarity<td class=centered>" <<
        fixed << setprecision(2) << info.correctedJaccard() <<
        "<tr><th class=left>Estimated offset in markers<td class=centered>" << info.offsetInMarkers <<
        "<tr><th class=left>Estimated offset in bases<td class=centered>" << info.offsetInBases <<
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
        "<th class=centered colspan=3>Length"
        "<th colspan=3>Anchor A"
        "<th colspan=3>Anchor B"
        "<th colspan=3>Offset"
        "<th rowspan=2>Classification"
        "<tr>"
        "<th>Bases"
        "<th>Markers"
        "<th>Anchors"
        "<th>Base<br>Position"
        "<th>Marker<br>ordinal"
        "<th>Position<br>in journey"
        "<th>Base<br>Position"
        "<th>Marker<br>ordinal"
        "<th>Position<br>in journey"
        "<th>Base<br>Position"
        "<th>Marker<br>ordinal"
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
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(reads.getReadSequenceLength(orientedReadId.getReadId()));
            const auto journey = journeys[orientedReadId];

            // Get the positions of Anchor A in this oriented read.
            const uint32_t ordinalA = itA->ordinal;
            const int64_t positionA = int64_t(orientedReadMarkers[ordinalA].position);

            // Find the hypothetical positions of Anchor B, assuming the estimated base offset.
            const int64_t positionB = positionA + info.offsetInBases;
            const bool isShort = positionB<0 or positionB >= lengthInBases;

            html <<
                "<tr><td class=centered>"
                "<a href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&strand=" << orientedReadId.getStrand() << "'>" << orientedReadId << "</a>"
                "<td class=centered>" << lengthInBases <<
                "<td class=centered>" << orientedReadMarkers.size() <<
                "<td class=centered>" << journey.size() <<
                "<td class=centered>" << positionA <<
                "<td class=centered>" << ordinalA <<
                "<td class=centered>" << itA->positionInJourney <<
                "<td class=centered style='color:Red'>" << positionB <<
                "<td>"
                "<td class=centered style='color:Red'>" << "<td><td><td>"
                "<td class=centered>OnlyA, " << (isShort ? "short" : "missing");

            ++itA;
            continue;
        }

        else if(itA == endA or ((itB!=endB) and (itB->orientedReadId < itA->orientedReadId))) {
            // This oriented read only appears in Anchor B.
            const OrientedReadId orientedReadId = itB->orientedReadId;
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(reads.getReadSequenceLength(orientedReadId.getReadId()));
            const auto journey = journeys[orientedReadId];

            // Get the positions of Anchor B in this oriented read.
            const uint32_t ordinalB = itB->ordinal;
            const int64_t positionB = int64_t(orientedReadMarkers[ordinalB].position);

            // Find the hypothetical positions of edge A, assuming the estimated base offset.
            const int64_t positionA = positionB - info.offsetInBases;
            const bool isShort = positionA<0 or positionA >= lengthInBases;

            html <<
                "<tr><td class=centered>"
                "<a href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&strand=" << orientedReadId.getStrand() << "'>" << orientedReadId << "</a>"
                "<td class=centered>" << lengthInBases <<
                "<td class=centered>" << orientedReadMarkers.size() <<
                "<td class=centered>" << journey.size() <<
                "<td class=centered style='color:Red'>" << positionA <<
                "<td><td>"
                "<td class=centered>" << positionB <<
                "<td class=centered>" << ordinalB <<
                "<td class=centered>" << itB->positionInJourney <<
                "<td class=centered>" << "<td><td>"
                "<td class=centered>OnlyB, " << (isShort ? "short" : "missing");

            ++itB;
            continue;
        }

        else {
            // This oriented read appears in both Anchors.
            const OrientedReadId orientedReadId = itA->orientedReadId;
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(reads.getReadSequenceLength(orientedReadId.getReadId()));
            const auto journey = journeys[orientedReadId];

            // Get the positions of Anchor A in this oriented read.
            const uint32_t ordinalA = itA->ordinal;
            const int64_t positionA = int64_t(orientedReadMarkers[ordinalA].position);

            // Get the positions of Anchor B in this oriented read.
            const uint32_t ordinalB = itB->ordinal;
            const int64_t positionB = int64_t(orientedReadMarkers[ordinalB].position);

            // Compute estimated offsets.
            const int64_t ordinalOffset = uint64_t(ordinalB) - uint64_t(ordinalA);
            const int64_t baseOffset = positionB - positionA;

            html <<
                "<tr><td class=centered>"
                "<a href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&strand=" << orientedReadId.getStrand() << "'>" << orientedReadId << "</a>"
                "<td class=centered>" << lengthInBases <<
                "<td class=centered>" << orientedReadMarkers.size() <<
                "<td class=centered>" << journey.size() <<
                "<td class=centered>" << positionA <<
                "<td class=centered>" << ordinalA <<
                "<td class=centered>" << itA->positionInJourney <<
                "<td class=centered>" << positionB <<
                "<td class=centered>" << ordinalB <<
                "<td class=centered>" << itB->positionInJourney <<
                "<td class=centered>" << baseOffset <<
                "<td class=centered>" << ordinalOffset <<
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

string shasta::anchorIdToString(AnchorId n)
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



AnchorId shasta::anchorIdFromString(const string& s)
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



// Get the ordinal for the AnchorMarkerInfo corresponding to a
// given AnchorId and OrientedReadId.
// This asserts if the given AnchorId does not contain an AnchorMarkerInfo
// for the requested OrientedReadId.
uint32_t Anchors::getOrdinal(AnchorId anchorId, OrientedReadId orientedReadId) const
{
    for(const auto& markerInfo: anchorMarkerInfos[anchorId]) {
        if(markerInfo.orientedReadId == orientedReadId) {
            return markerInfo.ordinal;
        }
    }

    SHASTA_ASSERT(0);
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



// Read following.
void Anchors::followOrientedReads(
    const Journeys& journeys,
    AnchorId anchorId0,
    uint64_t direction,                         // 0 = forward, 1 = backward
    uint64_t minCommonCount,
    double minJaccard,
    double minCorrectedJaccard,
    vector< pair<AnchorId, AnchorPairInfo> >& anchorInfos
    ) const
{
    const Anchor anchor0 = (*this)[anchorId0];

    // Gather all AnchorIds reached by the forward or backward portions of the
    // journeys of the oriented reads on this anchor.
    vector<AnchorId> anchorIds;
    for(const AnchorMarkerInfo& anchorMarkerInfo: anchor0) {
        const OrientedReadId orientedReadId = anchorMarkerInfo.orientedReadId;
        const uint64_t position0 = anchorMarkerInfo.positionInJourney;
        const auto journey = journeys[orientedReadId];

        // Figure out the forward or backward portion of the journey.
        uint64_t begin;
        uint64_t end;
        if(direction == 0) {
            begin = position0 + 1;
            end = journey.size();
        } else {
            begin = 0;
            end = position0;
        }

        // Copy the AnchorIds on this portion of the journey.
        copy(journey.begin() + begin, journey.begin() + end, back_inserter(anchorIds));
    }

    // Only keep the ones we saw at least minCommonCount times.
    vector<uint64_t> count;
    deduplicateAndCountWithThreshold(anchorIds, count, minCommonCount);



    // Gather the ones that satisfy our criteria.
    anchorInfos.clear();
    for(const AnchorId anchorId1: anchorIds) {
        AnchorPairInfo info;
        if(direction == 0) {
            analyzeAnchorPair(anchorId0, anchorId1, info);
        } else {
            analyzeAnchorPair(anchorId1, anchorId0, info);
        }
        if(info.common < minCommonCount) {
            continue;
        }
        if(info.common == 0) {
            continue;
        }
        if(info.jaccard() < minJaccard) {
            continue;
        }
        if(info.correctedJaccard() < minCorrectedJaccard) {
            continue;
        }
        anchorInfos.push_back(make_pair(anchorId1, info));
    }


    // Sort them by offset.
    class SortHelper {
    public:
        bool operator()(
            const pair<AnchorId, AnchorPairInfo>& x,
            const pair<AnchorId, AnchorPairInfo>& y
            ) const
        {
            return x.second.offsetInBases < y.second.offsetInBases;
        }
    };
    sort(anchorInfos.begin(), anchorInfos.end(), SortHelper());
}



Anchors::Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const Markers& markers,
    shared_ptr<MarkerKmers> markerKmers,
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage,
    uint64_t maxHomopolymerLength,
    uint64_t threadCount) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    k(k),
    markers(markers)
{
    kHalf = k / 2;

    performanceLog << timestamp << "Anchor creation begins." << endl;

    // Store arguments so all threads can see them.
    ConstructData& data = constructData;
    data.minAnchorCoverage = minAnchorCoverage;
    data.maxAnchorCoverage = maxAnchorCoverage;
    data.maxHomopolymerLength = maxHomopolymerLength;
    data.markerKmers = markerKmers;

    // During multithreaded pass 1 we loop over all marker k-mers
    // and for each one we find out if it can be used to generate
    // a pair of anchors or not. If it can be used,
    // we also fill in the coverage - that is,
    // the number of usable MarkerInfos that will go in each of the
    // two anchors.
    const uint64_t markerKmerCount = markerKmers->size();
    data.coverage.createNew(largeDataName("tmp-kmerToAnchorInfos"), largeDataPageSize);
    data.coverage.resize(markerKmerCount);
    const uint64_t batchSize = 1000;
    setupLoadBalancing(markerKmerCount, batchSize);
    runThreads(&Anchors::constructThreadFunctionPass1, threadCount);



    // Assign AnchorIds to marker k-mers and allocate space for
    // each anchor.
    anchorMarkerInfos.createNew(
            largeDataName("AnchorMarkerInfos"),
            largeDataPageSize);
    anchorInfos.createNew(largeDataName("AnchorInfos"), largeDataPageSize);
    kmerToAnchorTable.createNew(largeDataName("KmerToAnchorTable"), largeDataPageSize);
    kmerToAnchorTable.resize(markerKmerCount);
    AnchorId anchorId = 0;
    for(uint64_t kmerIndex=0; kmerIndex<markerKmerCount; kmerIndex++) {
        const uint64_t coverage = data.coverage[kmerIndex];
        if(coverage == 0) {
            // This k-mer does not generate any anchors.
            kmerToAnchorTable[kmerIndex] = invalid<AnchorId>;
        } else {
            // This k-mer generates two anchors.
            anchorMarkerInfos.appendVector(coverage);
            anchorMarkerInfos.appendVector(coverage);
            kmerToAnchorTable[kmerIndex] = anchorId;
            anchorInfos.push_back(AnchorInfo(kmerIndex));
            anchorInfos.push_back(AnchorInfo(kmerIndex));
            anchorId += 2;
        }
    }
    data.coverage.remove();
    const uint64_t anchorCount = anchorMarkerInfos.size();
    SHASTA_ASSERT(anchorId == anchorCount);
    SHASTA_ASSERT(anchorInfos.size() == anchorCount);



    // In pass 2 we fill in the AnchorMarkerInfos for each anchor.
    SHASTA_ASSERT((batchSize % 2) == 0);
    setupLoadBalancing(anchorCount, batchSize);
    runThreads(&Anchors::constructThreadFunctionPass2, threadCount);



    cout << "Number of anchors per strand: " << anchorCount / 2 << endl;
    performanceLog << timestamp << "Anchor creation from marker kmers ends." << endl;

}




void Anchors::constructThreadFunctionPass1(uint64_t /* threadId */)
{

    ConstructData& data = constructData;
    const uint64_t minAnchorCoverage = data.minAnchorCoverage;
    const uint64_t maxAnchorCoverage = data.maxAnchorCoverage;
    const uint64_t maxHomopolymerLength = data.maxHomopolymerLength;
    const MarkerKmers& markerKmers = *data.markerKmers;

    // Loop over batches of marker Kmers assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over marker k-mers assigned to this batch.
        for(uint64_t kmerIndex=begin; kmerIndex!=end; kmerIndex++) {

            SHASTA_ASSERT(data.coverage[kmerIndex] == 0);

            // Get the MarkerInfos for this marker Kmer.
            const span<const MarkerInfo> markerInfos = markerKmers[kmerIndex];

            // If the Kmer as a long homopolymer run, don't generate an anchor.
            const Kmer kmer = markerKmers.getKmer(markerInfos.front());
            if(kmer.maxHomopolymerLength(k) > maxHomopolymerLength) {
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

            // Check for low coverage using the usable marker infos.
            if(usableMarkerInfosCount >= minAnchorCoverage) {
                // If getting here, we will generate a pair of Anchors corresponding to this Kmer.
                data.coverage[kmerIndex] = usableMarkerInfosCount;
            }
        }
    }
}



void Anchors::constructThreadFunctionPass2(uint64_t /* threadId */)
{

    ConstructData& data = constructData;
    const uint64_t minAnchorCoverage = data.minAnchorCoverage;
    const uint64_t maxAnchorCoverage = data.maxAnchorCoverage;
    const MarkerKmers& markerKmers = *data.markerKmers;


    // A vector used below and defined here to reduce memory allocation activity.
    // It will contain the MarkerInfos for a marker Kmer, excluding
    // the ones for which the same ReadId appears more than once in the same Kmer.
    // There are the ones that will be used to generate anchors.
    vector<MarkerInfo> usableMarkerInfos;

    // Loop over batches of AnchorIds assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over marker k-mers assigned to this batch.
        for(AnchorId anchorId=0; anchorId!=end; anchorId+=2) {
            const uint64_t kmerIndex = anchorInfos[anchorId].kmerIndex;

            // Get the MarkerInfos for this marker Kmer.
            const span<const MarkerInfo> markerInfos = markerKmers[kmerIndex];

            // We already checked for high coverage during pass 1.
            SHASTA_ASSERT(markerInfos.size() <= maxAnchorCoverage);

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

            // We already checked for low coverage durign pass1.
            SHASTA_ASSERT(usableMarkerInfos.size() >= minAnchorCoverage);

            // Fill in the AnchorMarkerInfos for this anchor.
            const auto& anchorMarkerInfos0 = anchorMarkerInfos[anchorId];
            SHASTA_ASSERT(anchorMarkerInfos0.size() == usableMarkerInfos.size());
            copy(usableMarkerInfos.begin(), usableMarkerInfos.end(), anchorMarkerInfos0.begin());

            // Reverse complement the usableMarkerInfos, then
            // generate the second anchor in the pair.
            for(MarkerInfo& markerInfo: usableMarkerInfos) {
                markerInfo = markerInfo.reverseComplement(markers);
            }
            const auto& anchorMarkerInfos1 = anchorMarkerInfos[anchorId + 1];
            SHASTA_ASSERT(anchorMarkerInfos1.size() == usableMarkerInfos.size());
            copy(usableMarkerInfos.begin(), usableMarkerInfos.end(), anchorMarkerInfos1.begin());
        }
    }
}

