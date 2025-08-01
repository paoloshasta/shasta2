// Shasta.
#include "AnchorPair.hpp"
#include "Anchor.hpp"
#include "approximateTopologicalSort.hpp"
#include "color.hpp"
#include "deduplicate.hpp"
#include "graphvizToHtml.hpp"
#include "hcsClustering.hpp"
#include "html.hpp"
#include "HttpServer.hpp"
#include "Journeys.hpp"
#include "Markers.hpp"
#include "orderPairs.hpp"
#include "orderVectors.hpp"
#include "runCommandWithTimeout.hpp"
#include "shastaLapack.hpp"
#include "tmpDirectory.hpp"
#include "Reads.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/dynamic_bitset.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/iterator/function_output_iterator.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include <fstream.hpp>
#include <stdexcept.hpp>



AnchorPair::AnchorPair(
    const Anchors& anchors,
    AnchorId anchorIdA,
    AnchorId anchorIdB,
    bool adjacentInJourney) :
    anchorIdA(anchorIdA),
    anchorIdB(anchorIdB)
{
    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];

    // Loop over common oriented reads between these two anchors.
    // If adjacentInJourney is false, the journey offset is required to be positive.
    // If adjacentInJourney is true, the journey offset is required to be exactly 1.

    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

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

        const OrientedReadId orientedReadId = itA->orientedReadId;
        SHASTA_ASSERT(orientedReadId == itB->orientedReadId);

        if(adjacentInJourney) {
            if(itB->positionInJourney == itA->positionInJourney + 1) {
                orientedReadIds.push_back(orientedReadId);
            }
        } else {
            if(itB->positionInJourney > itA->positionInJourney) {
                orientedReadIds.push_back(orientedReadId);
            }
        }

        ++itA;
        ++itB;
    }
}



// Copy from another AnchorPair, but excluding some OrientedReadIds.
AnchorPair::AnchorPair(
    const AnchorPair& that,
    const vector<OrientedReadId>& excludedOrientedReadIds) :
    anchorIdA(that.anchorIdA),
    anchorIdB(that.anchorIdB)
{
    std::set_difference(
        that.orientedReadIds.begin(), that.orientedReadIds.end(),
        excludedOrientedReadIds.begin(), excludedOrientedReadIds.end(),
        back_inserter(orientedReadIds));
}


// Get positions in journey, ordinals, and base positions
// for each of the two reads and for each of the two anchors.
// The positions returned are the midpoint of the markers
// corresponding to anchorIdA and anchorIdB.
void AnchorPair::get(
    const Anchors& anchors,
    vector< pair<Positions, Positions> >& positions) const
{

    const uint32_t kHalf = uint32_t(anchors.markers.k / 2);
    positions.clear();

    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];

    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    auto itA = beginA;
    auto itB = beginB;
    auto it = orientedReadIds.begin();
    const auto itEnd = orientedReadIds.end();
    while(itA != endA and itB != endB and it != itEnd) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common OrientedReadId.
        const OrientedReadId orientedReadId = itA->orientedReadId;
        SHASTA_ASSERT(orientedReadId == itB->orientedReadId);

        // Only process is this is one of our OrientedReadIds;
        if(orientedReadId == *it) {
            ++it;

            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

            const uint32_t positionInJourneyA = itA->positionInJourney;
            const uint32_t positionInJourneyB = itB->positionInJourney;
            SHASTA_ASSERT(positionInJourneyB >= positionInJourneyA);    // Allow degenerate AnchorPair witn anchorIdA==anchorIdB
            const uint32_t ordinalA = itA->ordinal;
            const uint32_t ordinalB = itB->ordinal;
            const uint32_t positionA = orientedReadMarkers[ordinalA].position + kHalf;
            const uint32_t positionB = orientedReadMarkers[ordinalB].position + kHalf;

            positions.push_back(make_pair(
                Positions(positionInJourneyA, ordinalA, positionA),
                Positions(positionInJourneyB, ordinalB, positionB)
                ));
        }

        ++itA;
        ++itB;
    }
}



// Remove from the AnchorPair OrientedReadIds that have negative offsets.
void AnchorPair::removeNegativeOffsets(const Anchors& anchors)
{
    vector< pair<uint32_t, uint32_t> > ordinals;
    getOrdinals(anchors, ordinals);
    SHASTA_ASSERT(ordinals.size() == orientedReadIds.size());

    vector<OrientedReadId> newOrientedReadIds;
    for(uint64_t i=0; i<orientedReadIds.size(); i++) {
        const auto& p = ordinals[i];
        if(p.second >= p.first) {
            newOrientedReadIds.push_back(orientedReadIds[i]);
        }
    }

    orientedReadIds.swap(newOrientedReadIds);
}




// Same as the above, but also returns compute the sequences.
void AnchorPair::get(
    const Anchors& anchors,
    vector< pair<Positions, Positions> >& positions,
    vector< vector<Base> >& sequences) const
{
    const uint32_t kHalf = uint32_t(anchors.markers.k / 2);
    const Reads& reads = anchors.markers.reads;

    get(anchors, positions);

    sequences.clear();
    sequences.resize(orientedReadIds.size());

    for(uint64_t i=0; i<orientedReadIds.size(); i++) {
        const OrientedReadId orientedReadId = orientedReadIds[i];
        const auto& positionsAB = positions[i];
        vector<Base>& sequence = sequences[i];

        const uint32_t positionA = positionsAB.first.basePosition + kHalf;
        const uint32_t positionB = positionsAB.second.basePosition + kHalf;

        for(uint32_t position=positionA; position!=positionB; position++) {
            sequence.push_back(reads.getOrientedReadBase(orientedReadId, position));
        }
    }
}



// Same as the above, but only compute the ordinals.
void AnchorPair::getOrdinals(
    const Anchors& anchors,
    vector< pair<uint32_t, uint32_t> >& ordinals) const
{
    ordinals.clear();

    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];

    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    auto itA = beginA;
    auto itB = beginB;
    auto it = orientedReadIds.begin();
    const auto itEnd = orientedReadIds.end();
    while(itA != endA and itB != endB and it != itEnd) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common OrientedReadId.
        const OrientedReadId orientedReadId = itA->orientedReadId;
        SHASTA_ASSERT(orientedReadId == itB->orientedReadId);

        // Only process is this is one of our OrientedReadIds;
        if(orientedReadId == *it) {
            ++it;

            const uint32_t ordinalA = itA->ordinal;
            const uint32_t ordinalB = itB->ordinal;

            ordinals.push_back(make_pair(ordinalA, ordinalB));
        }

        ++itA;
        ++itB;
    }

    SHASTA_ASSERT(it == orientedReadIds.end());
}



// Just return the journey positions.
void AnchorPair::getPositionsInJourneys(
    const Anchors& anchors,
    vector< pair<uint32_t, uint32_t> >& positionsInJourneys) const
{
    positionsInJourneys.clear();

    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];

    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    auto itA = beginA;
    auto itB = beginB;
    auto it = orientedReadIds.begin();
    const auto itEnd = orientedReadIds.end();
    while(itA != endA and itB != endB and it != itEnd) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common OrientedReadId.
        const OrientedReadId orientedReadId = itA->orientedReadId;
        SHASTA_ASSERT(orientedReadId == itB->orientedReadId);

        // Only process is this is one of our OrientedReadIds;
        if(orientedReadId == *it) {
            ++it;

            const uint32_t positionInJourneyA = itA->positionInJourney;
            const uint32_t positionInJourneyB = itB->positionInJourney;

            positionsInJourneys.push_back(make_pair(positionInJourneyA, positionInJourneyB));
        }

        ++itA;
        ++itB;
    }

    SHASTA_ASSERT(it == orientedReadIds.end());

}



// This finds AnchorPairs as follows:
// - anchorIdA is as specified.
// - Coverage is at least minCoverage.
// - All oriented reads have a journey offset equal to 1.
void AnchorPair::createChildren(
    const Anchors& anchors,
    const Journeys& journeys,
    AnchorId anchorIdA,
    uint64_t minCoverage,
    vector<AnchorPair>& anchorPairs
    )
{
    // Find possible choices for anchorIdB.
    vector<AnchorId> anchorIdsB;
    vector<uint64_t> coverage;
    anchors.findChildren(journeys, anchorIdA, anchorIdsB, coverage, minCoverage);

    anchorPairs.clear();
    for(const AnchorId anchorIdB: anchorIdsB) {
        anchorPairs.emplace_back(anchors, anchorIdA, anchorIdB, true);
    }
}


// "Join" constructor from two AnchorPairs.
// This constructs a new AnchorPair as follows:
// - anchorIdA is the same as anchorPair0.anchorIdA.
// - anchorIdB is the same as anchorPair1.anchorIdB.
// - orientedReadIds are the intersection of
//   anchorPair0.orientedReadIds and anchorPair1.orientedReadIds,
//   with the additional requirement that the new journey offset
//   is positive. That is, each OrientedReadId of the new AnchorPair
//   visits anchorIdB after visiting the new anchorIdA.
// This constructor is used in detangling.
AnchorPair::AnchorPair(
    const Anchors& anchors,
    const AnchorPair& anchorPair0,
    const AnchorPair& anchorPair1) :
        anchorIdA(anchorPair0.anchorIdA),
        anchorIdB(anchorPair1.anchorIdB)
{
    vector< pair<Positions, Positions> > positions0;
    vector< pair<Positions, Positions> > positions1;
    anchorPair0.get(anchors, positions0);
    anchorPair1.get(anchors, positions1);

    const uint64_t n0 = anchorPair0.size();
    const uint64_t n1 = anchorPair1.size();
    SHASTA_ASSERT(positions0.size() == n0);
    SHASTA_ASSERT(positions1.size() == n1);

    // Joint loop over OrientedReadIds of the two AnchorPairs.
    uint64_t i0 = 0;
    uint64_t i1 = 0;
    while((i0 != n0) and (i1 != n1)) {
        const OrientedReadId orientedReadId0 = anchorPair0.orientedReadIds[i0];
        const OrientedReadId orientedReadId1 = anchorPair1.orientedReadIds[i1];

        if(orientedReadId0 < orientedReadId1) {
            ++i0;
            continue;
        }

        if(orientedReadId1 < orientedReadId0) {
            ++i1;
            continue;
        }

        SHASTA_ASSERT(orientedReadId0 == orientedReadId1);
        const OrientedReadId orientedReadId = orientedReadId0;

        // We found an OrientedReadId common between anchorPair0 and anchorPair1.
        // We also have to check the journey offsets.
        if(positions0[i0].first.positionInJourney < positions1[i1].second.positionInJourney) {
            orientedReadIds.push_back(orientedReadId);
        }

        ++i0;
        ++i1;
    }
}



uint32_t AnchorPair::getAverageOffset(const Anchors& anchors) const
{
    const uint32_t kHalf = uint32_t(anchors.markers.k / 2);

    uint64_t sumBaseOffset = 0;

    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];

    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    auto itA = beginA;
    auto itB = beginB;
    auto it = orientedReadIds.begin();
    const auto itEnd = orientedReadIds.end();
    while(itA != endA and itB != endB and it != itEnd) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common OrientedReadId.
        const OrientedReadId orientedReadId = itA->orientedReadId;
        SHASTA_ASSERT(orientedReadId == itB->orientedReadId);

        // Only process is this is one of our OrientedReadIds;
        if(orientedReadId == *it) {
            ++it;

            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

            const uint32_t ordinalA = itA->ordinal;
            const uint32_t ordinalB = itB->ordinal;
            if(ordinalB < ordinalA) {       // Degenerate AnchorPair with AnchorIdA==AnchorIdB is ok.
                throw runtime_error(
                    "Order violation at anchor pair " +
                    anchorIdToString(anchorIdA) + " " +
                    anchorIdToString(anchorIdB) + " " +
                    orientedReadId.getString() + " ordinals " +
                    to_string(ordinalA) + " " +
                    to_string(ordinalB));
            }
            const uint32_t positionA = orientedReadMarkers[ordinalA].position + kHalf;
            const uint32_t positionB = orientedReadMarkers[ordinalB].position + kHalf;
            SHASTA_ASSERT(positionB >= positionA);      // Degenerate AnchorPair with AnchorIdA==AnchorIdB is ok.

            const uint32_t offset = positionB - positionA;
            sumBaseOffset += offset;
        }

        ++itA;
        ++itB;
    }

    SHASTA_ASSERT(it == orientedReadIds.end());

    return uint32_t(std::round(double(sumBaseOffset) / double(size())));
}



void AnchorPair::getOffsets(
    const Anchors& anchors,
    uint32_t& averageBaseOffset,
    uint32_t& minBaseOffset,
    uint32_t& maxBaseOffset) const
{
    const uint32_t kHalf = uint32_t(anchors.markers.k / 2);

    uint64_t sumBaseOffset = 0;
    minBaseOffset = std::numeric_limits<uint32_t>::max();
    maxBaseOffset = 0;

    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];

    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    auto itA = beginA;
    auto itB = beginB;
    auto it = orientedReadIds.begin();
    const auto itEnd = orientedReadIds.end();
    while(itA != endA and itB != endB and it != itEnd) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common OrientedReadId.
        const OrientedReadId orientedReadId = itA->orientedReadId;
        SHASTA_ASSERT(orientedReadId == itB->orientedReadId);

        // Only process is this is one of our OrientedReadIds;
        if(orientedReadId == *it) {
            ++it;

            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

            const uint32_t ordinalA = itA->ordinal;
            const uint32_t ordinalB = itB->ordinal;
            if(ordinalB < ordinalA) {          // Degenerate AnchorPair with AnchorIdA==AnchorIdB is ok.
                throw runtime_error(
                    "Order violation at anchor pair " +
                    anchorIdToString(anchorIdA) + " " +
                    anchorIdToString(anchorIdB) + " " +
                    orientedReadId.getString() + " ordinals " +
                    to_string(ordinalA) + " " +
                    to_string(ordinalB));
            }
            const uint32_t positionA = orientedReadMarkers[ordinalA].position + kHalf;
            const uint32_t positionB = orientedReadMarkers[ordinalB].position + kHalf;
            SHASTA_ASSERT(positionB > positionA);

            const uint32_t offset = positionB - positionA;
            sumBaseOffset += offset;
            minBaseOffset = min(minBaseOffset, offset);
            maxBaseOffset = max(maxBaseOffset, offset);
        }

        ++itA;
        ++itB;
    }

    SHASTA_ASSERT(it == orientedReadIds.end());

    averageBaseOffset = uint32_t(std::round(double(sumBaseOffset) / double(size())));
}



// Split the AnchorPair into one or more AnchorPairs with consistent offsets.
// In the resulting AnchorPairs, if the position offsets are sorted in
// increasing order, any two adjacent offsets D0 and D1
// will satisfy D1 - D0 <= aDrift + bDrift * (D0 + D1) / 2.
void AnchorPair::splitByOffsets(
    const Anchors& anchors,
    double aDrift,
    double bDrift,
    vector<AnchorPair>& newAnchorPairs) const
{
    vector< pair<Positions, Positions> > positions;
    get(anchors, positions);
    const uint64_t n = orientedReadIds.size();
    SHASTA_ASSERT(positions.size() == n);

    // Gather pairs(index, offset) where index is the index
    // in the OrientedReadIds, vector.
    vector< pair<uint64_t, uint64_t> > offsets;
    for(uint64_t i=0; i<n; i++) {
        const uint32_t positionA = positions[i].first.basePosition;
        const uint32_t positionB = positions[i].second.basePosition;
        SHASTA_ASSERT(positionB >= positionA);  // Allow degenerate AnchorPair with anchorIdA = anchorIdB.
        const uint64_t offset = positionB - positionA;
        offsets.push_back(make_pair(i, offset));
    }
    sort(offsets.begin(), offsets.end(), OrderPairsBySecondOnly<uint64_t, uint64_t>());

    // Find places where we have to split.
    vector<uint64_t> splitPoints;
    splitPoints.push_back(0);
    for(uint64_t i1=1; i1<n; i1++) {
        const uint64_t i0 = i1 - 1;
        const double offset0 = double(offsets[i0].second);
        const double offset1 = double(offsets[i1].second);
        if(offset1 - offset0 > aDrift + .5 * bDrift  * (offset1 + offset0)) {
            splitPoints.push_back(i1);
        }
    }
    splitPoints.push_back(n);

    // Each interval between split points generates a new AnchorPair.
    newAnchorPairs.clear();
    newAnchorPairs.resize(splitPoints.size() - 1);
    for(uint64_t i=0; i<splitPoints.size() -1 ; i++) {
        const uint64_t j0 = splitPoints[i];
        const uint64_t j1 = splitPoints[i + 1];

        newAnchorPairs[i].anchorIdA = anchorIdA;
        newAnchorPairs[i].anchorIdB = anchorIdB;
        for(uint64_t j=j0; j!=j1; j++) {
            newAnchorPairs[i].orientedReadIds.push_back(orientedReadIds[offsets[j].first]);
        }
        sort(newAnchorPairs[i].orientedReadIds.begin(), newAnchorPairs[i].orientedReadIds.end());
    }


    // Sort them by decreasing coverage.
    class SortHelper {
    public:
        bool operator() (const AnchorPair& x, const AnchorPair& y) const
        {
            return x.orientedReadIds.size() > y.orientedReadIds.size();
        }
    };
    sort(newAnchorPairs.begin(), newAnchorPairs.end(), SortHelper());

}



// This returns true if a call to split with the same arguments would not split this Anchor.
// The second and third areguments are work vectors added as arguments for performancew,
bool AnchorPair::isConsistent(
    const Anchors& anchors,
    double aDrift,
    double bDrift,
    vector< pair<Positions, Positions> >& positions,
    vector<uint64_t>& offsets) const
{

    get(anchors, positions);
    const uint64_t n = orientedReadIds.size();
    SHASTA_ASSERT(positions.size() == n);

    // Gather offsets.
    offsets.clear();
    offsets.resize(n);
    for(uint64_t i=0; i<n; i++) {
        const uint32_t positionA = positions[i].first.basePosition;
        const uint32_t positionB = positions[i].second.basePosition;
        SHASTA_ASSERT(positionB > positionA);
        const uint64_t offset = positionB - positionA;
        offsets[i] = offset;
    }
    sort(offsets.begin(), offsets.end());

    for(uint64_t i1=1; i1<n; i1++) {
        const uint64_t i0 = i1 - 1;
        const double offset0 = double(offsets[i0]);
        const double offset1 = double(offsets[i1]);
        if(offset1 - offset0 > aDrift + .5 * bDrift  * (offset1 + offset0)) {
            return false;
        }
    }

    return true;
}



// Count OrientedReadIds in common with another AnchorPair.
uint64_t AnchorPair::countCommon(const AnchorPair& y) const
{
    const AnchorPair& x = *this;

    uint64_t n = 0;
    auto counter = [&n](auto){++n;};

    std::set_intersection(
        x.orientedReadIds.begin(),
        x.orientedReadIds.end(),
        y.orientedReadIds.begin(),
        y.orientedReadIds.end(),
        boost::make_function_output_iterator(counter)
    );

    return n;

}



// Split the AnchorPair using clustering of the oriented read journey portions
// within this AnchorPair.
// The new AnchorPairs are sorted by decreasing size.
void AnchorPair::splitByClustering(
    const Anchors& anchors,
    const Journeys& journeys,
    double clusteringMinJaccard,
    vector<AnchorPair>& newAnchorPairs
    ) const
{
    std::ostream html(0);
    vector< vector<uint64_t> > clusters;
    anchors.clusterAnchorPairOrientedReads(*this, journeys, clusteringMinJaccard, clusters, html);

    // Create an AnchorPair for each cluster.
    // The clusters are sorted by decreasing size, and so the new AnchorPairs will
    // also be sorted by decreasing size.
    newAnchorPairs.clear();
    newAnchorPairs.reserve(clusters.size());
    for(const vector<uint64_t>& cluster: clusters) {
        newAnchorPairs.emplace_back();
        AnchorPair& newAnchorPair = newAnchorPairs.back();
        newAnchorPair.anchorIdA = anchorIdA;
        newAnchorPair.anchorIdB = anchorIdB;
        for(const uint64_t i: cluster) {
            newAnchorPair.orientedReadIds.push_back(orientedReadIds[i]);
        }
    }
}



bool AnchorPair::contains(OrientedReadId orientedReadId) const
{
    const auto it = std::lower_bound(orientedReadIds.begin(), orientedReadIds.end(), orientedReadId);
    return it != orientedReadIds.end() and (*it == orientedReadId);
}



void AnchorPair::computeClusteringMatrix(
    const Journeys& journeys,
    const vector< pair<uint32_t, uint32_t> >& positionsInJourneys,  // As computed by getPositionsInJourneys.
    const vector<AnchorId>& internalAnchorIds,                      // As computed by getInternalAnchorIds.
    Matrix& clusteringMatrix
    ) const
{
    // Initialize the ClusterMatrix to all zeros.
    clusteringMatrix.resize(size(), internalAnchorIds.size());
    for(uint64_t j=0; j<internalAnchorIds.size(); j++) {
        for(uint64_t i=0; i<size(); i++) {
            clusteringMatrix(i, j) = 0.;
        }
    }

    // Loop over all OrientedReadIds to fill it in.
    for(uint64_t i=0; i<size(); i++) {
        const OrientedReadId orientedReadId = orientedReadIds[i];
        const Journey journey = journeys[orientedReadId];

        // Get the positions in the journey of anchorIdA and anchorIdB.
        const auto& p = positionsInJourneys[i];
        const uint32_t positionInJourneyA = p.first;
        const uint32_t positionInJourneyB = p.second;

        // Loop over this portion of the journey, excluding anchorIdA and anchorIdB.
        for(uint64_t position=positionInJourneyA+1; position<positionInJourneyB; position++) {
            const AnchorId anchorId = journey[position];

            // Find the index of this AnchorId in internalAnchorIds.
            const auto it = std::ranges::lower_bound(internalAnchorIds, anchorId);
            SHASTA_ASSERT(it != internalAnchorIds.end());
            SHASTA_ASSERT(*it == anchorId);
            const uint64_t j = it - internalAnchorIds.begin();

            // Set this element of the ClusteringMatrix.
            clusteringMatrix(i, j) = 1.;
        }
    }
}



// Singular value decomposition of the clustering matrix.
void AnchorPair::clusteringMatrixSvd(
    Matrix& clusteringMatrix,
    vector<double>& singularValues,
    Matrix& leftSingularVectors,
    Matrix& rightSingularVectors) const
{
    // Shift all the columns so they have zero average.
    for(uint64_t j=0; j<clusteringMatrix.size2(); j++) {
        double sum = 0.;
        for(uint64_t i=0; i<clusteringMatrix.size1(); i++) {
            sum += clusteringMatrix(i, j);
        }
        const double average = sum / double(clusteringMatrix.size1());
        for(uint64_t i=0; i<clusteringMatrix.size1(); i++) {
            clusteringMatrix(i, j) -= average;
        }
    }

    // Compute the SVD.
    dgesvd(clusteringMatrix, singularValues, leftSingularVectors, rightSingularVectors);

}



// Use the scaled left singular values to compute a distance matrix
// between oriented reads.
void AnchorPair::computeDistanceMatrix(
    uint64_t singularValueCount,    // Only use the first singular values
    const vector<double>& singularValues,
    const Matrix& leftSingularVectors,
    Matrix& distanceMatrix
    ) const
{
    singularValueCount = min(singularValueCount, singularValues.size());

    distanceMatrix.clear();
    distanceMatrix.resize(size(), size());

    vector<double> xA(singularValueCount);
    vector<double> xB(singularValueCount);

    for(uint64_t iA=0; iA<size(); iA++) {
        distanceMatrix(iA, iA) = 0.;

        for(uint64_t k=0; k<singularValueCount; k++) {
            xA[k] = singularValues[k] * leftSingularVectors(iA, k);
        }

        for(uint64_t iB=iA+1; iB<size(); iB++) {

            for(uint64_t k=0; k<singularValueCount; k++) {
                xB[k] = singularValues[k] * leftSingularVectors(iB, k);
            }

            double distance = 0.;
            for(uint64_t k=0; k<singularValueCount; k++) {
                const double d = xB[k] - xA[k];
                distance += d * d;
            }
            distance = sqrt(distance);
            distanceMatrix(iA, iB) = distance;
            distanceMatrix(iB, iA) = distance;
        }
    }
}



// Given the distance matrix, compute a similarity graph
// between OrientedReadIds in which each vertex represents an OrientedReadId and
// an edge is generated for OrientedReadId pairs
// with distance below the given threshold.
// Each vertex stores the id of the cluster it is assigned to.
AnchorPair::OrientedReadIdSimilarityGraph::OrientedReadIdSimilarityGraph(
    const Matrix& distanceMatrix,
    double maxDistance)
{
    using Graph = OrientedReadIdSimilarityGraph;
    Graph& graph = *this;

    const uint64_t n = distanceMatrix.size1();
    SHASTA_ASSERT(distanceMatrix.size2() == n);

    for(uint64_t i=0; i<n; i++) {
        add_vertex(graph);
    }

    for(uint64_t i=0; i<n; i++) {
        for(uint64_t j=i+1; j<n; j++) {
            if(distanceMatrix(i, j) <= maxDistance) {
                add_edge(i, j, graph);
            }
        }
    }

    // Clustering.
    hcsClustering(graph, clusters);

    // Sort the clusters by decreasing size.
    sort(clusters.begin(), clusters.end(), OrderVectorsByDecreasingSize<uint64_t>());

    // Store the clusterId of each vertex.
    for(uint64_t clusterId=0; clusterId<clusters.size(); clusterId++) {
        const vector<vertex_descriptor>& cluster = clusters[clusterId];
        for(const vertex_descriptor v: cluster) {
            graph[v] = clusterId;
        }
    }
}



void AnchorPair::OrientedReadIdSimilarityGraph::writeGraphviz(
    const string& fileName,
    const vector<OrientedReadId>& orientedReadIds) const
{
    ofstream dot(fileName);
    writeGraphviz(dot, orientedReadIds);
}



void AnchorPair::OrientedReadIdSimilarityGraph::writeGraphviz(
    ostream& dot,
    const vector<OrientedReadId>& orientedReadIds) const
{
    using Graph = OrientedReadIdSimilarityGraph;
    const Graph& graph = *this;

    dot << "graph G {\n";
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const OrientedReadId orientedReadId = orientedReadIds[v];
        const uint64_t clusterId = graph[v];
        const string color = hslToRgbString(double(clusterId) / double(clusters.size()), 0.75, 0.6);
        dot << v << " [label=\"" << orientedReadId << "\\n" << clusterId << "\""
            " fillcolor=\"" << color << "\"];\n";
    }
    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        dot << v0 << "--" << v1 << ";\n";
    }
    dot << "}\n";
}



void AnchorPair::OrientedReadIdSimilarityGraph::writeHtml(
    ostream& html,
    const vector<OrientedReadId>& orientedReadIds) const
{

    // Write it out in Graphviz format.
    const string uuid = "abc"; //to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    writeGraphviz(dotFileName, orientedReadIds);

    // Display it in html in svg format.
    const double timeout = 30.;
    const string options = "-Nshape=rectangle -Nstyle=filled -Goverlap=false -Gsplines=true -Gbgcolor=gray95";
    html << "<h3>Oriented read similarity graph</h3><p>";
    graphvizToHtml(dotFileName, "sfdp", timeout, options, html);
}



// Return the url for the exploreAnchorPair1 page for this AnchorPair.
string AnchorPair::url() const
{
    string s =
        "exploreAnchorPair1?"
        "anchorIdAString=" + HttpServer::urlEncode(anchorIdToString(anchorIdA)) +
        "&anchorIdBString=" + HttpServer::urlEncode(anchorIdToString(anchorIdB)) +
        "&orientedReadIdsString=";

    for(const OrientedReadId orientedReadId: orientedReadIds) {
        s += orientedReadId.getString();
        s += ",";
    }

    return s;
}



void AnchorPair::writeAllHtml(
    ostream& html,
    const Anchors& anchors,
    const Journeys& journeys) const
{
    // Write the summary and oriented reads.
    writeSummaryHtml(html, anchors);
    writeOrientedReadIdsHtml(html, anchors);

    // Get the positions of anchorIdA and anchorIdB in the
    // journeys of all OrientedReadids.
    vector< pair<uint32_t, uint32_t> > positionsInJourneys;
    getPositionsInJourneys(anchors, positionsInJourneys);

    // Write the journey portions between anchorIdA and anchorIdB.
    writeJourneysHtml(html, journeys, positionsInJourneys);

    // Get the internal AnchorIds.
    vector<AnchorId> internalAnchorIds;
    getInternalAnchorIds(journeys, positionsInJourneys, internalAnchorIds);

    // If there are no internalAnchorIds, stop here.
    if(internalAnchorIds.empty()) {
        html << "<p>There are no anchors internal to these journey portions.";
        return;
    }

    // Create the SimpleLocalAnchorGraph.
    // We create it here to get the approximate topological order for AnchorIds,
    // but we display it later.
    SimpleLocalAnchorGraph simpleLocalAnchorGraph(anchors, journeys, *this);

    // Get the internal AnchorIds in topological order;
    vector<AnchorId> internalAnchorIdsInTopologicalOrder;
    simpleLocalAnchorGraph.getInternalAnchorIdsInTopologicalOrder(internalAnchorIdsInTopologicalOrder);

    // Compute the clustering matrix.
    Matrix clusteringMatrix;
    computeClusteringMatrix(journeys, positionsInJourneys, internalAnchorIds, clusteringMatrix);

    // Write the clustering matrix, with the columns written in topological order.
    writeClusteringMatrix(html, internalAnchorIds, internalAnchorIdsInTopologicalOrder, clusteringMatrix);

    // SVD of the clustering matrix.
    vector<double> singularValues;
    Matrix leftSingularVectors;
    Matrix rightSingularVectors;
    clusteringMatrixSvd(clusteringMatrix, singularValues, leftSingularVectors, rightSingularVectors);

    // Write the SVD.
    const uint64_t singularValueCountForOutput = 6;
    writeClusteringMatrixSvd(html,
        internalAnchorIds, internalAnchorIdsInTopologicalOrder,
        singularValueCountForOutput,
        singularValues, leftSingularVectors, rightSingularVectors);

    // Use the scaled left singular vectors to compute a distance matrix between OrientedReadIds.
    const uint64_t singularValueCountForDistanceMatrix = 2;
    Matrix distanceMatrix;
    computeDistanceMatrix(singularValueCountForDistanceMatrix, singularValues, leftSingularVectors, distanceMatrix);

    // Write out the distance matrix.
    writeDistanceMatrixHtml(html, distanceMatrix);

    // Use this distance matrix to compute a similarity graph between oriented reads,
    // then write it out to html.
    const double maxDistance = 1.;
    const OrientedReadIdSimilarityGraph orientedReadIdSimilarityGraph(distanceMatrix, maxDistance);
    orientedReadIdSimilarityGraph.writeHtml(html, orientedReadIds);

    // Also write out the clusters.
    writeClustersHtml(html, orientedReadIdSimilarityGraph);

    // Write the simple local anchor graph.
    writeSimpleLocalAnchorGraphHtml(html, simpleLocalAnchorGraph);

}



void AnchorPair::writeSummaryHtml(ostream& html, const Anchors& anchors) const
{
    const uint64_t offset = getAverageOffset(anchors);
    html <<
        "<table>"
        "<tr><th>Anchor A<td class=centered>" << anchorIdToString(anchorIdA) <<
        "<tr><th>Anchor B<td class=centered>" << anchorIdToString(anchorIdB) <<
        "<tr><th>Coverage<td class=centered>" << size() <<
        "<tr><th>Average offset<td class=centered>" << offset <<
        "</table>";
}



void AnchorPair::writeOrientedReadIdsHtml(ostream& html, const Anchors& anchors) const
{
    vector< pair<Positions, Positions> > positions;
    get(anchors, positions);

    html << "<h3>Oriented reads</h3>";
    html <<
        "<p>"
        "<table>"
        "<tr><th>Oriented<br>read id"
        "<th>Position<br>in journey<br>A<th>Position<br>in journey<br>B<th>Journey<br>offset"
        "<th>OrdinalA<th>OrdinalB<th>Ordinal<br>offset"
        "<th>A middle<br>position"
        "<th>B middle<br>position"
        "<th>Sequence<br>length";

    for(uint64_t i=0; i<size(); i++) {
        const OrientedReadId orientedReadId = orientedReadIds[i];
        const auto& positionsAB = positions[i];

        const auto& positionsA = positionsAB.first;
        const auto& positionsB = positionsAB.second;

        html <<
            "<tr>"
            "<td class=centered>" << orientedReadId <<
            "<td class=centered>" << positionsA.positionInJourney <<
            "<td class=centered>" << positionsB.positionInJourney <<
            "<td class=centered>" << positionsB.positionInJourney - positionsA.positionInJourney <<
            "<td class=centered>" << positionsA.ordinal <<
            "<td class=centered>" << positionsB.ordinal <<
            "<td class=centered>" << positionsB.ordinal - positionsA.ordinal <<
            "<td class=centered>" << positionsA.basePosition <<
            "<td class=centered>" << positionsB.basePosition <<
            "<td class=centered>" << positionsB.basePosition - positionsA.basePosition;
    }
    html << "</table>";

}



void AnchorPair::writeJourneysHtml(
    ostream& html,
    const Journeys& journeys,
    const vector< pair<uint32_t, uint32_t> >& positionsInJourneys   // As computed by getPositionsInJourneys.
    ) const
{
    html << "<h3>Journey portions within this anchor pair</h3>";
    html << "<table>";

    for(uint64_t i=0; i<size(); i++) {
        const OrientedReadId orientedReadId = orientedReadIds[i];
        const Journey& journey = journeys[orientedReadId];

        const auto& positionInJourney = positionsInJourneys[i];
        const auto positionInJourneyA = positionInJourney.first;
        const auto positionInJourneyB = positionInJourney.second;

        if(html) {
            html << "<tr><th class=centered>" << orientedReadId;
        }
        for(auto position=positionInJourneyA+1; position<positionInJourneyB; position++) {
            const AnchorId anchorId = journey[position];
            html << "<td class=centered>" << anchorIdToString(anchorId);
        }
    }

    html << "</table>";

}



void AnchorPair::writeClusteringMatrix(
    ostream& html,
    // The internalAnchorIds as computed by getInternalAnchorIds.
    const vector<AnchorId>& internalAnchorIds,
    // The same AnchorIds, in the order in which the corresponding columns should be written out
    const vector<AnchorId>& internalAnchorIdsInOutputOrder,
    const Matrix& clusteringMatrix) const
{

    html <<
        "<h3>Clustering matrix</h3>"
        "<p><table>"
        "<tr><th>Oriented<br>read<br>id";
    for(const AnchorId anchorId: internalAnchorIdsInOutputOrder) {
        html << "<th>" << anchorIdToString(anchorId);
    }


    // Each OrientedReadId generates a row of the table.
    for(uint64_t i=0; i<size(); i++) {
        const OrientedReadId orientedReadId = orientedReadIds[i];
        html << "<tr><th>" << orientedReadId;

        // Loop over the AnchorIds in topological order.
        for(const AnchorId anchorId: internalAnchorIdsInOutputOrder) {

            // Find the index of this anchor in internalAnchorIds.
            // This is also the column index of the corresponding column
            // in the clustering matrix.
            auto it = std::ranges::lower_bound(internalAnchorIds, anchorId);
            SHASTA_ASSERT(it != internalAnchorIds.end());
            SHASTA_ASSERT(*it == anchorId);
            const int64_t j = it - internalAnchorIds.begin();

            html << "<td class=centered";
            if(clusteringMatrix(i, j) > 0.) {
                html << " style='background-color:LightGreen'>1";
            } else {
                html << ">0";
            }
        }
    }



    html << "</table>";

}



void AnchorPair::writeClusteringMatrixSvd(
    ostream& html,
    // The internalAnchorIds as computed by getInternalAnchorIds.
    const vector<AnchorId>& internalAnchorIds,
    // The same AnchorIds, in the order in which the corresponding columns should be written out
    const vector<AnchorId>& internalAnchorIdsInOutputOrder,
    uint64_t singularValueCount, // Number of singular values/vectors to be written,
    const vector<double>& singularValues,
    const Matrix& leftSingularVectors,
    const Matrix& rightSingularVectors) const
{
    // The actual number of singular values/vectors we will write out.
    singularValueCount = min(singularValueCount, singularValues.size());

    html << std::fixed << std::setprecision(3);

    // Singular values.
    html <<
        "<h3>Singular values</h3><table>";
    for(uint64_t j=0; j<singularValueCount; j++) {
        html << "<tr><th>S<sub>" << j << "</sub><td class=centered>" << singularValues[j];
    }
    html << "</table>";



    // Left singular vectors.
    html <<
        "<h3>Left singular vectors</h3><table>"
        "<tr><th>Oriented<br>read<br>id";
    for(uint64_t j=0; j<singularValueCount; j++) {
        html << "<th>L<sub>" << j << "</sub>";
    }
    for(uint64_t j=0; j<singularValueCount; j++) {
        html << "<th>S<sub>" << j << "</sub>L<sub>" << j << "</sub>";
    }
    html << "\n";
    for(uint64_t i=0; i<size(); i++) {
        const OrientedReadId orientedReadId = orientedReadIds[i];
        html << "<tr><th>" << orientedReadId;

        // Without scaling.
        for(uint64_t j=0; j<singularValueCount; j++) {
            html << "<td class=centered>" << leftSingularVectors(i, j);
        }

        // With scaling.
        for(uint64_t j=0; j<singularValueCount; j++) {
            html << "<td class=centered>" << singularValues[j] * leftSingularVectors(i, j);
        }

        html << "\n";
    }
    html << "</table>";



    // Right singular vectors.
    // The anchors are written in the order specified by internalAnchorIdsInOutputOrder.
    html <<
        "<h3>Right singular vectors</h3><table>"
        "<tr><th>Anchor";
    for(uint64_t j=0; j<singularValueCount; j++) {
        html << "<th>R<sub>" << j << "</sub>";
    }
    for(uint64_t j=0; j<singularValueCount; j++) {
        html << "<th>S<sub>" << j << "</sub>R<sub>" << j << "</sub>";
    }
    html << "\n";
    for(const AnchorId anchorId: internalAnchorIdsInOutputOrder) {

        // Find the index of this AnchorId in internalAnchorIds.
        // This is also the index of the corresponding column in rightSingularValues.
        const auto it = std::ranges::lower_bound(internalAnchorIds, anchorId);
        SHASTA_ASSERT(it != internalAnchorIds.end());
        SHASTA_ASSERT(*it == anchorId);
        const uint64_t j = it - internalAnchorIds.begin();

        html << "<tr><th>" << anchorIdToString(anchorId);

        // Without scaling.
        for(uint64_t i=0; i<singularValueCount; i++) {
            html << "<td class=centered>" << rightSingularVectors(i, j);
        }

        // With scaling.
        for(uint64_t i=0; i<singularValueCount; i++) {
            html << "<td class=centered>" << singularValues[i] * rightSingularVectors(i, j);
        }
        html << "\n";
    }
    html << "</table>";
}



void AnchorPair::writeDistanceMatrixHtml(
    ostream& html,
    const Matrix& distanceMatrix) const
{

    html << "<h3>Distance matrix</h3><table><tr><td>";
    for(const OrientedReadId orientedReadId: orientedReadIds) {
        html << "<th>" << orientedReadId;
    }

    for(uint64_t iA=0; iA<size(); iA++) {
        const OrientedReadId orientedReadIdA = orientedReadIds[iA];
        html << "<tr><th>" << orientedReadIdA;

        for(uint64_t iB=0; iB<size(); iB++) {
            html << "<td class=centered>" << distanceMatrix(iA, iB);
        }
    }

    html << "</table>";
}



void AnchorPair::writeClustersHtml(
    ostream& html, const OrientedReadIdSimilarityGraph& graph) const
{
    html << "<h3>Clusters</h3><table><tr><th>Id<th>Size";
    for(uint64_t clusterId=0; clusterId<graph.clusters.size(); clusterId++) {
        const auto& cluster = graph.clusters[clusterId];
        html << "<tr><th>" << clusterId << "<th>" << cluster.size();
        for(const uint64_t i: cluster) {
            html << "<td class=centered>" << orientedReadIds[i];
        }
    }
    html << "</table>";

}



void AnchorPair::writeSimpleLocalAnchorGraphHtml(
    ostream& html,
    const SimpleLocalAnchorGraph& simpleLocalAnchorGraph) const
{

    // Write it out in Graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    simpleLocalAnchorGraph.writeGraphviz(dotFileName);

    // Display it in html in svg format.
    const double timeout = 30.;
    const string options = "-Nshape=rectangle -Nstyle=filled -Nfillcolor=pink -Gbgcolor=gray95";
    html << "<h3>Local simple anchor graph</h3><p>";
    graphvizToHtml(dotFileName, "dot", timeout, options, html);
}



AnchorPair::SimpleLocalAnchorGraph::SimpleLocalAnchorGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    const AnchorPair& anchorPair)
{
    using Graph = SimpleLocalAnchorGraph;
    Graph& graph = *this;

    // Get the positions of anchorIdA and anchorIdB in the journeys
    // of the OrientedReadIds of this AnchorPair.
    vector< pair<uint32_t, uint32_t> > positionsInJourney;
    anchorPair.getPositionsInJourneys(anchors, positionsInJourney);


    // Create the vertices.
    // Vertex descriptors are the same as indexes into the AnchorIds vector.
    vector<AnchorId> anchorIds;
    vector<uint64_t> localCoverage;
    anchorPair.getAllAnchorIdsAndLocalCoverage(journeys, positionsInJourney, anchorIds, localCoverage);
    SHASTA_ASSERT(anchorIds.size() == localCoverage.size());
    for(uint64_t i=0; i<anchorIds.size(); i++) {
        add_vertex(SimpleLocalAnchorGraphVertex(anchorIds[i], localCoverage[i]), graph);
    }

    // Create the edges.
    for(uint64_t i=0; i<anchorPair.size(); i++) {
        const OrientedReadId orientedReadId = anchorPair.orientedReadIds[i];
        const Journey& journey = journeys[orientedReadId];

        const auto& positionsInJourneyAB = positionsInJourney[i];
        const auto positionInJourneyA = positionsInJourneyAB.first;
        const auto positionInJourneyB = positionsInJourneyAB.second;

        for(auto position=positionInJourneyA+1; position<=positionInJourneyB; position++) {
            const AnchorId anchorId1 = journey[position];
            const AnchorId anchorId0 = journey[position-1];
            const uint64_t i0 = std::ranges::find(anchorIds, anchorId0) - anchorIds.begin();
            const uint64_t i1 = std::ranges::find(anchorIds, anchorId1) - anchorIds.begin();
            Graph::edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = edge(i0, i1, graph);
            if(not edgeExists) {
                tie(e, edgeExists) = boost::add_edge(i0, i1, graph);
                SHASTA_ASSERT(edgeExists);
            }
            ++graph[e].localCoverage;
        }
    }

    approximateTopologicalSort();
    SHASTA_ASSERT(graph[approximateTopologicalOrder.front()].anchorId == anchorPair.anchorIdA);
    SHASTA_ASSERT(graph[approximateTopologicalOrder.back()].anchorId == anchorPair.anchorIdB);
}



void AnchorPair::SimpleLocalAnchorGraph::approximateTopologicalSort()
{
    using Graph = SimpleLocalAnchorGraph;
    Graph& graph = *this;

    // Sort the edges by decreasing coverage.
    vector< pair<Graph::edge_descriptor, uint64_t> > edgeTable;
    BGL_FORALL_EDGES(e, graph, Graph) {
        edgeTable.emplace_back(e, graph[e].localCoverage);
    }
    std::ranges::sort(edgeTable, OrderPairsBySecondOnlyGreater<Graph::edge_descriptor, uint64_t>());
    vector<Graph::edge_descriptor> edgesByCoverage;
    for(const auto& p: edgeTable) {
        edgesByCoverage.push_back(p.first);
    }

    // Do an approximate topological sort using edges in this order.
    shasta::approximateTopologicalSort(graph, edgesByCoverage);
    vector< pair<Graph::vertex_descriptor, uint64_t> > vertexTable;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        vertexTable.emplace_back(v, graph[v].rank);
    }
    std::ranges::sort(vertexTable, OrderPairsBySecondOnly<Graph::vertex_descriptor, uint64_t>());
    approximateTopologicalOrder.clear();
    for(const auto& p: vertexTable) {
        approximateTopologicalOrder.push_back(p.first);
    }
}



void AnchorPair::SimpleLocalAnchorGraph::getInternalAnchorIdsInTopologicalOrder(
    vector<AnchorId>& internalAnchorIdsInTopologicalOrder) const
{
    const SimpleLocalAnchorGraph& graph = *this;

    internalAnchorIdsInTopologicalOrder.clear();

    // We only want internal AnchorIds, so we have to skip the first and last vertex.
    for(uint64_t j=1; j<approximateTopologicalOrder.size()-1; j++) {
        const vertex_descriptor v = approximateTopologicalOrder[j];
        internalAnchorIdsInTopologicalOrder.push_back(graph[v].anchorId);
    }
}



void AnchorPair::SimpleLocalAnchorGraph::writeGraphviz(const string& fileName) const
{
    ofstream dot(fileName);
    writeGraphviz(dot);
}



void AnchorPair::SimpleLocalAnchorGraph::writeGraphviz(ostream& dot) const
{
    using Graph = SimpleLocalAnchorGraph;
    const Graph& graph = *this;

    dot << "digraph SimpleLocalAnchorGraph {\n";
    for(const vertex_descriptor v: approximateTopologicalOrder) {
        const SimpleLocalAnchorGraphVertex& vertex = graph[v];
         dot << v << " [label=\"" << anchorIdToString(vertex.anchorId) << "\\n" << vertex.localCoverage << "\"];\n";
    }
    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const SimpleLocalAnchorGraphEdge& edge = graph[e];
        dot << v0 << "->" << v1 << " [label=\"" << edge.localCoverage << "\"];\n";
    }
    dot << "}\n";

}


void AnchorPair::getAllAnchorIds(
    const Journeys& journeys,
    const vector< pair<uint32_t, uint32_t> >& positionsInJourneys,
    vector<AnchorId>& anchorIds) const
{
    anchorIds.clear();

    // Loop over OrientedReadIds in this AnchorPair.
    for(uint64_t i=0; i<size(); i++) {
        const OrientedReadId orientedReadId = orientedReadIds[i];
        const Journey journey = journeys[orientedReadId];

        // Get the positions in the journet of anchorIdA and anchorIdB.
        const auto& p = positionsInJourneys[i];
        const uint32_t positionInJourneyA = p.first;
        const uint32_t positionInJourneyB = p.second;

        // Loop over this portion of the journey, including anchorIdA and anchorIdB.
        for(uint64_t position=positionInJourneyA; position<=positionInJourneyB; position++) {
            anchorIds.push_back(journey[position]);
        }
    }

    deduplicate(anchorIds);
}



void AnchorPair::getInternalAnchorIds(
    const Journeys& journeys,
    const vector< pair<uint32_t, uint32_t> >& positionsInJourneys,
    vector<AnchorId>& anchorIds) const
{
    anchorIds.clear();

    // Loop over OrientedReadIds in this AnchorPair.
    for(uint64_t i=0; i<size(); i++) {
        const OrientedReadId orientedReadId = orientedReadIds[i];
        const Journey journey = journeys[orientedReadId];

        // Get the positions in the journey of anchorIdA and anchorIdB.
        const auto& p = positionsInJourneys[i];
        const uint32_t positionInJourneyA = p.first;
        const uint32_t positionInJourneyB = p.second;

        // Loop over this portion of the journey, excluding anchorIdA and anchorIdB.
        for(uint64_t position=positionInJourneyA+1; position<positionInJourneyB; position++) {
            anchorIds.push_back(journey[position]);
        }
    }

    deduplicate(anchorIds);
}



void AnchorPair::getAllAnchorIdsAndLocalCoverage(
    const Journeys& journeys,
    const vector< pair<uint32_t, uint32_t> >& positionsInJourneys,
    vector<AnchorId>& anchorIds,
    vector<uint64_t>& localCoverage) const
{
    // Loop over OrientedReadIds in this AnchorPair.
    anchorIds.clear();
    for(uint64_t i=0; i<size(); i++) {
        const OrientedReadId orientedReadId = orientedReadIds[i];
        const Journey journey = journeys[orientedReadId];

        // Get the positions in the journet of anchorIdA and anchorIdB.
        const auto& p = positionsInJourneys[i];
        const uint32_t positionInJourneyA = p.first;
        const uint32_t positionInJourneyB = p.second;
        for(uint64_t position=positionInJourneyA; position<=positionInJourneyB; position++) {
            anchorIds.push_back(journey[position]);
        }
    }

    // Deduplicate and count.
    deduplicateAndCount(anchorIds, localCoverage);

}
