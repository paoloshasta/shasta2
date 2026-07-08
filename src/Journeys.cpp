#include "Journeys.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "ReadId.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace shasta2;

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<Journeys>;



// Initial creation.
// This sets the positionInJourney for every AnchorMarkerInfo
// stored in the Anchors, and for this reason the Anchors
// are not passed in as const.
Journeys::Journeys(
    uint64_t orientedReadCount,
    shared_ptr<Anchors> anchorsPointer,
    uint64_t threadCount,
    const MappedMemoryOwner& mappedMemoryOwner) :
    MultithreadedObject<Journeys>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    anchorsPointer(anchorsPointer)

{
    performanceLog << timestamp << "Journeys creation begins." << endl;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t anchorCount = anchorsPointer->size();
    const uint64_t anchorBatchCount = 1000;
    const uint64_t orientedReadBatchCount = 1000;

    // Pass1: make space for the journeysWithPositions.
    journeysWithPositions.createNew(largeDataName("tmp-JourneysWithPositions"), largeDataPageSize);
    journeysWithPositions.beginPass1(orientedReadCount);
    setupLoadBalancing(anchorCount, anchorBatchCount);
    runThreads(&Journeys::threadFunction1, threadCount);

    // Pass2: store the unsorted journeysWithPositions.
    journeysWithPositions.beginPass2();
    setupLoadBalancing(anchorCount, anchorBatchCount);
    runThreads(&Journeys::threadFunction2, threadCount);
    journeysWithPositions.endPass2();
    journeysWithPositions.unreserve();

    // Pass 3:sort the journeysWithPostions and make space for the journeys.
    journeys.createNew(largeDataName("Journeys"), largeDataPageSize);
    journeys.beginPass1(orientedReadCount);
    setupLoadBalancing(orientedReadCount, orientedReadBatchCount);
    runThreads(&Journeys::threadFunction3, threadCount);

    // Pass 4: copy the sorted journeysWithPositions to the journeys.
    journeys.beginPass2();
    setupLoadBalancing(orientedReadCount, orientedReadBatchCount);
    runThreads(&Journeys::threadFunction4, threadCount);
    journeys.endPass2(false, true);
    journeys.unreserve();

    journeysWithPositions.remove();

    performanceLog << timestamp << "Journeys creation ends." << endl;
}



void Journeys::threadFunction1(uint64_t /* threadId */)
{
    threadFunction12(1);
}



void Journeys::threadFunction2(uint64_t /* threadId */)
{
    threadFunction12(2);
}



void Journeys::threadFunction12(uint64_t pass)
{
    const Anchors& anchors = *anchorsPointer;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all AnchorIds in this batch.
        for(AnchorId anchorId=begin; anchorId!=end; anchorId++) {
            Anchor anchor = anchors[anchorId];

            // Loop over the marker intervals of this Anchor.
            for(const auto& anchorMarkerInterval: anchor) {
                const auto orientedReadIdValue = anchorMarkerInterval.orientedReadId.getValue();

                if(pass == 1) {
                    journeysWithPositions.incrementCountMultithreaded(orientedReadIdValue);
                } else {
                    journeysWithPositions.storeMultithreaded(
                        orientedReadIdValue, {anchorId, anchorMarkerInterval.position});
                }
            }
        }
    }
}



void Journeys::threadFunction3(uint64_t /* threadId */)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all oriented reads assigned to this thread.
        for(uint64_t orientedReadValue=begin; orientedReadValue!=end; orientedReadValue++) {
            auto v = journeysWithPositions[orientedReadValue];
            sort(v.begin(), v.end(), OrderPairsBySecondOnly<uint64_t, uint32_t>());
            journeys.incrementCountMultithreaded(orientedReadValue, v.size());
        }
    }
}



void Journeys::threadFunction4(uint64_t /* threadId */)
{
    Anchors& anchors = *anchorsPointer;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all oriented reads assigned to this thread.
        for(uint64_t orientedReadValue=begin; orientedReadValue!=end; orientedReadValue++) {
            const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(orientedReadValue));

            // Check that the positions are monotonically increasing.
            const auto v = journeysWithPositions[orientedReadValue];
            for(uint64_t i1=1; i1<v.size(); i1++) {
                const uint64_t i0 = i1 - 1;
                if(not (v[i0].second < v[i1].second)) {
                    throw runtime_error(
                        "Invalid Journey for oriented read " + orientedReadId.getString() +
                        " detected near anchors " +
                        anchorIdToString(v[i0].first) + " " + anchorIdToString(v[i1].first));
                }
            }

            // Copy the journeysWithPositions to the journeys.
            const auto journey = journeys[orientedReadValue];
            SHASTA2_ASSERT(journey.size() == v.size());
            for(uint64_t i=0; i<v.size(); i++) {
                journey[i] = v[i].first;
            }

            // Store journey information for this oriented read in the marker interval.
            for(uint64_t position=0; position<journey.size(); position++) {
                const AnchorId anchorId = journey[position];
                span<AnchorMarkerInfo> markerInfos = anchors.anchorMarkerInfos[anchorId];
                bool found = false;
                for(AnchorMarkerInfo& markerInfo: markerInfos) {
                    if(markerInfo.orientedReadId == orientedReadId) {
                        markerInfo.positionInJourney = uint32_t(position);
                        found = true;
                        break;
                    }
                }
                SHASTA2_ASSERT(found);
            }
        }
    }
}



// Access from binary data.
Journeys::Journeys(const MappedMemoryOwner& mappedMemoryOwner) :
    MultithreadedObject<Journeys>(*this),
    MappedMemoryOwner(mappedMemoryOwner)
{
    journeys.accessExistingReadOnly(largeDataName("Journeys"));
}
