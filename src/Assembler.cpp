#include "Assembler.hpp"
#include "Anchor.hpp"
#include "AnchorGraph.hpp"
#include "AssemblyGraph.hpp"
#include "deduplicate.hpp"
#include "Journeys.hpp"
#include "KmerCheckerFactory.hpp"
#include "Markers.hpp"
#include "MarkerKmers.hpp"
#include "MurmurHash2.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
#include "ReadGraph.hpp"
#include "Reads.hpp"
#include "ReadSummary.hpp"
using namespace shasta2;

#include "MultithreadedObject.tpp"
template class MultithreadedObject<Assembler>;


// Construct a new Assembler.
Assembler::Assembler(
    const string& largeDataFileNamePrefix,
    size_t largeDataPageSize) :
    MultithreadedObject(*this),
    MappedMemoryOwner(largeDataFileNamePrefix, largeDataPageSize)
{


    assemblerInfo.createNew(largeDataName("Info"), largeDataPageSize);
    assemblerInfo->largeDataPageSize = largeDataPageSize;

    readsPointer = make_shared<Reads>();
    readsPointer->createNew(
        largeDataName("Reads"),
        largeDataName("ReadNames"),
        largeDataName("ReadIdsSortedByName"),
        largeDataPageSize
    );
}



// Construct an Assembler from binary data. This accesses the AssemblerInfo and the Reads.
Assembler::Assembler(const string& largeDataFileNamePrefix) :
    MultithreadedObject(*this),
    MappedMemoryOwner(largeDataFileNamePrefix, 0)
{

    assemblerInfo.accessExistingReadWrite(largeDataName("Info"));
    largeDataPageSize = assemblerInfo->largeDataPageSize;

    readsPointer = make_shared<Reads>();
    readsPointer->access(
        largeDataName("Reads"),
        largeDataName("ReadNames"),
        largeDataName("ReadIdsSortedByName")
    );

}



// This runs the entire assembly, under the following assumptions:
// - The current directory is the run directory.
// - The Data directory has already been created and set up, if necessary.
// - The input file names are either absolute,
//   or relative to the run directory, which is the current directory.
void Assembler::assemble(
    const Options& options,
    vector<string> inputFileNames)
{
    cout << "Number of threads: " << options.threadCount << endl;

    addReads(
        inputFileNames,
        options.minReadLength,
        options.threadCount);
    createReadSummaries();

    createKmerChecker(options.k, options.markerDensity);
    createMarkers(options.threadCount);

    findPalindromicReads();

    createMarkerKmers(options.maxMarkerErrorRate, options.threadCount);

    if(options.externalAnchorsName.empty()) {
        createAnchors(
            options.minAnchorCoverage,
            options.maxAnchorCoverage,
            options.maxAnchorRepeatLength,
            options.minAnchorDistinctSubkmerCount,
            options.threadCount);
    } else {
        readExternalAnchors(options.externalAnchorsName);
    }

    createJourneys(options.threadCount);
    storeAnchorGaps();

    createAnchorGraph(options);
    anchorGraphTransitiveReduction(options);
    saveAnchorGraph();

    createAssemblyGraph(options);

    writeReadSummaries();
}





void Assembler::createKmerChecker(
    uint64_t k,
    double markerDensity)
{
    assemblerInfo->k = k;
    assemblerInfo->markerDensity = markerDensity;
    kmerChecker = KmerCheckerFactory::createNew(
        k,
        markerDensity);
}



// Generate Anchors from MarkerKmers.
void Assembler::createAnchors(
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage,
    const vector<uint64_t>& maxAnchorRepeatLength,
    const vector<uint64_t>& minAnchorDistinctSubkmerCount,
    uint64_t threadCount)
{
    const bool doAnchorCleanup = false;

    if(doAnchorCleanup) {

        // Generate the initial Anchors from the MarkerKmers.
        shared_ptr<Anchors> initialAnchors = make_shared<Anchors>(
            "InitialAnchors",
            MappedMemoryOwner(*this),
            reads(),
            assemblerInfo->k,
            markers(),
            *markerKmers,
            minAnchorCoverage,
            maxAnchorCoverage,
            maxAnchorRepeatLength,
            minAnchorDistinctSubkmerCount,
            threadCount);
        anchorsPointer = initialAnchors;

        // The Anchor cleanup process will create the final Anchors by removing
        // some AnchorMarkerInfos from the initial anchors.
        // The ReadGraph decides which AnchorMarkerInfos should be kept.
        vector<bool> keep;
        createReadGraph(threadCount);
        readGraph().anchorCleanup(keep);

        // Create the new anchors.
        shared_ptr<Anchors> finalAnchors = make_shared<Anchors>(anchors(), "Anchors", keep);

        cout << "Before anchor cleanup, there were " << initialAnchors->size() <<
            " anchors with a total " << initialAnchors->anchorMarkerInfos.totalSize() <<
            " AnchorMarkerInfos." << endl;
        cout << "After anchor cleanup, there are " << finalAnchors->size() <<
            " anchors with a total " << finalAnchors->anchorMarkerInfos.totalSize() <<
            " AnchorMarkerInfos." << endl;

        // Keep the final anchors.
       anchorsPointer = finalAnchors;

       // Remove the initial anchors.
       initialAnchors->remove();
       initialAnchors = 0;
    }



    // Simple Anchor generation without cleanup.
    else {
       anchorsPointer = make_shared<Anchors>(
            "Anchors",
            MappedMemoryOwner(*this),
            reads(),
            assemblerInfo->k,
            markers(),
            *markerKmers,
            minAnchorCoverage,
            maxAnchorCoverage,
            maxAnchorRepeatLength,
            minAnchorDistinctSubkmerCount,
            threadCount);
    }

    cout << "There are " << anchors().size() << " anchors." << endl;
}



// Read Anchors from ExternalAnchors.
void Assembler::readExternalAnchors(const string& externalAnchorsName)
{
    anchorsPointer = make_shared<Anchors>(
        "Anchors",
        MappedMemoryOwner(*this),
        reads(),
        assemblerInfo->k,
        markers(),
        *markerKmers,
        externalAnchorsName);
}



// Access existing Anchors.
void Assembler::accessAnchors(bool writeAccess)
{
     anchorsPointer = make_shared<Anchors>("Anchors",
         MappedMemoryOwner(*this), reads(), assemblerInfo->k, markers(), *markerKmers, writeAccess);
}



void Assembler::createJourneys(uint64_t threadCount)
{
    const MappedMemoryOwner& mappedMemoryOwner = *this;

    journeysPointer = make_shared<Journeys>(
        2 * reads().readCount(),
        anchorsPointer,
        threadCount,
        mappedMemoryOwner);

}



void Assembler::accessJourneys()
{
    journeysPointer = make_shared<Journeys>(*this);
}



// Store anchor gaps information in ReadSummary for each read.
void Assembler::storeAnchorGaps()
{
    const uint32_t kHalf = uint32_t(anchors().kHalf);

    // Loop over all Reads.
    for(ReadId readId=0; readId<reads().readCount(); readId++) {
        ReadSummary& readSummary = readSummaries[readId];
        const uint32_t readLength = uint32_t(reads().getReadSequenceLength(readId));

        // Put it on strand 0.
        const OrientedReadId orientedReadId(readId, 0);

        // Get the markers and the journey of this oriented read.
        const auto orientedReadMarkers = markers()[orientedReadId.getValue()];
        const auto journey = journeys()[orientedReadId];

        if(journey.empty()) {
            readSummary.initialAnchorGap = readLength;
            readSummary.middleAnchorGap = readLength;
            readSummary.finalAnchorGap = readLength;
            continue;
        }

        // Compute the largest gap between adjacent anchors on the journey.
        uint32_t maxGap = 0;
        for(uint64_t i1=1; i1<journey.size(); i1++) {
            const uint64_t i0 = i1 - 1;

            const AnchorId anchorId0 = journey[i0];
            const AnchorId anchorId1 = journey[i1];

            const uint32_t ordinal0 = anchors().getOrdinal(anchorId0, orientedReadId);
            const uint32_t ordinal1 = anchors().getOrdinal(anchorId1, orientedReadId);

            const uint32_t position0 = orientedReadMarkers[ordinal0].position + kHalf;
            const uint32_t position1 = orientedReadMarkers[ordinal1].position + kHalf;

            const uint32_t gap = position1 - position0;
            maxGap = max(maxGap, gap);
        }
        readSummary.middleAnchorGap = maxGap;

        // Compute the number of bases preceding the first anchor on the journey.
        const AnchorId anchorId0 = journey.front();
        const uint32_t ordinal0 = anchors().getOrdinal(anchorId0, orientedReadId);
        readSummary.initialAnchorGap = orientedReadMarkers[ordinal0].position + kHalf;

        // Compute the number of bases following the last anchor on the journey.
        const AnchorId anchorId1 = journey.back();
        const uint32_t ordinal1 = anchors().getOrdinal(anchorId1, orientedReadId);
        readSummary.finalAnchorGap = readLength - orientedReadMarkers[ordinal1].position - kHalf;

    }

}



void Assembler::createAnchorGraph(const Options& options)
{
    anchorGraphPointer = make_shared<AnchorGraph>(
        anchors(), journeys(),
        options.minAnchorGraphEdgeCoverage);
}



void Assembler::createCompleteAnchorGraph()
{
    completeAnchorGraphPointer = make_shared<AnchorGraph>(anchors(), journeys(), 1);
    completeAnchorGraphPointer->save("CompleteAnchorGraph");
}



void Assembler::createAssemblyGraph(const Options& options)
{
    writePerformanceStatistics("Assembler::createAssemblyGraph begins");

    AssemblyGraph assemblyGraph(
        anchors(),
        journeys(),
        *anchorGraphPointer,
        options);
    assemblyGraph.simplifyAndAssemble();

    writePerformanceStatistics("Assembler::createAssemblyGraph ends");
}



void Assembler::accessAnchorGraph()
{
    const MappedMemoryOwner& mappedMemoryOwner = *this;
    anchorGraphPointer = make_shared<AnchorGraph>(mappedMemoryOwner, "AnchorGraph");
}


void Assembler::saveAnchorGraph()
{
    anchorGraphPointer->save("AnchorGraph");
}



void Assembler::anchorGraphTransitiveReduction(
    const Options& options)
{
    anchorGraphPointer->transitiveReduction(
        options.transitiveReductionMaxEdgeCoverage,
        options.transitiveReductionMaxDistance);
}



void Assembler::accessCompleteAnchorGraph()
{
    const MappedMemoryOwner& mappedMemoryOwner = *this;
    completeAnchorGraphPointer = make_shared<AnchorGraph>(mappedMemoryOwner, "CompleteAnchorGraph");
}



void Assembler::createReadSummaries()
{
    readSummaries.createNew(largeDataName("ReadSummaries"), largeDataPageSize);
    readSummaries.resize(reads().readCount());
}



void Assembler::accessReadSummaries()
{
    readSummaries.accessExistingReadWrite(largeDataName("ReadSummaries"));
}



void Assembler::writeReadSummaries() const
{
    ofstream csv("ReadSummary.csv");
    csv <<
        "ReadId,"
        "Length,"
        "Use for assembly,"
        "Is palindromic,"
        "Has high error rare,"
        "Palindromic rate,"
        "Initial marker error rate,"
        "Marker error rate,"
        "Initial anchor gap,"
        "Middle anchor gap,"
        "Final anchor gap,"
        "\n";

    for(ReadId readId=0; readId<readSummaries.size(); readId++) {
        const ReadSummary& readSummary = readSummaries[readId];

        csv <<
            readId << "," <<
            reads().getReadSequenceLength(readId) << "," <<
            (readSummary.isInUse() ? "Yes" : "No") << "," <<
            (readSummary.isPalindromic ? "Yes" : "No") << "," <<
            (readSummary.hasHighErrorRate ? "Yes" : "No") << "," <<
            readSummary.palindromicRate << "," <<
            readSummary.initialMarkerErrorRate << "," <<
            readSummary.markerErrorRate << "," <<
            readSummary.initialAnchorGap << "," <<
            readSummary.middleAnchorGap << "," <<
            readSummary.finalAnchorGap << "," <<
            "\n";
    }
}



void Assembler::createReadGraph(uint64_t threadCount)
{
    readGraphPointer = make_shared<ReadGraph>(anchors(), threadCount);
}



void Assembler::accessReadGraph()
{
    readGraphPointer = make_shared<ReadGraph>(anchors());
}
