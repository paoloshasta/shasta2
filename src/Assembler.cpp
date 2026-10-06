#include "Assembler.hpp"
#include "Anchor.hpp"
#include "AnchorGraph.hpp"
#include "AssemblyGraph.hpp"
#include "deduplicate.hpp"
#include "HomopolymerModel.hpp"
#include "HomopolymerModelTable.hpp"
#include "Journeys.hpp"
#include "KmerCheckerFactory.hpp"
#include "Markers.hpp"
#include "MarkerKmers.hpp"
#include "memoryInformation.hpp"
#include "MurmurHash2.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
#include "ReadLoader.hpp"
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
    const vector<string>& inputFileNames,
    const string& externalAnchorsNameAbsolutePath,
    const string& externalAnchorGraphNameAbsolutePath)
{
    cout << "Number of threads: " << options.threadCount << endl;

    // Create the HomopolymerModel.
    createHomopolymerModel(options.homopolymerModelName);

    // Load the reads.
    addReads(
        inputFileNames,
        options.minReadLength,
        options.threadCount);
    createReadSummaries();
    writeReadSummaries(true);



    // Generates Anchors from MarkerKmers or from ExternalAnchors.
    if(externalAnchorsNameAbsolutePath.empty()) {

        // Create the Markers.
        createKmerChecker(options.k, options.markerDensity);
        createMarkers(options.threadCount);
        kmerChecker = 0;

        // Flag palindromic reads. They will be excluded from the rest of
        // the assembly process.
        findPalindromicReadsMultithreaded(options.threadCount);

        // Create the MarkerKmers.
        createMarkerKmers(options.maxMarkerErrorRate, options.threadCount);
        removeIfAllowed(options, *markersPointer);
        markersPointer = 0;

        // Create the Anchors.
        createAnchors(
            options.minAnchorCoverage,
            options.maxAnchorCoverage,
            options.maxAnchorRepeatLength,
            options.minAnchorDistinctSubkmerCount,
            options.threadCount);
        removeIfAllowed(options, *markerKmers);
        markerKmers = 0;

    } else {

        assemblerInfo->k = options.k;
        readExternalAnchors(externalAnchorsNameAbsolutePath);
    }



    // Create the Journeys.
    createJourneys(options.threadCount);
    storeAnchorGaps();



    // Create the AnchorGraph using the Journeys or read it in.
    if(externalAnchorsNameAbsolutePath.empty() or externalAnchorGraphNameAbsolutePath.empty()) {
        createAnchorGraph(options);
    } else {
        accessAnchorGraph(externalAnchorGraphNameAbsolutePath);
    }
    anchorGraphTransitiveReduction(options);
    if((options.memoryMode == "filesystem") and options.keepBinaryData) {
        saveAnchorGraph();
    }



    // Create the AssemblyGraph.
    createAssemblyGraph(options, true);

    writeReadSummaries(false);
}




void Assembler::createHomopolymerModel(const string& homopolymerModelName)
{
    // If the homopolymerModelName is empty, homopolymerModelPointer is set to 0,
    // and a median scheme is used to determine homopolymer lengths
    // instead of a homopolymer model.
    if(homopolymerModelName.empty()) {
        homopolymerModelPointer = 0;
        cout << "Not using a homopolymer model." << endl;
        return;
    }

    // If homopolymerModelName is in the homopolymerModelTable, create the model from there.
    const auto it = homopolymerModelTable.find(homopolymerModelName);
    if(it != homopolymerModelTable.end()) {
        const string& csvString = it->second;
        std::istringstream csv(csvString);
        try {
            homopolymerModelPointer = make_shared<const HomopolymerModel>(csv);
            cout << "Using homopolymer model " << homopolymerModelName << endl;
        } catch(std::exception& e) {
            cout << e.what() << endl;
            throw runtime_error("The above error occurred while reading homopolymer model " + homopolymerModelName);
        }
        return;
    }



    // If getting here, the homopolymerModelName was not in homopolymerModelTable.
    // In this case, homopolymerModelName must be an absolute path to the csv
    // file that defines the Homopolymer model.
    if(homopolymerModelName[0] != '/') {
        throw runtime_error(
            "Options --homopolymer-model must specify an absolute path but the following was used: " +
            homopolymerModelName);
    }

    // Open the csv file that defines the HomopolymerModel.
    ifstream file(homopolymerModelName);
    if(not file) {
        throw runtime_error("Could not open " + homopolymerModelName);
    }

    // Create the homopolymer model.
    try {
        homopolymerModelPointer = make_shared<const HomopolymerModel>(file);
        cout << "Using homopolymer model " << homopolymerModelName << endl;
    } catch(std::exception& e) {
        cout << e.what() << endl;
        throw runtime_error("The above error occurred while reading homopolymer model " + homopolymerModelName);
    }
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
    anchorsPointer = make_shared<Anchors>(
        "Anchors",
        MappedMemoryOwner(*this),
        reads(),
        assemblerInfo->k,
        *markerKmers,
        minAnchorCoverage,
        maxAnchorCoverage,
        maxAnchorRepeatLength,
        minAnchorDistinctSubkmerCount,
        threadCount);
}



// Read Anchors from ExternalAnchors.
void Assembler::readExternalAnchors(const string& externalAnchorsName)
{
    anchorsPointer = make_shared<Anchors>(
        "Anchors",
        MappedMemoryOwner(*this),
        reads(),
        assemblerInfo->k,
        externalAnchorsName);
}



// Access existing Anchors.
void Assembler::accessAnchors()
{
     anchorsPointer = make_shared<Anchors>("Anchors",
         MappedMemoryOwner(*this), reads(), assemblerInfo->k);
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

    // Loop over all Reads.
    for(ReadId readId=0; readId<reads().readCount(); readId++) {
        ReadSummary& readSummary = readSummaries[readId];
        const uint32_t readLength = uint32_t(reads().getReadSequenceLength(readId));

        // Put it on strand 0.
        const OrientedReadId orientedReadId(readId, 0);

        // Get the markers and the journey of this oriented read.
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

            const uint32_t position0 = anchors().getPosition(anchorId0, orientedReadId);
            const uint32_t position1 = anchors().getPosition(anchorId1, orientedReadId);

            const uint32_t gap = position1 - position0;
            maxGap = max(maxGap, gap);
        }
        readSummary.middleAnchorGap = maxGap;

        // Compute the number of bases preceding the first anchor on the journey.
        const AnchorId anchorId0 = journey.front();
        readSummary.initialAnchorGap = anchors().getPosition(anchorId0, orientedReadId);

        // Compute the number of bases following the last anchor on the journey.
        const AnchorId anchorId1 = journey.back();
        readSummary.finalAnchorGap = readLength - anchors().getPosition(anchorId1, orientedReadId);

    }

}



void Assembler::createAnchorGraph(const Options& options)
{
    anchorGraphPointer = make_shared<AnchorGraph>(
        anchors(), journeys(),
        options.minAnchorGraphEdgeCoverage,
        options.minAnchorGraphEdgeCoverageFraction);
}



void Assembler::createCompleteAnchorGraph()
{
    completeAnchorGraphPointer = make_shared<AnchorGraph>(anchors(), journeys(), 1, 0.);
    completeAnchorGraphPointer->save("CompleteAnchorGraph");
}




void Assembler::createAssemblyGraph(const Options& options, bool removeAnchorGraph)
{
    writeMemoryStatistics("Assembler::createAssemblyGraph begins");

    AssemblyGraph assemblyGraph(
        anchors(),
        journeys(),
        *anchorGraphPointer,
        options,
        homopolymerModelPointer);

    writeMemoryStatistics("Before removing AnchorGraph");

    if(removeAnchorGraph) {
        anchorGraphPointer = 0;
    }

    assemblyGraph.simplifyAndAssemble();

    writeMemoryStatistics("Assembler::createAssemblyGraph ends");
}



void Assembler::accessAnchorGraph(string name)
{
    const MappedMemoryOwner& mappedMemoryOwner = *this;
    if(name.empty()) {
        name = largeDataName("AnchorGraph");
    } else {
        cout << "Loading AnchorGraph from " << name << endl;
    }
    anchorGraphPointer = make_shared<AnchorGraph>(mappedMemoryOwner, name);
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
    completeAnchorGraphPointer = make_shared<AnchorGraph>(mappedMemoryOwner,
        largeDataName("CompleteAnchorGraph"));
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



void Assembler::writeReadSummaries(bool partial) const
{
    ofstream csv("ReadSummary.csv");
    csv <<
        "ReadId,"
        "Name,"
        "Length,";
    if(not partial) {
        csv <<
            "Use for assembly,"
            "Is palindromic,"
            "Has high error rare,"
            "Palindromic rate,"
            "Initial marker error rate,"
            "Marker error rate,"
            "Initial anchor gap,"
            "Middle anchor gap,"
            "Final anchor gap,";
        }
    csv << "\n";

    for(ReadId readId=0; readId<readSummaries.size(); readId++) {
        const ReadSummary& readSummary = readSummaries[readId];

        csv <<
            readId << "," <<
            reads().getReadName(readId) << "," <<
            reads().getReadSequenceLength(readId) << ",";
        if(not partial) {
            csv <<
                (readSummary.isInUse() ? "Yes" : "No") << "," <<
                (readSummary.isPalindromic ? "Yes" : "No") << "," <<
                (readSummary.hasHighErrorRate ? "Yes" : "No") << "," <<
                readSummary.palindromicRate << "," <<
                readSummary.initialMarkerErrorRate << "," <<
                readSummary.markerErrorRate << "," <<
                readSummary.initialAnchorGap << "," <<
                readSummary.middleAnchorGap << "," <<
                readSummary.finalAnchorGap << ",";
        }
        csv << "\n";
    }
}



void Assembler::analyzeAnchors(const Options& options) const
{
    SHASTA2_ASSERT(anchorsPointer);
    SHASTA2_ASSERT(journeysPointer);
    const uint64_t k = assemblerInfo->k;

    const Anchors& anchors = *anchorsPointer;
    const Journeys& journeys = *journeysPointer;

    cout << "Assembler::analyzeAnchors begins." << endl;



    // Read a reference to be used for this analysis.
    cout << "Loading the reference." << endl;
    Reads reference;
    reference.createNew(
        largeDataName("Reference"),
        largeDataName("ReferenceNames"),
        largeDataName("ReferenceIdsSortedByName"),
        largeDataPageSize
    );
    ReadLoader readLoader(
        "reference.fasta", 0, options.actualThreadCount(),
        largeDataFileNamePrefix, largeDataPageSize,
        reference);



    // Gather all k-mers present in the reference.
    // Store both strands.
    vector<Kmer> referenceKmers;
    cout << "Gathering reference k-mers." << endl;
    for(uint32_t i=0; i<reference.readCount(); i++) {
        const LongBaseSequenceView referenceContig = reference.getRead(i);
        if(referenceContig.baseCount < k) {
            continue;
        }

        // Loop over k-mers of this reference contig.
        Kmer kmer;
        for(size_t position=0; position<k; position++) {
            kmer.set(position, referenceContig[position]);
        }
        for(uint32_t position=0; position+k < referenceContig.baseCount; position++) {
            referenceKmers.push_back(kmer);
            referenceKmers.push_back(kmer.reverseComplement(k));

            // Update the k-mer.
            kmer.shiftLeft();
            kmer.set(k-1, referenceContig[position+k]);
        }
    }

    cout << "Deduplicating reference k-mers." << endl;
    vector<uint64_t> kmerFrequency;
    deduplicateAndCount(referenceKmers, kmerFrequency);


    class Key {
    public:
        uint64_t maxForwardCoverage;
        uint64_t maxBackwardCoverage;
        auto operator<=>(const Key&) const = default;
    };
    class Value {
    public:
        uint64_t inReferenceCount = 0;
        uint64_t notInReferenceCount = 0;
        auto operator<=>(const Value&) const = default;
    };
    std::map<Key, Value> m;



    // Loop over positive (even) anchors.
    cout << "Analyzing anchors." << endl;
    vector<AnchorId> nextOrPrevious;
    vector<uint64_t> count;
    ofstream csv("AnalyzeAnchors.csv");
    csv << "AnchorId,Coverage,K-mer,Frequency in reference,Max forward coverage,Max backward coverage,\n";
    for(AnchorId anchorId=0; anchorId<anchors.size(); anchorId+=2) {
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

        // Get the Kmer and find out how many times it is present in the reference (both strands).
        const Kmer kmer = anchors.anchorKmer(anchorId);
        const auto it = std::lower_bound(referenceKmers.begin(), referenceKmers.end(), kmer);
        uint64_t frequency = 0;
        if(it != referenceKmers.end()) {
            if(*it == kmer) {
                frequency = kmerFrequency[it - referenceKmers.begin()];
            }
        }

        // Update our map.
        const Key key = Key({maxForwardCoverage, maxBackwardCoverage});
        auto jt = m.find(key);
        if(jt == m.end()) {
            tie(jt, ignore) = m.insert({key, Value()});
        }
        Value& value = jt->second;
        if(frequency == 0) {
            ++value.notInReferenceCount;
        } else {
            ++value.inReferenceCount;
        }


        csv << anchorIdToString(anchorId) << ",";
        csv << anchor.size() << ",";
        kmer.write(csv, k);
        csv << ",";
        csv << frequency << ",";
        csv << maxForwardCoverage << ",";
        csv << maxBackwardCoverage << ",";
        csv << "\n";
    }
    reference.remove();



    // Write out the map.
    {
        ofstream csv("AnalyzeAnchors-Summary.csv");
        csv << "Max forward coverage,Max backward coverage,"
            "In reference,Not in reference,Not in reference ratio\n";
        for(const auto&[key, value]: m) {
            csv << key.maxForwardCoverage << ",";
            csv << key.maxBackwardCoverage << ",";
            csv << value.inReferenceCount << ",";
            csv << value.notInReferenceCount << ",";
            csv << double(value.notInReferenceCount) / double(value.inReferenceCount + value.notInReferenceCount) << ",";
            csv << "\n";
        }
    }

    cout << "Assembler::analyzeAnchors ends." << endl;
}
