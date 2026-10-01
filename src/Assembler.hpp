#pragma once

// Shasta.
#include "Options.hpp"
#include "HttpServer.hpp"
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedObject.hpp"
#include "MemoryMappedVector.hpp"
#include "MultithreadedObject.hpp"
#include "shastaTypes.hpp"

// Standard library.
#include "memory.hpp"
#include "string.hpp"
#include "utility.hpp"

namespace shasta2 {

    class Assembler;

    class AssemblyGraphPostprocessor;
    class AssemblyGraphPostprocessor;
    class AnchorGraph;
    class Anchors;
    class AssemblerInfo;
    class FastaLoader;
    class HomopolymerModel;
    class Journeys;
    class KmerChecker;
    class KmersOptions;
    class LongBaseSequences;
    class Markers;
    class MarkerKmers;
    class ReadGraph;
    class Reads;
    class ReadSummary;


    // Write an html form to select strand.
    void writeStrandSelection(
        ostream&,               // The html stream to write the form to.
        const string& name,     // The selection name.
        bool select0,           // Whether strand 0 is selected.
        bool select1);          // Whether strand 1 is selected.


    extern template class MultithreadedObject<Assembler>;
}



// Class used to store various pieces of assembler information in shared memory.
class shasta2::AssemblerInfo {
public:

    // The length of k-mers used to define markers.
    uint64_t k;

    // The marker density.
    double markerDensity;

    // The page size in use for this run.
    uint64_t largeDataPageSize;
};



class shasta2::Assembler :
    public MultithreadedObject<Assembler>,
    public MappedMemoryOwner,
    public HttpServer {
public:


    /***************************************************************************

    The constructors specify the file name prefix for binary data files.
    If this is a directory name, it must include the final "/".

    The constructor for a new assembler also specifies the page size for binary data files.
    Typically, for a large run binary data files will reside in a huge page
    file system backed by 2MB pages.
    The page sizes specified here must be equal to, or be an exact multiple of,
    the actual size of the pages backing the data.

    ***************************************************************************/

    // Construct a new Assembler.
    Assembler(
        const string& largeDataFileNamePrefix,
        size_t largeDataPageSize);

    // Construct an Assembler from binary data. This accesses the AssemblerInfo and the Reads.
    Assembler(const string& largeDataFileNamePrefix);



    // Various pieces of assembler information stored in shared memory.
    MemoryMapped::Object<AssemblerInfo> assemblerInfo;

    // This runs the entire assembly, under the following assumptions:
    // - The current directory is the run directory.
    // - The Data directory has already been created and set up, if necessary.
    // - The input file names are either absolute,
    //   or relative to the run directory, which is the current directory.
    void assemble(
        const Options& options,
        const vector<string>& inputFileNames,
        const string& externalAnchorsNameAbsolutePath,
        const string& externalAnchorGraphNameAbsolutePath);

    // The homopolymer model used by msaRepair, if one was specified with
    // --homopolymer-model. Null otherwise.
    shared_ptr<const HomopolymerModel> homopolymerModelPointer;
    void createHomopolymerModel(const string& homopolymerName);


    // Reads.
    shared_ptr<Reads> readsPointer;
    const Reads& reads() const {
        SHASTA2_ASSERT(readsPointer);
        return *readsPointer;
    }
    void computeReadIdsSortedByName();
    void addReads(
        const vector<string>& fileNames,
        uint64_t minReadLength,
        size_t threadCount);
    void addReads(
        const string& fileName,
        uint64_t minReadLength,
        size_t threadCount,
        FastaLoader&);
    void histogramReadLength(const string& fileName);

    void findPalindromicReads();
    void findPalindromicReadsMultithreaded(uint64_t threadCount);
    void findPalindromicReadsThreadFunction(uint64_t threadId);
    double analyzeStrandReversal(ReadId, bool debug) const;



    // Read summary information.
    MemoryMapped::Vector<ReadSummary> readSummaries;
    void createReadSummaries();
    void accessReadSummaries();
    void writeReadSummaries(bool partial) const;



    // The KmerChecker is used to find out if a given Kmer is a marker.
    shared_ptr<KmerChecker> kmerChecker;
    public:
    void createKmerChecker(uint64_t k, double markerDensity);



    // The markers on all oriented reads.
    shared_ptr<Markers> markersPointer;
    const Markers& markers() const
    {
        SHASTA2_ASSERT(markersPointer);
        return *markersPointer;
    }
    void checkMarkersAreOpen() const;
    void createMarkers(size_t threadCount);
    void accessMarkers();



    // The MarkerKmers keep track of the locations in the oriented reads
    // where each marker k-mer appears.
    shared_ptr<MarkerKmers> markerKmers;
    void createMarkerKmers(double maxMarkerErrorRate, uint64_t threadCount);
    void accessMarkerKmers();

    // Compute marker error rates for each read.
    // This computes the number of low frequency markers for each ReadId
    // and stores it in the lowFrequencyMarkerCount vector.
    void computeMarkerErrorRates();
    void computeMarkerErrorRatesMultithreaded(uint64_t threadCount);
    void computeMarkerErrorRatesThreadFunction(uint64_t threadId);
    vector<uint64_t> lowFrequencyMarkerCount;



    // Anchors.
    shared_ptr<Anchors> anchorsPointer;
    const Anchors& anchors() const
    {
        SHASTA2_ASSERT(anchorsPointer);
        return *anchorsPointer;
    }

    // Generate Anchors from MarkerKmers.
    void createAnchors(
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage,
        const vector<uint64_t>& maxAnchorRepeatLength,
        const vector<uint64_t>& minAnchorDistinctSubkmerCount,
        uint64_t threadCount);

    // Read Anchors from ExternalAnchors.
    void readExternalAnchors(const string& name);

    // Access existing Anchors.
    void accessAnchors();



    // Journeys.
    shared_ptr<Journeys> journeysPointer;
    const Journeys& journeys() const
    {
        SHASTA2_ASSERT(journeysPointer);
        return *journeysPointer;
    }
    void createJourneys(uint64_t threadCount);
    void accessJourneys();

    // Store anchor gaps information in ReadSummary for each read.
    void storeAnchorGaps();



    // AnchorGraph.
    shared_ptr<AnchorGraph> anchorGraphPointer;
    void createAnchorGraph(const Options&);
    void accessAnchorGraph(string name = "");
    void saveAnchorGraph();
    void anchorGraphTransitiveReduction(const Options&);



    // The complete AnchorGraph, which includes all possible edges
    // generated from the Journeys.
    // This is only used in the Python API and in the http server.
    // It is not used in the standard assembly process.
    shared_ptr<AnchorGraph> completeAnchorGraphPointer;
    void createCompleteAnchorGraph();
    void accessCompleteAnchorGraph();



    // AssemblyGraph.
    void createAssemblyGraph(const Options&, bool removeAnchorGraph);



    // Data and functions used for the http server.
    // This function puts the server into an endless loop
    // of processing requests.
    void writeHtmlBegin(ostream&) const;
    void writeHtmlEnd(ostream&) const;
    static void writeStyle(ostream& html);


    void writeNavigation(ostream&) const;
    void writeNavigation(
        ostream& html,
        const string& title,
        const vector<pair <string, string> >&) const;

    static void writePngToHtml(
        ostream& html,
        const string& pngFileName,
        const string useMap = ""
        );
    static void writeGnuPlotPngToHtml(
        ostream& html,
        int width,
        int height,
        const string& gnuplotCommands);

    void fillServerFunctionTable();
    void processRequest(
        const vector<string>& request,
        ostream&,
        const BrowserInformation&) override;
    void exploreSummary(const vector<string>&, ostream&);
    void exploreReadRaw(const vector<string>&, ostream&);
    void exploreLookupRead(const vector<string>&, ostream&);
    void exploreReadSequence(const vector<string>&, ostream&);
    void exploreReadMarkers(const vector<string>&, ostream&);
    void exploreMarkerKmer(const vector<string>&, ostream&);
    void exploreFindMarkerKmers(const vector<string>&, ostream&);
    static void addScaleSvgButtons(ostream&, uint64_t sizePixels);

    class HttpServerData {
    public:

        using ServerFunction = void (Assembler::*) (
            const vector<string>& request,
            ostream&);
        std::map<string, ServerFunction> functionTable;

        const Options* options = 0;
    };
    HttpServerData httpServerData;

    // Access all available assembly data, without thorwing an exception
    // on failures.
    void accessAllSoft();


    void exploreAnchor(const vector<string>&, ostream&);
    void exploreAnchorPair2(const vector<string>&, ostream&);
    void exploreJourney(const vector<string>&, ostream&);
    void exploreLocalReadAnchorGraph(const vector<string>&, ostream&);
    void exploreLocalAnchorGraph(const vector<string>&, ostream&);
    void exploreSegments(const vector<string>&, ostream&);
    void exploreSegmentSequence(const vector<string>&, ostream&);
    void exploreSegmentSteps(const vector<string>&, ostream&);
    void exploreSegmentStepSupport(const vector<string>&, ostream&);
    void exploreSegmentStep(const vector<string>&, ostream&);
    void exploreTangleMatrix(const vector<string>&, ostream&);
    void exploreSegmentPair(const vector<string>&, ostream&);
    void exploreSimilarSequences(const vector<string>&, ostream&);

    // Get the AssemblyGraph for a given assembly stage.
    AssemblyGraphPostprocessor& getAssemblyGraph(
        const string& assemblyStage,
        const Options&);
    std::map<string, shared_ptr<AssemblyGraphPostprocessor> > assemblyGraphTable;



    // Remove the memory mapped objects owned by an object of type T,
    // if the options allow it. Type T must implement close()
    // and remove().
    template<class T> void removeIfAllowed(const Options& options, T& t) const
    {
        if(options.memoryMode == "filesystem") {
            if(options.keepBinaryData) {
                t.close();
            } else {
                t.remove();
            }
        } else {
            t.close();
        }
    }

};

