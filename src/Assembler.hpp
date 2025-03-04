#pragma once

// Shasta.
#include "AssemblerOptions.hpp"
#include "HttpServer.hpp"
#include "invalid.hpp"
#include "Kmer.hpp"
#include "MappedMemoryOwner.hpp"
#include "Marker.hpp"
#include "MemoryMappedObject.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "ReadId.hpp"
#include "shastaTypes.hpp"

// Standard library.
#include "memory.hpp"
#include "string.hpp"
#include "utility.hpp"

namespace shasta {

    class Assembler;
    class AssemblerInfo;
    class AssemblerOptions;
    class KmerChecker;
    class KmersOptions;
    class LongBaseSequences;
    class MarkerKmers;
    class Mode3Assembler;
    class Reads;

    namespace mode3 {
        class Anchors;
    }

    // Write an html form to select strand.
    void writeStrandSelection(
        ostream&,               // The html stream to write the form to.
        const string& name,     // The selection name.
        bool select0,           // Whether strand 0 is selected.
        bool select1);          // Whether strand 1 is selected.


    extern template class MultithreadedObject<Assembler>;
}



// Class used to store various pieces of assembler information in shared memory.
class shasta::AssemblerInfo {
public:

    // The length of k-mers used to define markers.
    size_t k;

    // The page size in use for this run.
    size_t largeDataPageSize;
};



class shasta::Assembler :
    public MultithreadedObject<Assembler>,
    public MappedMemoryOwner,
    public HttpServer {
public:


    /***************************************************************************

    The constructors specify the file name prefix for binary data files.
    If this is a directory name, it must include the final "/".

    The constructor also specifies the page size for binary data files.
    Typically, for a large run binary data files will reside in a huge page
    file system backed by 2MB pages.
    The page sizes specified here must be equal to, or be an exact multiple of,
    the actual size of the pages backing the data.

    ***************************************************************************/

    // Constructor.
    Assembler(
        const string& largeDataFileNamePrefix,
        bool createNew,
        size_t largeDataPageSize);



    // Various pieces of assembler information stored in shared memory.
    MemoryMapped::Object<AssemblerInfo> assemblerInfo;



    // Reads.
    shared_ptr<Reads> reads;
    const Reads& getReads() const {
        SHASTA_ASSERT(reads);
        return *reads;
    }
    void computeReadIdsSortedByName();
    void addReads(
        const string& fileName,
        uint64_t minReadLength,
        size_t threadCount);
    void histogramReadLength(const string& fileName);



    // The KmerChecker is used to find out if a given KmerId is a marker.
    shared_ptr<KmerChecker> kmerChecker;
    public:
    void createKmerChecker(
        const KmersOptions& kmersOptions,
        uint64_t threadCount);
    void accessKmerChecker();



    // The markers on all oriented reads. Indexed by OrientedReadId::getValue().
    MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t> markers;
    void checkMarkersAreOpen() const;

    // Get markers sorted by KmerId for a given OrientedReadId.
    void getMarkersSortedByKmerId(
        OrientedReadId,
        vector<MarkerWithOrdinal>&) const;

    void findMarkers(size_t threadCount);
    void accessMarkers();
    void writeMarkers(ReadId, Strand, const string& fileName);

    // Given a marker by its OrientedReadId and ordinal,
    // return the corresponding global marker id.
    MarkerId getMarkerId(OrientedReadId, uint32_t ordinal) const;

    MarkerId getReverseComplementMarkerId(OrientedReadId, uint32_t ordinal) const;
    MarkerId getMarkerId(const MarkerDescriptor& m) const
    {
        return getMarkerId(m.first, m.second);
    }
    MarkerId getReverseComplementMarkerId(const MarkerDescriptor& m) const
    {
        return getReverseComplementMarkerId(m.first, m.second);
    }

    // Inverse of the above: given a global marker id,
    // return its OrientedReadId and ordinal.
    // This requires a binary search in the markers toc.
    // This could be avoided, at the cost of storing
    // an additional 4 bytes per marker.
    pair<OrientedReadId, uint32_t> findMarkerId(MarkerId) const;

    // Low level functions to get marker Kmers/KmerIds of an oriented read.
    // They are obtained from the reads and not from CompressedMarker::kmerId,
    // which will soon go away.

    // Get the marker Kmer for an oriented read and ordinal.
    Kmer getOrientedReadMarkerKmer(OrientedReadId, uint32_t ordinal) const;
    Kmer getOrientedReadMarkerKmerStrand0(ReadId, uint32_t ordinal) const;
    Kmer getOrientedReadMarkerKmerStrand1(ReadId, uint32_t ordinal) const;

    // Get the marker KmerId for an oriented read and ordinal.
    KmerId getOrientedReadMarkerKmerId(OrientedReadId, uint32_t ordinal) const;

    // Get all marker Kmers for an oriented read.
    void getOrientedReadMarkerKmers(OrientedReadId, const span<Kmer>&) const;
    void getOrientedReadMarkerKmersStrand0(ReadId, const span<Kmer>&) const;
    void getOrientedReadMarkerKmersStrand1(ReadId, const span<Kmer>&) const;

    // Get all marker KmerIds for an oriented read.
    void getOrientedReadMarkerKmerIds(OrientedReadId, const span<KmerId>&) const;
    void getOrientedReadMarkerKmerIdsStrand0(ReadId, const span<KmerId>&) const;
    void getOrientedReadMarkerKmerIdsStrand1(ReadId, const span<KmerId>&) const;

    // Get all MarkerWithOrdinals for an oriented read (includes position, KmerId, and ordinal).
    void getOrientedReadMarkers(OrientedReadId, const span<MarkerWithOrdinal>&) const;
    void getOrientedReadMarkersStrand0(ReadId, const span<MarkerWithOrdinal>&) const;
    void getOrientedReadMarkersStrand1(ReadId, const span<MarkerWithOrdinal>&) const;

    // Get all marker Kmers/KmerIds for a read in both orientations.
    void getReadMarkerKmers(
        ReadId,
        const span<Kmer>& Kmers0,
        const span<Kmer>& Kmers1) const;
    void getReadMarkerKmerIds(
        ReadId,
        const span<KmerId>& kmerIds0,
        const span<KmerId>& kmerIds1) const;

    // Get the Kmer/KmerId for an oriented read at a given marker ordinal.
    Kmer getOrientedReadMarkerKmer(OrientedReadId, uint64_t ordinal) const;
    KmerId getOrientedReadMarkerKmerId(OrientedReadId, uint64_t ordinal) const;

    // The MarkerKmers keep track of the locations in the oriented reads
    // where each marker k-mer appears.
    shared_ptr<MarkerKmers> markerKmers;
    void createMarkerKmers(uint64_t threadCount);
    void accessMarkerKmers();

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
    void exploreMarkerKmers(const vector<string>&, ostream&);
    static void addScaleSvgButtons(ostream&, uint64_t sizePixels);

    class HttpServerData {
    public:

        using ServerFunction = void (Assembler::*) (
            const vector<string>& request,
            ostream&);
        std::map<string, ServerFunction> functionTable;

        const AssemblerOptions* assemblerOptions = 0;
    };
    HttpServerData httpServerData;




    void writeMakeAllTablesCopyable(ostream&) const;


    // Access all available assembly data, without thorwing an exception
    // on failures.
    void accessAllSoft();


    // Mode 3 assembly.
    shared_ptr<Mode3Assembler> mode3Assembler;
    void accessMode3Assembler();


    // Http server functions related to Mode 3 assembly.
    void exploreAnchor(const vector<string>&, ostream&);
    void exploreAnchorPair(const vector<string>&, ostream&);
    void exploreJourney(const vector<string>&, ostream&);
    void exploreReadFollowing(const vector<string>&, ostream&);
    void exploreLocalAssembly(const vector<string>&, ostream&);
    void exploreLocalAnchorGraph(const vector<string>&, ostream&);

};

