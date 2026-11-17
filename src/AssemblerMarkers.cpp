// Shasta.
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "extractKmer128.hpp"
#include "htmlLog.hpp"
#include "Markers.hpp"
#include "memoryInformation.hpp"
#include "MarkerKmers.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "ReadSummary.hpp"
#include "timestamp.hpp"
using namespace shasta2;

// Standard library.
#include "fstream.hpp"


void Assembler::createMarkers(size_t threadCount)
{
    readsPointer->checkReadsAreOpen();
    SHASTA2_ASSERT(kmerChecker);

    markersPointer = make_shared<Markers>(
        *this,
        assemblerInfo->k,
        readsPointer,
        kmerChecker,
        threadCount);
}



void Assembler::accessMarkers()
{
    markersPointer = make_shared<Markers>(*this, assemblerInfo->k, readsPointer);
}



void Assembler::checkMarkersAreOpen() const
{
    if(not(markersPointer and !markersPointer->isOpen())) {
        throw runtime_error("Markers are not accessible.");
    }
}



void Assembler::createMarkerKmers(double maxMarkerErrorRate, uint64_t threadCount)
{
    writeMemoryStatistics("Assembler::createMarkerKmers begins");

    const MappedMemoryOwner& mappedMemoryOwner = *this;
    const uint64_t readCount = reads().readCount();

    // First, create marker k-mers using all reads.
    markerKmers = make_shared<MarkerKmers>(
        assemblerInfo->k,
        mappedMemoryOwner,
        reads(),
        readSummaries,
        markers(),
        threadCount);

    // Compute marker error rates using all reads and store them in the ReadSummaries.
    computeMarkerErrorRatesMultithreaded(threadCount);
    for(ReadId readId=0; readId<readCount; readId++) {
        const auto readMarkers = markers()[OrientedReadId(readId, 0).getValue()];
        const uint64_t readMarkerCount = readMarkers.size();
        const double markerErrorRate = double(lowFrequencyMarkerCount[readId]) / double(readMarkerCount);
        ReadSummary& readSummary = readSummaries[readId];
        readSummary.initialMarkerErrorRate = markerErrorRate;
        readSummary.hasHighErrorRate = (markerErrorRate > maxMarkerErrorRate);
    }


    // Count the reads and bases that will be used for assembly.
    uint64_t keepCount = 0;
    uint64_t baseKeepCount = 0;
    for(ReadId readId=0; readId<readCount; readId++) {
        if(readSummaries[readId].isInUse()) {
            ++keepCount;
            baseKeepCount += reads().getRead(readId).baseCount;
        }
    }
    cout << "Kept " << keepCount << " reads with low marker error rate out of " <<
        readCount << " total (" << baseKeepCount <<
        " bases out of " << reads().getTotalBaseCount() << " total)." << endl;

    // Do it again, this time using only the reads with low marker error rate.
    markerKmers->remove();
    markerKmers = make_shared<MarkerKmers>(
        assemblerInfo->k,
        mappedMemoryOwner,
        reads(),
        readSummaries,
        markers(),
        threadCount);

    // Compute marker error rates.
    computeMarkerErrorRatesMultithreaded(threadCount);
    for(ReadId readId=0; readId<readCount; readId++) {
        const auto readMarkers = markers()[OrientedReadId(readId, 0).getValue()];
        const uint64_t readMarkerCount = readMarkers.size();
        const double markerErrorRate = double(lowFrequencyMarkerCount[readId]) / double(readMarkerCount);
        ReadSummary& readSummary = readSummaries[readId];
        readSummary.markerErrorRate = markerErrorRate;
    }

    writeMemoryStatistics("Assembler::createMarkerKmers ends");
}



void Assembler::writeMarkerKmersInformationToHtmlLog(uint64_t minAnchorCoverage, uint64_t maxAnchorCoverage)
{
    const string shade = 
        "set object 1 rectangle from " + to_string(minAnchorCoverage) + 
        ", graph 0 to " + to_string(maxAnchorCoverage) + 
        ", graph 1 behind fillcolor rgb \"#ffcccc\" fillstyle solid 0.5 noborder";

    const string gnuplotCommands = shade +  R"(
        set datafile separator ","
        set key noautotitle

        set xtics out nomirror
        set ytics out nomirror
        set grid xtics mxtics ytics linestyle 1 linewidth 1 linecolor rgb '#e0e0e0'

        set xlabel "Coverage"
        set ylabel "Frequency"
        set title "Marker k-mer frequency histogram"

        set logscale xy

        plot "MarkerKmers-Histogram.csv" using 1:2 with points pointtype 7 pointsize 0.5 linecolor rgb "#4c6691"
    )";

    const uint64_t width = 800;
    const uint64_t height = 600;

    try {   
        writeGnuPlotPngToHtml(htmlLog, width, height,gnuplotCommands);
    } catch (std::exception& e) {
        htmlLog << "<br>The marker k-mer histogram could not be shown, "
            "possibly because the gnuplot and/or base64 commands are not available.";
    }
    htmlLog << 
        "<br>Anchors are created from a subset of the marker k-mers in the shaded region."
        "<br>See <a href='MarkerKmers-Histogram.csv'>MarkerKmers-Histogram.csv</a>." << flush;
}



void Assembler::accessMarkerKmers()
{
    const MappedMemoryOwner& mappedMemoryOwner = *this;

    markerKmers = make_shared<MarkerKmers>(
        assemblerInfo->k,
        mappedMemoryOwner,
        reads(),
        markers());
}



// This computes the number of low frequency markers for each ReadId
// and stores it in the lowFrequencyMarkerCount vector.
void Assembler::computeMarkerErrorRates()
{
    SHASTA2_ASSERT(markerKmers);
    const uint64_t k = assemblerInfo->k;
    const uint64_t readCount = reads().readCount();

    lowFrequencyMarkerCount.clear();
    lowFrequencyMarkerCount.resize(readCount, 0);

    // Loop over all reads.
    for(ReadId readId=0; readId<readCount; readId++) {
        const auto readSequence = reads().getRead(readId);
        const auto readMarkers = markers()[OrientedReadId(readId, 0).getValue()];

        // Count frequency 1 marker kmers in this read.
        for(const Marker& marker: readMarkers) {
            const uint32_t position = marker.position;
            Kmer kmer;
            extractKmer128(readSequence, position, k, kmer);
            const uint64_t frequency = markerKmers->getFrequency(kmer);
            // Also count frequency 0.
            // Possible, for reads for which useRead is false;
            if(frequency <= 1) {
                ++lowFrequencyMarkerCount[readId];
            }
        }
    }
}


// This computes the number of low frequency markers for each ReadId
// and stores it in the lowFrequencyMarkerCount vector.
void Assembler::computeMarkerErrorRatesMultithreaded(uint64_t threadCount)
{
    const uint64_t readCount = reads().readCount();
    lowFrequencyMarkerCount.clear();
    lowFrequencyMarkerCount.resize(readCount, 0);

    const uint64_t batchCount = 10;
    setupLoadBalancing(readCount, batchCount);
    runThreads(&Assembler::computeMarkerErrorRatesThreadFunction, threadCount);

}



void Assembler::computeMarkerErrorRatesThreadFunction([[maybe_unused]] uint64_t threadId)
{
    const uint64_t k = assemblerInfo->k;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all ReadIds in this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {

            const auto readSequence = reads().getRead(readId);
            const auto readMarkers = markers()[OrientedReadId(readId, 0).getValue()];

            // Count frequency 1 marker kmers in this read.
            uint64_t lowFrequencyCount = 0;
            for(const Marker& marker: readMarkers) {
                const uint32_t position = marker.position;
                Kmer kmer;
                extractKmer128(readSequence, position, k, kmer);
                const uint64_t frequency = markerKmers->getFrequency(kmer);
                // Also count frequency 0.
                // Possible, for reads for which useRead is false;
                if(frequency <= 1) {
                    ++lowFrequencyCount;
                }
            }
            lowFrequencyMarkerCount[readId] = lowFrequencyCount;

        }
    }

}



void Assembler::findPalindromicReads()
{
    // EXPOSE WHEN CODE STABILIZES.
    const double palindromicRateThreshold = 0.2;


    const ReadId readCount = reads().readCount();

    // Loop over all reads.
    uint64_t palindromicReadCount = 0;
    for(ReadId readId=0; readId<readCount; readId++) {
        ReadSummary& readSummary = readSummaries[readId];

        readSummary.palindromicRate = analyzeStrandReversal(readId, false);

        if(readSummary.palindromicRate > palindromicRateThreshold) {
            readSummary.isPalindromic = true;
            ++palindromicReadCount;
        }
    }

    cout << "Flagged " << palindromicReadCount << " reads as palindromic." << endl;
}



void Assembler::findPalindromicReadsMultithreaded(uint64_t threadCount)
{
    const ReadId readCount = reads().readCount();
    const uint64_t batchSize = 10;
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::findPalindromicReadsThreadFunction, threadCount);

    uint64_t palindromicReadCount = 0;
    for(ReadId readId=0; readId<readCount; readId++) {
        if(readSummaries[readId].isPalindromic) {
            ++palindromicReadCount;
        }
    }
    cout << "Flagged " << palindromicReadCount << " reads as palindromic." << endl;
}



void Assembler::findPalindromicReadsThreadFunction([[maybe_unused]] uint64_t threadId)
{
    // EXPOSE WHEN CODE STABILIZES.
    const double palindromicRateThreshold = 0.2;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all ReadIds in this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {
            ReadSummary& readSummary = readSummaries[readId];

            readSummary.palindromicRate = analyzeStrandReversal(readId, false);

            if(readSummary.palindromicRate > palindromicRateThreshold) {
                readSummary.isPalindromic = true;
            }

        }
    }

}





double Assembler::analyzeStrandReversal(ReadId readId, bool debug) const
{
    // EXPOSE WHEN CODE STABILIZES.
    const double driftRate = 0.01;

    const uint64_t k = assemblerInfo->k;

    // Access the sequence of this read (without reverse complementing).
    const LongBaseSequenceView sequence0 = reads().getRead(readId);
    const uint32_t length = uint32_t(sequence0.baseCount);
    const uint32_t maxDrift = uint32_t(std::round(driftRate * double(length)));

    if(debug) {
        cout << "Length " << length << endl;
        cout << "Max drift " << maxDrift << endl;
    }

    // Access the markers of this read (without reverse complementing.
    const auto markers0 = markers()[OrientedReadId(readId, 0).getValue()];


    // Gather all marker Kmers of this read and its reverse complement,
    // and the positions in which each Kmer appears.
    std::map<Kmer, vector<uint32_t> > kmerTable0;
    std::map<Kmer, vector<uint32_t> > kmerTable1;

    for(const Marker& marker0: markers0) {
        const uint32_t position0 = marker0.position;
        Kmer kmer0;
        extractKmer128(sequence0, position0, k, kmer0);
        kmerTable0[kmer0].push_back(position0);

        const uint32_t position1 = length - 1 - position0;
        const Kmer kmer1 = kmer0.reverseComplement(k);
        kmerTable1[kmer1].push_back(position1);
    }

    for(auto& p: kmerTable1) {
        vector<uint32_t>& positions1 = p.second;
        std::ranges::reverse(positions1);
    }


    if(debug) {
        // Loop over common marker k-mers.
        ofstream csv("AnalyzeStrandReversal.csv");
        for(const auto& p0: kmerTable0) {
            const Kmer& kmer = p0.first;
            auto it1 = kmerTable1.find(kmer);
            if(it1 != kmerTable1.end()) {
                const auto& p1 = *it1;
                const vector<uint32_t>& positions0 = p0.second;
                const vector<uint32_t>& positions1 = p1.second;
                for(const uint32_t position0: positions0) {
                    for(const uint32_t position1: positions1) {
                        csv <<
                            position0 << "," <<
                            position1 << "," <<
                            int32_t(position1) - int32_t(position0) << "\n";
                    }
                }
            }
        }
    }


    // To decide if this read is palindromic, loop over all its marker kmers.
    // Find how many times we find the same k-mer at a nearby position
    // in the reverse complemented read.
    uint64_t totalCount = 0;
    uint64_t successCount = 0;
    for(const auto& p0: kmerTable0) {
        const Kmer& kmer = p0.first;
        const vector<uint32_t>& positions0 = p0.second;
        totalCount += positions0.size();

        const auto it1 = kmerTable1.find(kmer);
        if(it1 != kmerTable1.end()) {
            const vector<uint32_t>& positions1 = it1->second;
            SHASTA2_ASSERT(not positions1.empty());

            // This can be made more efficient.
            for(uint64_t position0: positions0) {
                for(uint64_t position1: positions1) {
                    if(abs(int32_t(position0) - int32_t(position1)) < maxDrift) {
                        ++successCount;
                        break;
                    }
                }
            }
        }
    }

    const double successRate = double(successCount) / double(totalCount);

    if(debug) {
        cout << "Total " << totalCount << ", success " << successCount <<
            ", success rate " << successRate << endl;
    }

    return successRate;

}
