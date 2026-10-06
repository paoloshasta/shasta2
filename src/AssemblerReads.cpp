// Shasta.
#include "Assembler.hpp"
#include "FastaLoader.hpp"
#include "memoryInformation.hpp"
#include "performanceLog.hpp"
#include "ReadLoader.hpp"
#include "timestamp.hpp"
using namespace shasta2;

// Standard libraries.
#include "algorithm.hpp"
#include "chrono.hpp"
#include <cmath>
#include "iterator.hpp"



void Assembler::addReads(
    const vector<string>& fileNames,
    uint64_t minReadLength,
    size_t threadCount)
{
    writeMemoryStatistics("Assembler::addReads begins");

    const auto t0 = steady_clock::now();

    FastaLoader fastaLoader(minReadLength, threadCount, *readsPointer);
    for(const string& fileName: fileNames) {
        addReads(
            fileName,
            minReadLength,
            threadCount,
            fastaLoader);
    }

    // Free up unused allocated memory.
    readsPointer->unreserve();

    reads().checkSanity();
    readsPointer->computeReadLengthHistogram();

    if(reads().readCount() == 0) {
        throw runtime_error("There are no input reads.");
    }
    computeReadIdsSortedByName();
    histogramReadLength("ReadLengthHistogram.csv");

    // Check that no reads are too long.
    // This limitation is because Marker::position is Uint24.
    const uint64_t maxReadLength = (2 << 24) - 1;
    for(ReadId readId=0; readId<reads().readCount(); readId++) {
        if(reads().getReadSequenceLength(readId) > maxReadLength ) {
            std::ostringstream message;
            message << "Read " << to_string(readId) << " " << reads().getReadName(readId) <<
                " is too long. Maximum length allowed is " << maxReadLength << endl;
            throw runtime_error(message.str());
        }
    }

    const auto t1 = steady_clock::now();
    performanceLog << "Read loading took " << seconds(t1-t0) << "s." << endl;

    writeMemoryStatistics("Assembler::addReads ends");
}



// Add reads.
// The reads are added to those already previously present.
void Assembler::addReads(
    const string& fileName,
    uint64_t minReadLength,
    const size_t threadCount,
    FastaLoader& fastaLoader)
{
    const string extension = std::filesystem::path(fileName).extension();

    if(extension == ".fasta") {
        fastaLoader.read(fileName);
    } else {
        ReadLoader readLoader(
            fileName,
            minReadLength,
            threadCount,
            largeDataFileNamePrefix,
            largeDataPageSize,
            *readsPointer);
    }
}


// Create a histogram of read lengths.
// All lengths here are raw sequence lengths
// (length of the original read), not lengths
// in run-length representation.
void Assembler::histogramReadLength(const string& fileName)
{
    readsPointer->computeReadLengthHistogram();
    readsPointer->writeReadLengthHistogram(fileName);


    cout << "Number of reads: " << readsPointer->readCount() << endl;
    cout << "Number of bases: " << readsPointer->getTotalBaseCount() << endl;
    cout << "Average read length: " <<
        uint64_t(std::round(double(readsPointer->getTotalBaseCount()) / double(readsPointer->readCount()))) << endl;
    cout << "Read N50: " << readsPointer->getN50() << endl;

}



void Assembler::computeReadIdsSortedByName()
{
    readsPointer->computeReadIdsSortedByName();
}

