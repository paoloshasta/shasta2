#pragma once

/*******************************************************************************

Among all 4^k k-mers of length k, we choose a subset that we call "markers".
The markers are selected at the beginning of an assembly
and never changed, and selected in such a way that,
if (and only if) a k-mer is a marker, its reverse complement
is also a marker.

*******************************************************************************/

#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "Uint.hpp"

namespace shasta {

    class Marker;
    class Markers;

    class KmerChecker;
    class Reads;
}



// Markers in shared memory are stored using class Marker.
class shasta::Marker {
public:

    // The position of this marker in the oriented read.
    // This limits the length of a read to 2^24=16Mib bases.
    Uint24 position;
};



// The markers on all oriented reads. Indexed by OrientedReadId::getValue().
class shasta::Markers:
    public MappedMemoryOwner,
    public MemoryMapped::VectorOfVectors<Marker, uint64_t> {
public:

    Markers(
        const MappedMemoryOwner&,
        size_t k,
        const KmerChecker&,
        const Reads& reads,
        size_t threadCount);

    Markers(const MappedMemoryOwner&);
};
