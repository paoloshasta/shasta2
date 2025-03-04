#pragma once

// Shasta.
#include "mode3-Anchor.hpp"
#include "AssemblerOptions.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"
#include "ReadId.hpp"
#include "shastaTypes.hpp"

// Standard library.
#include "memory.hpp"
#include "span.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    class Mode3Assembler;

    class MarkerInterval;
    class Mode3AssemblyOptions;

    namespace mode3 {
        class AssemblyGraph;
        class AssemblyGraphPostprocessor;
    }

    namespace MemoryMapped {
        template<class T, class Int> class VectorOfVectors;
    }
}


class shasta::Mode3Assembler :
    public MultithreadedObject<Mode3Assembler>,
    public MappedMemoryOwner {
public:

    // This constructor runs the assembly.
    Mode3Assembler(
        const MappedMemoryOwner&,
        uint64_t k,
        const Reads&,
        const MemoryMapped::VectorOfVectors<Marker, uint64_t>& markers,
        shared_ptr<mode3::Anchors>,
        uint64_t threadCount,
        const Mode3AssemblyOptions&,
        bool debug);

    // This constructor just accesses binary data.
    Mode3Assembler(
        const MappedMemoryOwner&,
        uint64_t k,
        const Reads&,
        const MemoryMapped::VectorOfVectors<Marker, uint64_t>& markers,
        shared_ptr<mode3::Anchors>,
        const Mode3AssemblyOptions&);

private:
    uint64_t k;
    const Reads& reads;
    const MemoryMapped::VectorOfVectors<Marker, uint64_t>& markers;
    bool debug;

    // The main input to Mode3Assembler is a set of anchors.
    // Each anchor consists of a span of MarkerIntervals, with the following requirements:
    // - All MarkerIntervals correspond to exactly the same sequence in the corresponding oriented reads, and:
    //      * Those portions of the oriented reads are believed to be aligned.
    //      * They apear in a low number of copies in the genome being sequenced.
    // - There are no duplicate oriented reads in an anchor.
    // - The anchor coverage (number of oriented reads) is isn [minPrimaryCoverage, maxPrimaryCoverage].
    // For now the anchors are simply a reference to assembler.markerGraph.edgeMarkerIntervals,
    // butit might be possible to construct the anchors by other means.
    shared_ptr<mode3::Anchors> anchorsPointer;
public:
    mode3::Anchors& anchors()
    {
        return *anchorsPointer;
    }
    const mode3::Anchors& anchors() const
    {
        return *anchorsPointer;
    }

    const Mode3AssemblyOptions& options;



    // The oriented reads present in each anchor
    // define a bipartite graph. We want to compute connected components
    // of this bipartite graph and process them one at a time.
    // These are also connected components of the global anchor graph
    // (with one vertex for each anchor, and edges created by following the reads).

    // The OrientedReadIds in each component. Indexed by componentId.
    MemoryMapped::VectorOfVectors<OrientedReadId, uint64_t> componentOrientedReadIds;

    // The anchorIds in each component. Indexed by componentId.
    MemoryMapped::VectorOfVectors<mode3::AnchorId, uint64_t> componentAnchorIds;

public:

    // Http server functions for Mode 3 assembly.
    void exploreAnchor(const vector<string>&, ostream&);
    void exploreAnchorPair(const vector<string>&, ostream&);
    void exploreJourney(const vector<string>&, ostream&);
    void exploreReadFollowing(const vector<string>&, ostream&);
    void exploreLocalAssembly(
        const vector<string>&, ostream&);
    void exploreLocalAnchorGraph(
        const vector<string>&, ostream&);

};
