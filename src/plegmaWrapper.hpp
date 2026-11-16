#pragma once

// C++ wrapper to plegma (https://github.com/kokyriakidis/plegma),
// a partial order aligner used by LocalAssembly7::Method::Plegma.
// plegma is only available if shasta2 was built with
// cmake option BUILD_WITH_PLEGMA (see scripts/BuildPlegma.py).

#include "cstdint.hpp"
#include "memory.hpp"
#include "vector.hpp"

namespace shasta2 {

    class Base;
    class HomopolymerModel;

    // Where a sequence given to plegma is fixed:
    // - Both: on both anchors of the local assembly.
    // - Left: on the left anchor (anchor A) only.
    // - Right: on the right anchor (anchor B) only.
    // plegma aligns a sequence fixed on one side only
    // from that side and lets it end anywhere.
    enum class PlegmaAnchor {
        Both,
        Left,
        Right
    };

    // One oriented read given to plegma.
    class PlegmaRead {
    public:
        // The sequence of the oriented read. Oriented reads with the same
        // sequence can point to the same vector<Base>. It must stay valid
        // while plegma runs.
        const vector<Base>* sequence;

        PlegmaAnchor anchor;

        // The OrientedReadId strand (0 or 1). plegma uses it to look up
        // the homopolymer model, which has a separate error profile
        // for each strand.
        uint64_t strand;
    };

    // Returns true if this build includes plegma.
    bool plegmaIsAvailable();

    // Compute the consensus of the given oriented reads with plegma.
    // Each oriented read is entered once, in the order given.
    // If homopolymerModelPointer is not null, plegma uses that model
    // to choose the length of each homopolymer run in the consensus.
    // Returns false if plegma could not join the two anchors (this can only
    // happen if no oriented read is on both anchors), in which case
    // the consensus is left empty.
    // Throws if this build does not include plegma, or if plegma fails.
    bool plegma(
        const vector<PlegmaRead>&,
        const shared_ptr<const HomopolymerModel>& homopolymerModelPointer,
        vector<Base>& consensus);

}
