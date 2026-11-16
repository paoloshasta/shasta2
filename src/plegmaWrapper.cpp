// Shasta.
#include "plegmaWrapper.hpp"
#include "Base.hpp"
#include "HomopolymerModel.hpp"
using namespace shasta2;

// Standard library.
#include <map>
#include <mutex>
#include "stdexcept.hpp"
#include "string.hpp"

// plegma.
#ifdef SHASTA2_PLEGMA
#include "plegma.hpp"
#endif



bool shasta2::plegmaIsAvailable()
{
#ifdef SHASTA2_PLEGMA
    return true;
#else
    return false;
#endif
}



#ifdef SHASTA2_PLEGMA

// plegma reads the homopolymer model from the same csv text as shasta2.
// This makes the plegma model once for each HomopolymerModel and shares
// it between threads (plegma only reads it). The weak_ptr detects a
// HomopolymerModel that was destroyed and another one later created
// at the same address.
static const plegma::Model* getPlegmaModel(
    const shared_ptr<const HomopolymerModel>& homopolymerModelPointer)
{
    if(not homopolymerModelPointer) {
        return nullptr;
    }

    class Entry {
    public:
        std::weak_ptr<const HomopolymerModel> homopolymerModel;
        std::unique_ptr<plegma::Model> plegmaModel;
    };
    static std::mutex mutex;
    static std::map<const HomopolymerModel*, Entry> entries;

    std::lock_guard<std::mutex> lock(mutex);
    Entry& entry = entries[homopolymerModelPointer.get()];
    if(not entry.plegmaModel or entry.homopolymerModel.expired()) {
        entry.homopolymerModel = homopolymerModelPointer;
        entry.plegmaModel = std::make_unique<plegma::Model>(
            plegma::Model::parse(homopolymerModelPointer->csv()));
    }
    return entry.plegmaModel.get();
}



bool shasta2::plegma(
    const vector<PlegmaRead>& reads,
    const shared_ptr<const HomopolymerModel>& homopolymerModelPointer,
    vector<Base>& consensus)
{
    consensus.clear();

    // plegma takes the sequences as characters.
    // Convert each distinct sequence once.
    std::map<const vector<Base>*, string> texts;
    for(const PlegmaRead& read: reads) {
        string& text = texts[read.sequence];
        if(text.empty()) {
            text.reserve(read.sequence->size());
            for(const Base b: *read.sequence) {
                text.push_back(b.character());
            }
        }
    }

    vector<plegma::Read> plegmaReads;
    plegmaReads.reserve(reads.size());
    for(const PlegmaRead& read: reads) {
        plegma::Read& plegmaRead = plegmaReads.emplace_back();
        plegmaRead.sequence = texts.at(read.sequence);
        switch(read.anchor) {
        case PlegmaAnchor::Both:
            plegmaRead.anchor = plegma::Anchor::Both;
            break;
        case PlegmaAnchor::Left:
            plegmaRead.anchor = plegma::Anchor::Left;
            break;
        case PlegmaAnchor::Right:
            plegmaRead.anchor = plegma::Anchor::Right;
            break;
        }
        plegmaRead.strand = (read.strand == 0) ? plegma::Strand::Forward : plegma::Strand::Reverse;
    }

    // Unpinned: all the oriented reads of a local assembly cover the same region,
    // so each is aligned to the whole graph. Pinning each read to the consensus
    // so far is for long, divergent sequences, and here it only costs time
    // (2.4 times slower on a chr21 assembly, with the same consensus).
    const vector<string> parts = plegma::consensus(
        plegmaReads, getPlegmaModel(homopolymerModelPointer), plegma::Options().pinned(false));

    // plegma returns two parts, the left side and the right side,
    // if no read is on both anchors and it could not join them.
    if(parts.size() != 1) {
        return false;
    }

    consensus.reserve(parts.front().size());
    for(const char c: parts.front()) {
        consensus.push_back(Base::fromCharacter(c));
    }
    return true;
}



#else



bool shasta2::plegma(
    const vector<PlegmaRead>&,
    const shared_ptr<const HomopolymerModel>&,
    vector<Base>&)
{
    throw runtime_error("This build of shasta2 does not include plegma. "
        "To use the Plegma local assembly method, run scripts/BuildPlegma.py, "
        "then build shasta2 with cmake option -DBUILD_WITH_PLEGMA=ON.");
}

#endif
