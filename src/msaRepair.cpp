// Shasta2.
#include "msaRepair.hpp"
#include "HomopolymerModel.hpp"
#include "invalid.hpp"
#include "SHASTA2_ASSERT.hpp"
using namespace shasta2;

// Theseus. Used directly, not through theseusWrapper: see msaRepairAlignExtended
// below for why.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#include "theseus/theseus_msa_aligner.h"
#pragma GCC diagnostic pop

// Standard library.
#include "algorithm.hpp"
#include "iostream.hpp"
#include <boost/tokenizer.hpp>
#include <sstream>



namespace shasta2 {

    // A row of an alignment is either bare AlignedBase symbols, or
    // AlignedExtendedBase symbols paired with a run length once it has been
    // through the extended alphabet - an AlignedExtendedSequence never
    // appears bare. These let the handful of functions that only care where
    // the gaps are work on both without being written twice.
    inline bool msaRepairIsGap(AlignedBase b) { return b.isGap(); }
    template<class Symbol> inline bool msaRepairIsGap(const pair<Symbol, uint64_t>& p)
    {
        return p.first.isGap();
    }

    // The columns of a row that hold its first and last base, as a half open
    // interval. An all gap row gives an empty interval.
    //
    // Which end of a row is padding and which is a deletion is decided from
    // this, in more than one place, so it has to mean the same thing in all of
    // them.
    template<class Row> pair<uint64_t, uint64_t> msaRepairBaseSpan(const Row& row)
    {
        uint64_t first = row.size();
        uint64_t last = 0;
        for(uint64_t j=0; j<row.size(); j++) {
            if(not msaRepairIsGap(row[j])) {
                if(first == row.size()) {
                    first = j;
                }
                last = j;
            }
        }
        if(first == row.size()) {
            return make_pair(0UL, 0UL);
        }
        return make_pair(first, last + 1);
    }

    // Count the columns of a range that hold more than one distinct base, gaps
    // ignored.
    //
    // This is the measure of how badly a stretch of alignment is misplaced: a
    // column holding both an A and a G means one of them is in the wrong place,
    // because a read that has neither contributes a gap, not a base. The repair
    // guards itself with this and the tests assert on it, so both must be the
    // same code or the guard and the test can disagree.
    uint64_t msaRepairImpureColumnCount(
        const vector< vector<AlignedBase> >& rows,
        uint64_t begin,
        uint64_t end)
    {
        uint64_t impureColumnCount = 0;
        for(uint64_t j=begin; j<end; j++) {
            array<bool, 4> seen = {false, false, false, false};
            for(const vector<AlignedBase>& row: rows) {
                const AlignedBase b = row[j];
                if(not b.isGap()) {
                    seen[b.value] = true;
                }
            }
            if(std::ranges::count(seen, true) > 1) {
                ++impureColumnCount;
            }
        }
        return impureColumnCount;
    }

    uint64_t msaRepairImpureColumnCount(const vector< vector<AlignedBase> >& rows)
    {
        return rows.empty() ? 0 :
                msaRepairImpureColumnCount(rows, 0, rows.front().size());
    }

    // An alignment written as one string per row, '-' for a gap. This is how the
    // tests and the stored example alignments spell an alignment out.
    vector< vector<AlignedBase> > msaRepairAlignmentFromStrings(const vector<string>& rows)
    {
        vector< vector<AlignedBase> > alignment;
        alignment.reserve(rows.size());
        for(const string& row: rows) {
            alignment.push_back(vectorOfAlignedBasesFromString(row));
        }
        return alignment;
    }

    uint64_t msaRepairImpureColumnCount(const vector<string>& rows)
    {
        return msaRepairImpureColumnCount(msaRepairAlignmentFromStrings(rows));
    }
}



// The look up tables used by the fromCharacter functions.
// Note these differ from BaseInitializer and AlignedBaseInitializer, which map
// lower case characters to the same values as upper case. Here lower case means
// a poly symbol, so 'a' and 'A' are deliberately different.
namespace shasta2 {

    // The eight symbols, shared by both tables. The gap belongs to only one of
    // them and is added by its initializer.
    static void fillExtendedBaseTable(std::array<uint8_t, 256>& table)
    {
        fill(table.begin(), table.end(), 255);
        table['A'] = 0;
        table['C'] = 1;
        table['G'] = 2;
        table['T'] = 3;
        table['a'] = 4;
        table['c'] = 5;
        table['g'] = 6;
        table['t'] = 7;
    }
}

ExtendedBaseInitializer ExtendedBaseInitializer::singleton;
std::array<uint8_t, 256> ExtendedBaseInitializer::table;
ExtendedBaseInitializer::ExtendedBaseInitializer()
{
    fillExtendedBaseTable(table);
}

AlignedExtendedBaseInitializer AlignedExtendedBaseInitializer::singleton;
std::array<uint8_t, 256> AlignedExtendedBaseInitializer::table;
AlignedExtendedBaseInitializer::AlignedExtendedBaseInitializer()
{
    fillExtendedBaseTable(table);
    table['-'] = AlignedExtendedBase::gapValue;
}



// See msaRepair.hpp for comments.
void shasta2::encodeExtended(
    const vector<Base>& sequence,
    uint64_t threshold,
    ExtendedSequence& encoded)
{
    // A threshold of 0 is allowed and means every run, including a run of one
    // base, becomes a poly symbol. The symbol string is then the run length
    // encoding of the sequence and carries no length information at all.
    //
    // That is not the same as plain RLE, which is what rle.hpp does and what
    // this alphabet exists to improve on. Plain RLE throws the lengths away.
    // Here they move into the second half of the pair and are voted on, so the
    // representation stays lossless.

    encoded.clear();

    // The encoding never has more symbols than the sequence has bases.
    encoded.reserve(sequence.size());

    const uint64_t n = sequence.size();
    for(uint64_t runBegin=0; runBegin<n; /* incremented below */) {
        const Base base = sequence[runBegin];

        // Find the end of the homopolymer run that begins here.
        uint64_t runEnd = runBegin + 1;
        while((runEnd < n) and (sequence[runEnd] == base)) {
            ++runEnd;
        }
        const uint64_t runLength = runEnd - runBegin;

        if(runLength > threshold) {

            // A long run becomes a single poly symbol standing for the whole
            // run. Its length is carried alongside and plays no part in any
            // alignment done on the symbols.
            encoded.push_back(make_pair(ExtendedBase::polyFromBase(base), runLength));

        } else {

            // A short run keeps one plain symbol per base, so its length is
            // still visible to an aligner. Short run lengths are reliable and
            // are real signal, unlike long ones.
            for(uint64_t i=0; i<runLength; i++) {
                encoded.push_back(make_pair(ExtendedBase::fromBase(base), 1UL));
            }
        }

        runBegin = runEnd;
    }
}



// See msaRepair.hpp for comments.
void shasta2::decodeExtended(
    const ExtendedSequence& encoded,
    vector<Base>& sequence)
{
    sequence.clear();

    uint64_t n = 0;
    for(const auto& [extendedBase, runLength]: encoded) {
        n += runLength;
    }
    sequence.reserve(n);

    for(const auto& [extendedBase, runLength]: encoded) {
        const Base base = extendedBase.base();
        for(uint64_t j=0; j<runLength; j++) {
            sequence.push_back(base);
        }
    }
    SHASTA2_ASSERT(sequence.size() == n);
}



// See msaRepair.hpp for comments.
void shasta2::attachRunLengths(
    const vector<AlignedExtendedBase>& alignedRow,
    const ExtendedSequence& encoded,
    Anchoring anchoring,
    AlignedExtendedSequence& row)
{
    // Count the non-gap symbols of the row.
    uint64_t rowSymbolCount = 0;
    for(const AlignedExtendedBase e: alignedRow) {
        if(not e.isGap()) {
            ++rowSymbolCount;
        }
    }
    SHASTA2_ASSERT(rowSymbolCount <= encoded.size());

    // Locate the part of the encoding that this row covers.
    // Which end theseus may have trimmed is known from the anchoring and is not
    // inferred: see the comment in msaRepair.hpp for why searching the encoding for
    // the row's symbols gives the wrong answer on a repeat.
    uint64_t offset = 0;
    switch(anchoring) {
    case Anchoring::BothSides:
        SHASTA2_ASSERT(rowSymbolCount == encoded.size());
        offset = 0;
        break;
    case Anchoring::LeftOnly:
        // Anchored on the left, trimmed on the right: a prefix.
        offset = 0;
        break;
    case Anchoring::RightOnly:
        // Anchored on the right, trimmed on the left: a suffix.
        offset = encoded.size() - rowSymbolCount;
        break;
    }

    row.clear();
    row.reserve(alignedRow.size());
    uint64_t i = 0;
    for(const AlignedExtendedBase e: alignedRow) {
        if(e.isGap()) {
            row.push_back(make_pair(e, 0UL));
        } else {
            const auto& [encodedBase, runLength] = encoded[offset + i];
            SHASTA2_ASSERT(e == AlignedExtendedBase(encodedBase));
            row.push_back(make_pair(e, runLength));
            ++i;
        }
    }
    SHASTA2_ASSERT(i == rowSymbolCount);
}



namespace shasta2 {

    // The vote over the TRUE run lengths of the rows sharing one poly symbol
    // at one alignment column: mode, median and mean, computed together in a
    // single pass over lengthWeight so that every RunLengthEstimator shares
    // one scan instead of recomputing its own.
    class MsaRepairLengthVote {
    public:
        uint64_t modeLength = 0;
        uint64_t modeWeight = 0;
        uint64_t medianLength = 0;
        uint64_t cumulativeAtMedian = 0;
        uint64_t weightAtMedian = 0;
        uint64_t meanLength = 0;

        MsaRepairLengthVote(const vector<uint64_t>& lengthWeight, uint64_t maxObserved, uint64_t totalWeight)
        {
            uint64_t cumulative = 0;
            uint64_t weightedSum = 0;
            bool medianFound = false;
            for(uint64_t length=1; length<=maxObserved; length++) {
                const uint64_t weight = lengthWeight[length];

                // The mode: the most frequent length. Ties go to the shorter
                // run, which the strict > gives us scanning in increasing
                // order.
                if(weight > modeWeight) {
                    modeWeight = weight;
                    modeLength = length;
                }

                // The weighted median: the smallest length whose running
                // weight crosses 50%.
                cumulative += weight;
                if((not medianFound) and (2 * cumulative >= totalWeight)) {
                    medianLength = length;
                    cumulativeAtMedian = cumulative;
                    weightAtMedian = weight;
                    medianFound = true;
                }

                weightedSum += length * weight;
            }
            SHASTA2_ASSERT(medianFound);

            // The weighted mean, rounded to the nearest integer, ties rounding up.
            meanLength = (weightedSum + (totalWeight / 2)) / totalWeight;
            if(meanLength == 0) {
                meanLength = 1;
            }
        }
    };
}



// See msaRepair.hpp for comments. This overload infers the spans from the gaps.
void shasta2::extendedConsensus(
    const vector<AlignedExtendedSequence>& alignment,
    const vector<uint64_t>& weights,
    RunLengthEstimator estimator,
    vector< pair<Base, uint64_t> >& consensus,
    AlignedExtendedSequence& alignedConsensus)
{
    // The span of each row is the columns between its first and last non-gap
    // symbol. Everything outside that is read as padding.
    vector< pair<uint64_t, uint64_t> > spans;
    spans.reserve(alignment.size());
    for(const AlignedExtendedSequence& row: alignment) {
        spans.push_back(msaRepairBaseSpan(row));
    }

    extendedConsensus(alignment, weights, estimator, spans,
        consensus, alignedConsensus);
}



void shasta2::extendedConsensus(
    const vector<AlignedExtendedSequence>& alignment,
    const vector<uint64_t>& weights,
    RunLengthEstimator estimator,
    const vector< pair<uint64_t, uint64_t> >& spans,
    vector< pair<Base, uint64_t> >& consensus,
    AlignedExtendedSequence& alignedConsensus,
    shared_ptr<const HomopolymerModel> homopolymerModelPointer,
    const vector< array<uint64_t, 2> >& strandWeights)
{
    consensus.clear();
    alignedConsensus.clear();

    const uint64_t n = alignment.size();
    if(n == 0) {
        return;
    }
    SHASTA2_ASSERT(weights.size() == n);
    SHASTA2_ASSERT(spans.size() == n);
    if(homopolymerModelPointer) {
        SHASTA2_ASSERT(strandWeights.size() == n);
        for(uint64_t i=0; i<n; i++) {
            SHASTA2_ASSERT(strandWeights[i][0] + strandWeights[i][1] == weights[i]);
        }
    }

    const uint64_t alignmentLength = alignment.front().size();
    for(uint64_t i=0; i<n; i++) {
        SHASTA2_ASSERT(alignment[i].size() == alignmentLength);
        SHASTA2_ASSERT(spans[i].first <= spans[i].second);
        SHASTA2_ASSERT(spans[i].second <= alignmentLength);
    }

    alignedConsensus.resize(alignmentLength);

    // Not a guaranteed bound in either direction - a gap column contributes no
    // base, while a poly column usually contributes several - but the two
    // effects are both small relative to alignmentLength in practice, so this
    // is a reasonable estimate that avoids most reallocations without
    // claiming to know the exact final size.
    consensus.reserve(alignmentLength);

    // Scratch, reused by the run length vote at every poly column.
    vector<uint64_t> lengthWeight;
    array<vector<uint64_t>, 2> observedLengths;

    // The consensus symbol of column j, the one with the most weight over the
    // rows that cover it, and that weight. A column no row covers gets a gap
    // with weight 0.
    const auto voteSymbol = [&](uint64_t j) {

        // Total weight for each symbol at this column, including the gap.
        array<uint64_t, AlignedExtendedBase::gapValue + 1> symbolWeight;
        fill(symbolWeight.begin(), symbolWeight.end(), 0);
        uint64_t coveringWeight = 0;
        for(uint64_t i=0; i<n; i++) {
            if((j < spans[i].first) or (j >= spans[i].second)) {
                continue;
            }
            symbolWeight[alignment[i][j].first.value] += weights[i];
            coveringWeight += weights[i];
        }
        if(coveringWeight == 0) {
            return make_pair(AlignedExtendedBase::gap(), 0UL);
        }
        const auto it = std::ranges::max_element(symbolWeight);
        return make_pair(AlignedExtendedBase::fromInteger(uint64_t(it - symbolWeight.begin())), *it);
    };

    // The homopolymer model also uses the bases on each side of a run: for
    // each column, the base of the nearest consensus symbol to its left and to
    // its right that is not a gap, or HomopolymerModel::unknownFlank if there
    // is none.
    vector<uint64_t> leftFlank;
    vector<uint64_t> rightFlank;
    if(homopolymerModelPointer) {
        vector<AlignedExtendedBase> symbols(alignmentLength);
        for(uint64_t j=0; j<alignmentLength; j++) {
            symbols[j] = voteSymbol(j).first;
        }
        leftFlank.assign(alignmentLength, HomopolymerModel::unknownFlank);
        rightFlank.assign(alignmentLength, HomopolymerModel::unknownFlank);
        uint64_t flank = HomopolymerModel::unknownFlank;
        for(uint64_t j=0; j<alignmentLength; j++) {
            leftFlank[j] = flank;
            if(not symbols[j].isGap()) {
                flank = symbols[j].base().value;
            }
        }
        flank = HomopolymerModel::unknownFlank;
        for(uint64_t j=alignmentLength; j>0; j--) {
            rightFlank[j - 1] = flank;
            if(not symbols[j - 1].isGap()) {
                flank = symbols[j - 1].base().value;
            }
        }
    }

    // Loop over alignment columns.
    for(uint64_t j=0; j<alignmentLength; j++) {
        const auto [consensusSymbol, symbolCoverage] = voteSymbol(j);

        // A column no row covers contributes nothing, and neither does a gap.
        if(consensusSymbol.isGap()) {
            alignedConsensus[j] = make_pair(consensusSymbol, 0UL);
            continue;
        }

        // The weight supporting this symbol.
        uint64_t coverage = symbolCoverage;

        // The consensus run length. A plain symbol always stands for one base.
        // A poly symbol gets its length by voting on the TRUE run lengths of the
        // rows that carry the same poly symbol here. This is the whole point of
        // the extended alphabet: the length is decided here, not by counting how
        // many columns the aligner happened to fill.
        uint64_t consensusRunLength = 1;
        if(consensusSymbol.isPoly()) {

            // Gather the observed lengths with their weights.
            // Lengths are bounded by the sequence length, so a vector indexed by
            // length is both simpler and faster than a map. It is reused across
            // columns rather than reallocated at each one, because a long
            // sequence has many poly columns and each allocation would be
            // proportional to a run length.
            // Only rows that cover this column may vote on the length, for the
            // same reason they may not vote on the symbol.
            uint64_t maxObserved = 0;
            uint64_t totalWeight = 0;
            for(uint64_t i=0; i<n; i++) {
                if((j < spans[i].first) or (j >= spans[i].second)) {
                    continue;
                }
                if(alignment[i][j].first == consensusSymbol) {
                    maxObserved = max(maxObserved, alignment[i][j].second);
                    totalWeight += weights[i];
                }
            }
            SHASTA2_ASSERT(totalWeight > 0);

            if(homopolymerModelPointer) {

                // The most likely length under the model, given the length
                // each covering row reports on each strand. A row that
                // stands for several reads enters its length once per read.
                for(vector<uint64_t>& v: observedLengths) {
                    v.clear();
                }
                for(uint64_t i=0; i<n; i++) {
                    if((j < spans[i].first) or (j >= spans[i].second)) {
                        continue;
                    }
                    if(alignment[i][j].first == consensusSymbol) {
                        for(uint64_t strand=0; strand<2; strand++) {
                            observedLengths[strand].insert(observedLengths[strand].end(),
                                strandWeights[i][strand], alignment[i][j].second);
                        }
                    }
                }
                consensusRunLength = homopolymerModelPointer->mostLikelyLength(
                    consensusSymbol.base(), leftFlank[j], rightFlank[j], observedLengths);
                coverage = totalWeight;

            } else {

                lengthWeight.assign(maxObserved + 1, 0);
                for(uint64_t i=0; i<n; i++) {
                    if((j < spans[i].first) or (j >= spans[i].second)) {
                        continue;
                    }
                    if(alignment[i][j].first == consensusSymbol) {
                        lengthWeight[alignment[i][j].second] += weights[i];
                    }
                }

                // One shared pass computes mode, median and mean together; the
                // estimator just picks which of them (or which simple function of
                // them) to use.
                const MsaRepairLengthVote vote(lengthWeight, maxObserved, totalWeight);

                switch(estimator) {

                case RunLengthEstimator::Mode:
                    consensusRunLength = vote.modeLength;
                    coverage = vote.modeWeight;
                    break;

                case RunLengthEstimator::Median:
                    consensusRunLength = vote.medianLength;
                    coverage = vote.weightAtMedian;
                    break;

                case RunLengthEstimator::MedianMarginGated:
                    // Nudge the median up by one only when its cumulative support
                    // is not comfortably (>=60%) above 50% - a near-tie between
                    // "at least this long" and "shorter" - capped at maxObserved:
                    // never invent a length longer than what some row actually
                    // reports. See the comment on this estimator in the header.
                    consensusRunLength = (10 * vote.cumulativeAtMedian >= 6 * totalWeight) ?
                        vote.medianLength : min(vote.medianLength + 1, maxObserved);
                    coverage = vote.weightAtMedian;
                    break;

                case RunLengthEstimator::Average:
                    // The weighted mean. Every observed length contributes, so
                    // this uses the whole distribution rather than a single
                    // order statistic.
                    consensusRunLength = vote.meanLength;
                    coverage = totalWeight;
                    break;
                }
            }
            SHASTA2_ASSERT(consensusRunLength > 0);
        }

        alignedConsensus[j] = make_pair(consensusSymbol, consensusRunLength);

        // Store the run in the consensus.
        const Base base = consensusSymbol.base();
        for(uint64_t k=0; k<consensusRunLength; k++) {
            consensus.push_back(make_pair(base, coverage));
        }
    }
}



// See msaRepair.hpp for comments.
void shasta2::expandExtendedAlignment(
    const vector<AlignedExtendedSequence>& alignment,
    const AlignedExtendedSequence& alignedConsensus,
    vector< vector<AlignedBase> >& expandedAlignment,
    vector<AlignedBase>& expandedAlignedConsensus)
{
    expandedAlignment.clear();
    expandedAlignedConsensus.clear();

    const uint64_t n = alignment.size();
    if(n == 0) {
        return;
    }
    const uint64_t alignmentLength = alignment.front().size();
    SHASTA2_ASSERT(alignedConsensus.size() == alignmentLength);

    // The width of each column after expansion: the longest run seen there,
    // including the consensus, and at least 1.
    vector<uint64_t> width(alignmentLength, 1);
    uint64_t expandedLength = 0;
    for(uint64_t j=0; j<alignmentLength; j++) {
        for(uint64_t i=0; i<n; i++) {
            width[j] = max(width[j], alignment[i][j].second);
        }
        width[j] = max(width[j], alignedConsensus[j].second);
        expandedLength += width[j];
    }

    // Expand one row, padding on the right with gaps.
    // Every expanded row has the same length, known here, so it is reserved
    // once instead of being grown a base at a time.
    const auto expand = [&width, alignmentLength, expandedLength](
        const AlignedExtendedSequence& row,
        vector<AlignedBase>& expanded)
    {
        expanded.clear();
        expanded.reserve(expandedLength);
        for(uint64_t j=0; j<alignmentLength; j++) {
            const auto& [e, runLength] = row[j];
            for(uint64_t k=0; k<runLength; k++) {
                expanded.push_back(AlignedBase(e.base()));
            }
            for(uint64_t k=runLength; k<width[j]; k++) {
                expanded.push_back(AlignedBase::gap());
            }
        }
        SHASTA2_ASSERT(expanded.size() == expandedLength);
    };

    expandedAlignment.resize(n);
    for(uint64_t i=0; i<n; i++) {
        expand(alignment[i], expandedAlignment[i]);
    }
    expand(alignedConsensus, expandedAlignedConsensus);

    // All rows and the aligned consensus must have the same length.
    for(const vector<AlignedBase>& row: expandedAlignment) {
        SHASTA2_ASSERT(row.size() == expandedAlignedConsensus.size());
    }
}



// See msaRepair.hpp for comments.
bool shasta2::msaRepairPatternPresent(
    const vector<Base>& sequence,
    uint64_t threshold)
{
    // Walk the runs, keeping only the two previous run lengths. The pattern is a
    // run longer than the threshold, then a run of exactly one base, then
    // another run longer than the threshold. Nothing is allocated and the walk
    // stops at the first hit, which is what lets this be called on every local
    // assembly before deciding whether an alignment is even needed.
    uint64_t lengthBeforeLast = 0;
    uint64_t lastLength = 0;

    const uint64_t n = sequence.size();
    for(uint64_t runBegin=0; runBegin<n; /* incremented below */) {
        const Base base = sequence[runBegin];
        uint64_t runEnd = runBegin + 1;
        while((runEnd < n) and (sequence[runEnd] == base)) {
            ++runEnd;
        }
        const uint64_t runLength = runEnd - runBegin;

        if( (lengthBeforeLast > threshold) and
            (lastLength == 1) and
            (runLength > threshold)) {
            return true;
        }

        lengthBeforeLast = lastLength;
        lastLength = runLength;
        runBegin = runEnd;
    }
    return false;
}



// See msaRepair.hpp for comments.
bool shasta2::msaRepairLongRunPresent(
    const vector<Base>& sequence,
    uint64_t threshold)
{
    const uint64_t n = sequence.size();
    for(uint64_t runBegin=0; runBegin<n; /* incremented below */) {
        const Base base = sequence[runBegin];
        uint64_t runEnd = runBegin + 1;
        while((runEnd < n) and (sequence[runEnd] == base)) {
            ++runEnd;
        }
        if(runEnd - runBegin > threshold) {
            return true;
        }
        runBegin = runEnd;
    }
    return false;
}



// See msaRepair.hpp for comments.
bool shasta2::msaRepairTriggerPresent(
    const vector<Base>& sequence,
    MsaRepairTrigger trigger,
    uint64_t threshold)
{
    if(trigger == MsaRepairTrigger::PatternOnly) {
        return msaRepairPatternPresent(sequence, threshold);
    }
    return msaRepairLongRunPresent(sequence, threshold);
}



// See msaRepair.hpp for comments.
bool shasta2::msaRepairTriggerPresent(
    const vector< vector<Base> >& sequences,
    MsaRepairTrigger trigger,
    uint64_t threshold)
{
    return std::ranges::any_of(sequences, [&](const vector<Base>& sequence) {
        return msaRepairTriggerPresent(sequence, trigger, threshold);
    });
}



// See msaRepair.hpp for comments.
bool shasta2::msaRepairTriggerPresent(
    const vector< pair<vector<Base>, uint64_t> >& sequences,
    MsaRepairTrigger trigger,
    uint64_t threshold)
{
    return std::ranges::any_of(sequences, [&](const auto& sequenceAndCoverage) {
        return msaRepairTriggerPresent(sequenceAndCoverage.first, trigger, threshold);
    });
}



// See msaRepair.hpp for comments.
string shasta2::msaRepairTriggerDescription(MsaRepairTrigger trigger)
{
    return (trigger == MsaRepairTrigger::PatternOnly) ?
        "a long homopolymer run bordered by a single base" :
        "a long homopolymer run";
}



namespace shasta2 {

    // Return true if a region boundary at this column falls inside a homopolymer
    // run of any row, that is, if the last base of some row before the boundary
    // is the same as its next base at or after it.
    //
    // A boundary that cuts a run is a problem for everything downstream: the
    // pattern test would see a run shorter than it is and might not recognise
    // the pattern at all, and the encoding would collapse a partial run.
    static bool msa1BoundaryCutsRun(
        const vector< vector<AlignedBase> >& alignment,
        uint64_t column)
    {
        for(const vector<AlignedBase>& row: alignment) {

            // The last base before the boundary.
            AlignedBase before = AlignedBase::gap();
            for(uint64_t k=column; k>0; k--) {
                if(not row[k - 1].isGap()) {
                    before = row[k - 1];
                    break;
                }
            }
            if(before.isGap()) {
                continue;
            }

            // The first base at or after the boundary.
            AlignedBase after = AlignedBase::gap();
            for(uint64_t k=column; k<row.size(); k++) {
                if(not row[k].isGap()) {
                    after = row[k];
                    break;
                }
            }
            if(after.isGap()) {
                continue;
            }

            if(before == after) {
                return true;
            }
        }
        return false;
    }
}



// See msaRepair.hpp for comments.
void shasta2::msaRepairRowCoverage(
    const vector< vector<AlignedBase> >& alignment,
    const vector<Anchoring>& anchoring,
    vector< pair<uint64_t, uint64_t> >& coverage)
{
    coverage.clear();
    if(alignment.empty()) {
        return;
    }
    const uint64_t alignmentLength = alignment.front().size();

    // Exactly one entry is pushed per row.
    coverage.reserve(alignment.size());

    // An empty anchoring argument means every row is fixed on both sides.
    SHASTA2_ASSERT(anchoring.empty() or (anchoring.size() == alignment.size()));

    for(uint64_t i=0; i<alignment.size(); i++) {
        const Anchoring rowAnchoring =
            anchoring.empty() ? Anchoring::BothSides : anchoring[i];

        // A row fixed on both sides covers every column whatever its bases are,
        // so the scan is only worth doing for the other two. On the abpoa path
        // every row is fixed on both sides and this skips a full pass.
        if(rowAnchoring == Anchoring::BothSides) {
            // Covers every column, gaps included. A row with no base at all
            // still covers them: a read that is empty here has the whole thing
            // deleted, and that is evidence, not absence.
            coverage.push_back(make_pair(0UL, alignmentLength));
            continue;
        }
        const pair<uint64_t, uint64_t> baseSpan = msaRepairBaseSpan(alignment[i]);
        const bool noBases = (baseSpan.first == baseSpan.second);

        switch(rowAnchoring) {

        case Anchoring::BothSides:
            // Handled above.
            SHASTA2_ASSERT(false);
            break;

        case Anchoring::LeftOnly:
            // Covers up to its last base. Where its last base is has no answer
            // if it has none, and a read with no base is not constrained to
            // anything, so it covers nothing.
            coverage.push_back(noBases ?
                make_pair(0UL, 0UL) : make_pair(0UL, baseSpan.second));
            break;

        case Anchoring::RightOnly:
            coverage.push_back(noBases ?
                make_pair(0UL, 0UL) : make_pair(baseSpan.first, alignmentLength));
            break;
        }
    }
}



void shasta2::msaRepairFindBadRegions(
    const vector< vector<AlignedBase> >& alignment,
    const vector<AlignedBase>& alignedConsensus,
    MsaRepairTrigger trigger,
    uint64_t threshold,
    uint64_t flank,
    uint64_t mergeDistance,
    const vector< pair<uint64_t, uint64_t> >& coverageArgument,
    vector<MsaRepairRegion>& regions)
{
    regions.clear();
    if(alignment.empty()) {
        return;
    }
    const uint64_t alignmentLength = alignedConsensus.size();
    for(const vector<AlignedBase>& row: alignment) {
        SHASTA2_ASSERT(row.size() == alignmentLength);
    }

    // An empty coverage argument means every row covers the whole alignment.
    vector< pair<uint64_t, uint64_t> > wholeAlignment;
    if(coverageArgument.empty()) {
        wholeAlignment.assign(alignment.size(), make_pair(0UL, alignmentLength));
    }
    const vector< pair<uint64_t, uint64_t> >& coverage =
        coverageArgument.empty() ? wholeAlignment : coverageArgument;
    SHASTA2_ASSERT(coverage.size() == alignment.size());

    // Seed a region on each long homopolymer run of each row.
    //
    // The seeds are the runs themselves, not the columns where the rows disagree
    // with the consensus. Disagreement is not a selective signal here: reads
    // differ in run length, so their gaps land in different columns, and on real
    // reads 95% of columns end up with some row differing from the consensus.
    // Seeding on disagreement and then filtering by the pattern worked only
    // because the pattern is rare; with the broader trigger it selected almost
    // the whole alignment. Seeding on the runs keeps a region the size of what
    // it is repairing whichever trigger is used.
    //
    // A run of a row may have gap columns interleaved, where other rows have an
    // insertion, so the run is tracked by its first and last non-gap column.
    for(const vector<AlignedBase>& row: alignment) {
        AlignedBase runBase = AlignedBase::gap();
        uint64_t runBegin = 0;
        uint64_t runLast = 0;
        uint64_t runLength = 0;

        const auto closeRun = [&]() {
            if((runLength > threshold) and (not runBase.isGap())) {
                MsaRepairRegion seed;
                seed.begin = runBegin;
                seed.end = runLast + 1;
                regions.push_back(seed);
            }
        };

        for(uint64_t j=0; j<alignmentLength; j++) {
            const AlignedBase b = row[j];
            if(b.isGap()) {
                continue;
            }
            if(b == runBase) {
                ++runLength;
                runLast = j;
            } else {
                closeRun();
                runBase = b;
                runBegin = j;
                runLast = j;
                runLength = 1;
            }
        }
        closeRun();
    }
    if(regions.empty()) {
        return;
    }

    // Sort the seeds and merge those that overlap or are close together. Two
    // runs separated by a single base merge here, which is what lets the pattern
    // trigger see both of them in one window.
    sort(regions.begin(), regions.end());
    vector<MsaRepairRegion> merged;
    for(const MsaRepairRegion& seed: regions) {
        if((not merged.empty()) and (seed.begin <= merged.back().end + mergeDistance)) {
            merged.back().end = max(merged.back().end, seed.end);
        } else {
            merged.push_back(seed);
        }
    }

    // Widen by the flank, then widen further until neither boundary falls inside
    // a homopolymer run of any row. A window that clips a run would show the
    // pattern test and the encoding a run shorter than it is.
    for(MsaRepairRegion& region: merged) {
        region.begin = (region.begin > flank) ? (region.begin - flank) : 0;
        region.end = min(alignmentLength, region.end + flank);
        while((region.begin > 0) and msa1BoundaryCutsRun(alignment, region.begin)) {
            --region.begin;
        }
        while((region.end < alignmentLength) and msa1BoundaryCutsRun(alignment, region.end)) {
            ++region.end;
        }
    }

    // Widening may have made neighbours touch.
    regions.clear();
    for(const MsaRepairRegion& region: merged) {
        if((not regions.empty()) and (region.begin <= regions.back().end)) {
            regions.back().end = max(regions.back().end, region.end);
        } else {
            regions.push_back(region);
        }
    }

    // Keep only the regions the trigger accepts. For the pattern trigger this is
    // what does the filtering; for the broader trigger every seed already has a
    // long run in it, so this keeps them all.
    //
    // A region in which no row differs from the consensus is dropped whatever
    // the trigger says. Every row already agrees there, so every vote would
    // return what is already in place and the repair could only be a no-op. The
    // seeds are the runs rather than the disagreements, so unlike before this
    // has to be checked rather than being true by construction.
    vector<MsaRepairRegion> confirmed;
    vector<Base> windowSequence;
    for(const MsaRepairRegion& region: regions) {

        // A row that does not reach this column is not disagreeing with the
        // consensus there, it is simply absent. Counting the padding of a read
        // constrained on one side only as disagreement would make every region
        // beyond the end of some read look worth repairing.
        bool anyDiscordant = false;
        for(uint64_t j=region.begin; (j<region.end) and (not anyDiscordant); j++) {
            for(uint64_t i=0; i<alignment.size(); i++) {
                if(alignment[i][j] == alignedConsensus[j]) {
                    continue;
                }
                if((j < coverage[i].first) or (j >= coverage[i].second)) {
                    continue;
                }
                anyDiscordant = true;
                break;
            }
        }
        if(not anyDiscordant) {
            continue;
        }

        bool found = false;
        for(const vector<AlignedBase>& row: alignment) {
            windowSequence.clear();
            windowSequence.reserve(region.end - region.begin);
            for(uint64_t j=region.begin; j<region.end; j++) {
                if(not row[j].isGap()) {
                    windowSequence.push_back(Base(row[j]));
                }
            }
            if(msaRepairTriggerPresent(windowSequence, trigger, threshold)) {
                found = true;
                break;
            }
        }
        if(found) {
            confirmed.push_back(region);
        }
    }
    regions = std::move(confirmed);
}



namespace shasta2 {

    // Align sequences encoded in the extended alphabet, using theseus
    // directly rather than through theseusWrapper.cpp: duplicates a small
    // amount of its plumbing (build a TheseusMSA aligner, align each
    // sequence, read the alignment back out of print_as_msa) so the whole
    // msaRepair feature stays readable and changeable from this one file, and a
    // change here cannot affect theseus()'s own callers or vice versa.
    // Theseus needs no change to accept the extended alphabet: it only ever
    // compares two characters with ==, so it passes through as the string
    // ACGTacgt (see ExtendedBase in msaRepair.hpp).
    //
    // Only the alignment is computed, not the consensus: a consensus over
    // the extended alphabet also needs the run lengths, which theseus does
    // not have. See extendedConsensus for that.
    static void msaRepairAlignExtended(
        const vector< pair<ExtendedSequence, uint64_t> >& fixedSequences,
        const vector< pair<ExtendedSequence, uint64_t> >& leftFixedSequences,
        const vector< pair<ExtendedSequence, uint64_t> >& rightFixedSequences,
        vector< vector<AlignedExtendedBase> >& alignment)
    {
        // Pericles default penalties, the same theseus() uses.
        const int match = 0;
        const int mismatch = 2;
        const int gapo = 3;
        const int gape = 1;
        const theseus::Penalties penalties(match, mismatch, gapo, gape);
        const theseus::Heuristics heuristics;

        // Create the aligner, passing in the first sequence fixed on both sides.
        SHASTA2_ASSERT(not fixedSequences.empty());
        const auto& [firstSequence, firstSequenceWeight] = fixedSequences.front();
        theseus::TheseusMSA aligner(penalties, heuristics,
            toString(firstSequence), int(firstSequenceWeight), false);

        // Pass in the remaining sequences fixed on both sides.
        for(uint64_t i=1; i<fixedSequences.size(); i++) {
            const auto& [sequence, weight] = fixedSequences[i];
            const bool densityDrop = false;
            const bool lagPruning = false;
            aligner.align(toString(sequence), int(weight), false, false, densityDrop, lagPruning);
        }

        // Pass in the sequences fixed on the left only.
        for(const auto& [sequence, weight]: leftFixedSequences) {
            const bool densityDrop = true;
            const bool lagPruning = false;
            aligner.align(toString(sequence), int(weight), false, true, densityDrop, lagPruning);
        }

        // Pass in the sequences fixed on the right only.
        for(const auto& [sequence, weight]: rightFixedSequences) {
            const bool densityDrop = true;
            const bool lagPruning = false;
            aligner.align(toString(sequence), int(weight), true, true, densityDrop, lagPruning);
        }

        // Read the alignment back. print_as_msa writes a header line for each
        // sequence, then its alignment row, and finally a row for the aligned
        // consensus, which is not wanted here.
        alignment.clear();
        std::ostringstream s;
        aligner.print_as_msa(s);
        const string alignmentString = std::move(s).str();

        boost::tokenizer< boost::char_separator<char> > tokenizer(
            alignmentString, boost::char_separator<char>("\n"));
        for(const string& line: tokenizer) {
            SHASTA2_ASSERT(not line.empty());
            if(line[0] != '>') {

                // Parse with AlignedExtendedBase, NOT AlignedBase.
                // AlignedBase::fromCharacter maps 'a' to A, which would
                // silently discard the poly distinction and turn every long
                // run back into a plain base.
                alignment.push_back(vectorOfAlignedExtendedBasesFromString(line));
            }
        }

        SHASTA2_ASSERT(not alignment.empty());
        alignment.pop_back();
        SHASTA2_ASSERT(alignment.size() ==
            fixedSequences.size() + leftFixedSequences.size() + rightFixedSequences.size());

        // All rows must have the same length.
        for(const vector<AlignedExtendedBase>& row: alignment) {
            SHASTA2_ASSERT(row.size() == alignment.front().size());
        }
    }



    // Recompute one region of an alignment using the extended alphabet, and
    // return the replacement columns.
    //
    // The rows of the region are ungapped and encoded. They are repaired only if
    // their symbol strings all agree, in which case they stack without any
    // alignment and the usual two votes run, symbols first and run lengths
    // second, before expanding back to bases.
    //
    // Returns false and leaves the outputs alone if the region cannot be
    // How one row of the alignment sits in the window being repaired.
    //
    // Everything the repair needs to know about a row follows from these, so
    // they are kept together rather than as index lists that have to be built
    // in one place and read in three more.
    class MsaRepairRepairRow {
    public:

        // The window's bases for this row, ungapped and encoded.
        ExtendedSequence encoding;

        // Which ends of the window the row reaches. Meaningless when isEmpty.
        Anchoring anchoring = Anchoring::BothSides;

        // The row has no base at all in the window.
        bool isEmpty = false;

        // Only for an empty row: it covers the window, so its absence is a
        // deletion it may vote for rather than a read that is not here.
        bool coversWindow = false;
    };

    // improved, in which case the caller must leave it as it was.
    static bool msa1RepairRegion(
        const vector< vector<AlignedBase> >& alignment,
        const MsaRepairRegion& region,
        const vector<uint64_t>& weights,
        const vector< pair<uint64_t, uint64_t> >& coverage,
        uint64_t encodeThreshold,
        RunLengthEstimator estimator,
        shared_ptr<const HomopolymerModel> homopolymerModelPointer,
        const vector< array<uint64_t, 2> >& strandWeights,
        vector< vector<AlignedBase> >& newRows,
        vector<AlignedBase>& newAlignedConsensus,
        vector< pair<Base, uint64_t> >& newConsensus)
    {
        const uint64_t n = alignment.size();

        // Ungap each row inside the region, encode it, and record how it sits in
        // the window.
        //
        // A read constrained on one side only does not reach across the whole
        // alignment. A row can therefore reach one edge of the window and not
        // the other, and that has to be passed on: theseus is told which end a
        // sequence is anchored at, and telling it a right anchored fragment is
        // anchored on both sides invites it to place the fragment against the
        // wrong part of the window. On a window holding A12 G A11, with two
        // reads covering only the A11, calling them fixed on both sides made
        // theseus align them to the FIRST run and the vote returned 11 for it
        // instead of 12.
        // A row with no base at all in the window has nothing to align, so it is
        // left out of the realignment and comes back as the gaps it already was,
        // which loses nothing.
        //
        // Whether it votes is a different question, and the answer is not the
        // same for all of them. A row that covers the window and has no base in
        // it has the window deleted, and a deletion is a thing to vote for. A
        // row that does not reach the window is simply not here. Both look
        // identical in the alignment; the coverage tells them apart.
        vector<MsaRepairRepairRow> rows(n);
        vector<Base> windowSequence;
        windowSequence.reserve(region.end - region.begin);
        for(uint64_t i=0; i<n; i++) {
            windowSequence.clear();
            for(uint64_t j=region.begin; j<region.end; j++) {
                if(not alignment[i][j].isGap()) {
                    windowSequence.push_back(Base(alignment[i][j]));
                }
            }
            MsaRepairRepairRow& row = rows[i];
            encodeExtended(windowSequence, encodeThreshold, row.encoding);

            // Does the row reach each edge of the window? A row that stops short
            // of an edge is padded there and is not anchored at it.
            const bool reachesLeft = region.begin >= coverage[i].first;
            const bool reachesRight = region.end <= coverage[i].second;

            if(windowSequence.empty()) {
                row.isEmpty = true;
                row.coversWindow = reachesLeft and reachesRight;
                continue;
            }

            if(reachesLeft and reachesRight) {
                row.anchoring = Anchoring::BothSides;
            } else if(reachesLeft) {
                row.anchoring = Anchoring::LeftOnly;
            } else if(reachesRight) {
                row.anchoring = Anchoring::RightOnly;
            } else {

                // A row reaching neither edge of the window.
                //
                // Coverage as msa1RowCoverage computes it cannot produce one: a
                // row fixed on the left reaches every left edge, one fixed on
                // the right reaches every right edge, and one fixed on both
                // reaches both. So this only arises from a coverage argument
                // that did not come from there, and the window is left alone
                // rather than guessed at.
                //
                // There is no good way to handle it in any case. Theseus has no
                // group for a sequence free at both ends, and putting one in
                // with the spanning rows, where a global alignment would at
                // least keep all of its bases, makes things worse rather than
                // better: the fragment is stretched to reach both ends of the
                // graph, and the graph it corrupts is then used to place every
                // row after it. A fragment ATCGTTTTTGTTTTT came back as
                // AT--CGtGt--, and the next row lost six of its nine symbols.
                return false;
            }
        }

        // The rows of one group, in row order. Theseus is given the groups in
        // this order and returns them in the same order, so this is also the
        // order the run lengths are put back in.
        const auto group = [&](Anchoring anchoring) {
            vector<uint64_t> groupRows;
            for(uint64_t i=0; i<n; i++) {
                if((not rows[i].isEmpty) and (rows[i].anchoring == anchoring)) {
                    groupRows.push_back(i);
                }
            }
            return groupRows;
        };
        const vector<uint64_t> fixedRows = group(Anchoring::BothSides);
        const vector<uint64_t> leftFixedRows = group(Anchoring::LeftOnly);
        const vector<uint64_t> rightFixedRows = group(Anchoring::RightOnly);

        // Theseus needs at least one sequence anchored on both sides to build
        // the graph on. Without one there is nothing to align the rest against.
        if(fixedRows.empty()) {
            return false;
        }

        // Do all rows encode to the same symbol string?
        //
        // When they do, the homopolymer lengths are the only thing that differs
        // between the reads in this window. The rows stack without any alignment
        // at all, and that stacking is not an approximation: identical symbol
        // strings have exactly one gap free alignment and it is all matches, so
        // no alignment can score higher. The repair is then exact, because there
        // is no aligner left to make a placement decision and the defect being
        // repaired cannot recur.
        //
        // When they do not, the window holds real variation as well as the
        // homopolymers, and the rows have to be realigned. Leaving such a window
        // alone is not an option: a single substitution anywhere in it breaks
        // the symbol strings, and on reads with a 1% substitution rate that
        // abandons 97% of the regions found. The repair would almost never fire.
        // This can only be answered for the rows that span the window; a row
        // that covers part of it says nothing about the rest.
        bool uniform = leftFixedRows.empty() and rightFixedRows.empty() and
            std::ranges::none_of(rows, &MsaRepairRepairRow::isEmpty);
        // Compared symbol by symbol rather than by building strings. The run
        // lengths are deliberately not part of this comparison, and comparing
        // only the symbols states that, as well as allocating nothing in what is
        // the innermost loop of the repair.
        if(uniform) {
            const auto sameSymbols = [](const ExtendedSequence& x, const ExtendedSequence& y) {
                return std::ranges::equal(x, y, {},
                    [](const auto& p) { return p.first; },
                    [](const auto& p) { return p.first; });
            };
            for(uint64_t i=1; i<n; i++) {
                if(not sameSymbols(rows[i].encoding, rows[0].encoding)) {
                    uniform = false;
                    break;
                }
            }
        }

        vector<AlignedExtendedSequence> extendedAlignment(n);
        if(uniform) {
            for(uint64_t i=0; i<n; i++) {
                extendedAlignment[i].clear();
                extendedAlignment[i].reserve(rows[i].encoding.size());
                for(const auto& [symbol, runLength]: rows[i].encoding) {
                    extendedAlignment[i].push_back(
                        make_pair(AlignedExtendedBase(symbol), runLength));
                }
            }
        } else {

            // Realign the window in symbol space, where a long run is a single
            // symbol and a run length difference therefore cannot be traded
            // against a mismatch. This is what abpoa cannot do, being fixed at
            // four letters, and it is the whole reason the alphabet exists.
            //
            // Theseus is used because it compares characters and so accepts the
            // nine symbol alphabet unchanged. It is the weaker aligner of the
            // two on whole sequences, but the comparison that matters here is
            // not theseus against abpoa: it is a window realigned by theseus in
            // symbol space against the same window left as abpoa placed it, with
            // the misplaced base still in it. The window is also short and its
            // rows are nearly identical, which is the regime where the choice of
            // aligner matters least.
            //
            // Each row goes into the group matching how it sits in the window,
            // so a read that reaches only one edge is aligned as anchored at
            // that edge rather than as spanning.
            const auto encodedGroup = [&](const vector<uint64_t>& groupRows) {
                vector< pair<ExtendedSequence, uint64_t> > encoded;
                encoded.reserve(groupRows.size());
                for(const uint64_t i: groupRows) {
                    encoded.push_back(make_pair(rows[i].encoding, weights[i]));
                }
                return encoded;
            };
            vector< vector<AlignedExtendedBase> > alignedSymbols;
            msaRepairAlignExtended(encodedGroup(fixedRows), encodedGroup(leftFixedRows),
                encodedGroup(rightFixedRows), alignedSymbols);
            SHASTA2_ASSERT(alignedSymbols.size() ==
                fixedRows.size() + leftFixedRows.size() + rightFixedRows.size());

            // Theseus returns the rows in the order the groups were given.
            //
            // A sequence anchored on one side only is aligned semi globally, so
            // theseus is free to leave the far end out of the alignment when it
            // does not fit. Those symbols are then simply not in the row, and
            // writing the row back would delete bases from a read. The repair is
            // an alignment change and must never do that, so a trimmed row
            // abandons the whole region and the aligner's own placement stands.
            uint64_t k = 0;
            bool trimmed = false;
            const auto take = [&](const vector<uint64_t>& groupRows, Anchoring anchoring) {
                for(const uint64_t i: groupRows) {
                    const uint64_t symbolCount = uint64_t(std::ranges::count_if(
                        alignedSymbols[k],
                        [](AlignedExtendedBase e) { return not e.isGap(); }));
                    if(symbolCount != rows[i].encoding.size()) {
                        trimmed = true;
                        return;
                    }
                    attachRunLengths(alignedSymbols[k], rows[i].encoding, anchoring,
                        extendedAlignment[i]);
                    ++k;
                }
            };
            take(fixedRows, Anchoring::BothSides);
            if(not trimmed) {
                take(leftFixedRows, Anchoring::LeftOnly);
            }
            if(not trimmed) {
                take(rightFixedRows, Anchoring::RightOnly);
            }
            if(trimmed) {
                return false;
            }
            SHASTA2_ASSERT(k == alignedSymbols.size());

            // A row left out of the realignment is all gaps here. It had no base
            // in this window to begin with, so nothing is lost.
            const uint64_t length = extendedAlignment[fixedRows.front()].size();
            for(uint64_t i=0; i<n; i++) {
                if(rows[i].isEmpty) {
                    extendedAlignment[i].assign(length,
                        make_pair(AlignedExtendedBase::gap(), 0UL));
                }
            }
        }

        // The columns each row is entitled to vote in.
        //
        // This is known here and must not be left to be guessed from the gaps: a
        // row anchored at both ends covers the whole window even if it starts
        // with a gap, and that leading gap is a deletion the row is voting for.
        //
        // A row that is not anchored at an end stops there because the read
        // stops, and the read is free to stop in the middle of a homopolymer
        // run. The symbol it stops on is then a partial observation: its length
        // is the part of the run that happens to be inside the read, not the
        // length of the run. Letting it vote drags the consensus run length down
        // towards the shortest overlap. On a run of 16 with five reads ending
        // inside it, their 5, 12, 15, 7 and 7 outvoted the two reads that
        // actually spanned the run and the consensus came out at 7.
        //
        // So the boundary symbol of an unanchored end does not vote at all,
        // neither on the length nor on the symbol. It is still written into the
        // alignment; it just does not get a say in what the consensus says
        // there. Occasionally the read really does end exactly at a run boundary
        // and a good vote is thrown away, but the two cases are indistinguishable
        // from the alignment and the cost of the wrong guess is not symmetric.
        const uint64_t windowLength = extendedAlignment.front().size();

        // The half open span [first, end), with the boundary symbol dropped at
        // each end that is not anchored. An empty result means the row has
        // nothing left to say and abstains.
        const auto span = [&](uint64_t i, bool anchoredLeft, bool anchoredRight) {
            const auto [baseFirst, baseEnd] = msaRepairBaseSpan(extendedAlignment[i]);
            const uint64_t first = anchoredLeft ? 0 : (baseFirst + 1);
            const uint64_t end = anchoredRight ? windowLength :
                ((baseEnd > 0) ? (baseEnd - 1) : 0);
            if(first >= end) {
                return make_pair(0UL, 0UL);
            }
            return make_pair(first, end);
        };

        vector< pair<uint64_t, uint64_t> > spans(n, make_pair(0UL, windowLength));
        for(uint64_t i=0; i<n; i++) {
            const MsaRepairRepairRow& row = rows[i];
            if(row.isEmpty) {
                spans[i] = row.coversWindow ?
                    make_pair(0UL, windowLength) : make_pair(0UL, 0UL);
            } else if(row.anchoring == Anchoring::LeftOnly) {
                spans[i] = span(i, true, false);
            } else if(row.anchoring == Anchoring::RightOnly) {
                spans[i] = span(i, false, true);
            }
        }

        // Every row must have come out the same length.
        for(const AlignedExtendedSequence& row: extendedAlignment) {
            SHASTA2_ASSERT(row.size() == extendedAlignment.front().size());
        }

        // Vote, then expand.
        AlignedExtendedSequence alignedExtendedConsensus;
        extendedConsensus(extendedAlignment, weights, estimator, spans,
            newConsensus, alignedExtendedConsensus, homopolymerModelPointer, strandWeights);
        expandExtendedAlignment(extendedAlignment, alignedExtendedConsensus,
            newRows, newAlignedConsensus);

        // The bases of the region consensus must be the non-gap entries of the
        // aligned region consensus, in order. The caller splices both, so if
        // they disagreed the two would drift apart.
        uint64_t k = 0;
        for(const AlignedBase b: newAlignedConsensus) {
            if(not b.isGap()) {
                SHASTA2_ASSERT(k < newConsensus.size());
                SHASTA2_ASSERT(newConsensus[k].first == Base(b));
                ++k;
            }
        }
        SHASTA2_ASSERT(k == newConsensus.size());

        return true;
    }
}



// See msaRepair.hpp for comments.
uint64_t shasta2::msaRepair(
    vector< vector<AlignedBase> >& alignment,
    vector<AlignedBase>& alignedConsensus,
    vector< pair<Base, uint64_t> >& consensus,
    const vector<uint64_t>& weights,
    const vector<Anchoring>& anchoring,
    const MsaRepairOptions& options,
    const vector< array<uint64_t, 2> >& strandWeights)
{
    const auto& [trigger, threshold, encodeThreshold, estimator, homopolymerModelPointer,
        flank, mergeDistance] = options;
    const uint64_t n = alignment.size();
    if(n == 0) {
        return 0;
    }
    for(const vector<AlignedBase>& row: alignment) {
        SHASTA2_ASSERT(row.size() == alignedConsensus.size());
    }

    // An empty weights argument means every row weighs 1. Production always
    // passes weights, so the copy is only made when it is really needed.
    vector<uint64_t> unitWeights;
    if(weights.empty()) {
        unitWeights.assign(n, 1);
    }
    const vector<uint64_t>& rowWeights = weights.empty() ? unitWeights : weights;
    SHASTA2_ASSERT(rowWeights.size() == n);

    // The columns each row covers. An empty anchoring argument means every row
    // is fixed on both sides and so covers the whole alignment, which is what
    // abpoa always produces.
    vector< pair<uint64_t, uint64_t> > coverage;
    msaRepairRowCoverage(alignment, anchoring, coverage);

    // Find the regions worth repairing. Usually there are none, and then nothing
    // below runs and nothing is modified.
    vector<MsaRepairRegion> regions;
    msaRepairFindBadRegions(alignment, alignedConsensus, trigger, threshold, flank,
        mergeDistance, coverage, regions);
    if(regions.empty()) {
        return 0;
    }

    // The consensus must agree with the aligned consensus before anything is
    // spliced, because the base positions of a region are found by counting
    // non-gap entries of the aligned consensus.
    {
        uint64_t k = 0;
        for(const AlignedBase b: alignedConsensus) {
            if(not b.isGap()) {
                ++k;
            }
        }
        SHASTA2_ASSERT(k == consensus.size());
    }

    // Repair the regions from right to left, so that the column indices of the
    // regions still to be done, and the base positions of the consensus before
    // them, are not disturbed by the splices already made.
    uint64_t repairedCount = 0;
    for(uint64_t k=regions.size(); k>0; k--) {
        const MsaRepairRegion& region = regions[k - 1];

        vector< vector<AlignedBase> > newRows;
        vector<AlignedBase> newAlignedConsensus;
        vector< pair<Base, uint64_t> > newConsensus;
        if(not msa1RepairRegion(alignment, region, rowWeights, coverage, encodeThreshold,
            estimator, homopolymerModelPointer, strandWeights,
            newRows, newAlignedConsensus, newConsensus)) {
            continue;
        }
        SHASTA2_ASSERT(newRows.size() == n);

        // Keep the repair only if it did not make this window worse.
        //
        // The repair removes the misplaced base from a long run, but the
        // realignment that follows can introduce a smaller version of the same
        // fault elsewhere in the window: theseus prices two mismatches below two
        // gaps, exactly as abpoa does, and on the short runs that stay as plain
        // symbols that trade is still available to it. Measured over 20000
        // random alignments this happened 5 times, always on two row windows
        // where there is no majority to appeal to.
        //
        // Rather than argue about when it is safe, the window is scored before
        // and after by the thing the repair exists to fix, the number of columns
        // holding more than one base, and a repair that loses is discarded. That
        // makes never making the alignment worse a property of the code rather
        // than something the tests have to keep checking for.
        {
            const uint64_t before =
                    msaRepairImpureColumnCount(alignment, region.begin, region.end);
            const uint64_t after = msaRepairImpureColumnCount(newRows);
            if(after > before) {
                continue;
            }
        }

        // Find the range of consensus BASES this region covers, by counting the
        // non-gap entries of the aligned consensus. This has to be done before
        // the splice below, which changes the aligned consensus.
        uint64_t baseBegin = 0;
        for(uint64_t j=0; j<region.begin; j++) {
            if(not alignedConsensus[j].isGap()) {
                ++baseBegin;
            }
        }
        uint64_t baseEnd = baseBegin;
        for(uint64_t j=region.begin; j<region.end; j++) {
            if(not alignedConsensus[j].isGap()) {
                ++baseEnd;
            }
        }
        SHASTA2_ASSERT(baseEnd <= consensus.size());

        // Splice the alignment rows and the aligned consensus, in column space.
        // These must all be spliced together or they stop having equal length.
        for(uint64_t i=0; i<n; i++) {
            vector<AlignedBase>& row = alignment[i];
            row.erase(row.begin() + int64_t(region.begin),
                row.begin() + int64_t(region.end));
            row.insert(row.begin() + int64_t(region.begin),
                newRows[i].begin(), newRows[i].end());
        }
        alignedConsensus.erase(
            alignedConsensus.begin() + int64_t(region.begin),
            alignedConsensus.begin() + int64_t(region.end));
        alignedConsensus.insert(
            alignedConsensus.begin() + int64_t(region.begin),
            newAlignedConsensus.begin(), newAlignedConsensus.end());

        // Splice the consensus, in base space. The same erase and insert, just
        // indexed by base rather than by column.
        //
        // The consensus is spliced rather than recomputed from the repaired
        // alignment. Recomputing it would rewrite the coverage of every position
        // in the sequence, including all the positions outside any repaired
        // region, because the coverage abpoa reports is computed its own way and
        // is not necessarily what a recount of the rows would give. Rewriting
        // them would silently discard abpoa's numbers everywhere, which is
        // exactly what repairing locally is supposed to avoid.
        consensus.erase(
            consensus.begin() + int64_t(baseBegin),
            consensus.begin() + int64_t(baseEnd));
        consensus.insert(
            consensus.begin() + int64_t(baseBegin),
            newConsensus.begin(), newConsensus.end());

        ++repairedCount;
    }

    if(repairedCount == 0) {
        return 0;
    }

    // All rows, the aligned consensus and the consensus must still agree.
    for(const vector<AlignedBase>& row: alignment) {
        SHASTA2_ASSERT(row.size() == alignedConsensus.size());
    }
    {
        uint64_t k = 0;
        for(const AlignedBase b: alignedConsensus) {
            if(not b.isGap()) {
                SHASTA2_ASSERT(k < consensus.size());
                SHASTA2_ASSERT(consensus[k].first == Base(b));
                ++k;
            }
        }
        SHASTA2_ASSERT(k == consensus.size());
    }

    return repairedCount;
}



// The 19 bubble arm sequences this code was developed against.
// They are reads over a single locus with structure CTC A{n1} G A{n2} GTT,
// inside a 2x tandem duplication. Only n1 and n2 vary between them, and the
// TRUE values are known: n1 = 12 and n2 = 11.
namespace shasta2 {
    static const vector<string> msa1TestSequences = {
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAAAAGAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAGAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAGAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAGAAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAGAAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAGAAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
    };

    // The alignment abpoa actually produces for those sequences. Four of its
    // columns hold a mixture of A and G, because the G between the two long A
    // runs was placed in two different columns.
    static const vector<string> msa1BadAlignmentRows = {
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAA-G-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAG-A-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAG-A-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAG-A-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC----AAAAAAAAAAG-A-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAAAG-A-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAAAG-A-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAG-A-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAG-A-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-GAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-GAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAAAA-GAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAAAA-GAAAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---------AAAAAA-G--AAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---------AAAAAA-G--AAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAGAA-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAGAA-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC-------AAAAAAAGAA-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC-------AAAAAAAGAA-AAAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G--AAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTTTCCAGCCTGGGTGACAGAGCGAGACCCCAACTC---AAAAAAAAAAAA-G--AAAAAAAAAAGTTAAACTATAAAGTAAATTCCTCCCATAGTT",
    };

    static string msaRepairToString(const vector<AlignedBase>& v)
    {
        string s;
        for(const AlignedBase b: v) {
            s.push_back(b.character());
        }
        return s;
    }

    // The plain weighted column majority consensus: what an aligner reports and
    // what the repair is handed as input. The coverage of a base is the weight
    // that voted for it.
    //
    // This is deliberately the naive thing the repair exists to improve on, so
    // the tests start from the same place production does.
    static void msa1ColumnConsensus(
        const vector< vector<AlignedBase> >& alignment,
        const vector<uint64_t>& weights,
        vector<AlignedBase>& alignedConsensus,
        vector< pair<Base, uint64_t> >& consensus)
    {
        alignedConsensus.clear();
        consensus.clear();
        for(uint64_t j=0; j<alignment.front().size(); j++) {
            array<uint64_t, 5> weight;
            fill(weight.begin(), weight.end(), 0);
            for(uint64_t i=0; i<alignment.size(); i++) {
                weight[alignment[i][j].value] += weights[i];
            }
            const auto it = std::ranges::max_element(weight);
            const AlignedBase b = AlignedBase::fromInteger(uint64_t(it - weight.begin()));
            alignedConsensus.push_back(b);
            if(not b.isGap()) {
                consensus.push_back(make_pair(Base(b), *it));
            }
        }
    }

    // Options for a test, differing from the production defaults only in the
    // trigger. The tests exercise both triggers; everything else they take as
    // shipped, so a default that changes is exercised rather than bypassed.
    static MsaRepairOptions msa1TestOptions(MsaRepairTrigger trigger)
    {
        MsaRepairOptions options;
        options.trigger = trigger;
        return options;
    }

    // The bases of a row, with the gaps left out: the read it came from.
    static string msa1Ungap(const vector<AlignedBase>& v)
    {
        string s;
        for(const AlignedBase b: v) {
            if(not b.isGap()) {
                s.push_back(b.character());
            }
        }
        return s;
    }

    static string msa1ToString(const vector< pair<Base, uint64_t> >& consensus)
    {
        string s;
        for(const auto& [base, coverage]: consensus) {
            s.push_back(base.character());
        }
        return s;
    }
}



// Test the extended alphabet and its encoding.
void shasta2::testMsaRepairExtendedBase()
{
    // Check the symbol values are what the low-bits arithmetic in base() and
    // isPoly() assumes.
    for(uint64_t i=0; i<4; i++) {
        const Base base = Base::fromInteger(i);
        const ExtendedBase plain = ExtendedBase::fromBase(base);
        const ExtendedBase poly = ExtendedBase::polyFromBase(base);
        SHASTA2_ASSERT(plain.value == i);
        SHASTA2_ASSERT(poly.value == i + 4);
        SHASTA2_ASSERT(not plain.isPoly());
        SHASTA2_ASSERT(poly.isPoly());
        SHASTA2_ASSERT(plain.base() == base);
        SHASTA2_ASSERT(poly.base() == base);
        SHASTA2_ASSERT(plain.isValid());
        SHASTA2_ASSERT(poly.isValid());

        // The complement preserves the poly flag.
        SHASTA2_ASSERT(plain.complement().base() == base.complement());
        SHASTA2_ASSERT(poly.complement().base() == base.complement());
        SHASTA2_ASSERT(not plain.complement().isPoly());
        SHASTA2_ASSERT(poly.complement().isPoly());
        SHASTA2_ASSERT(plain.complement().complement() == plain);
        SHASTA2_ASSERT(poly.complement().complement() == poly);

        // AlignedExtendedBase agrees with ExtendedBase on the 8 real symbols.
        SHASTA2_ASSERT(AlignedExtendedBase(plain).value == plain.value);
        SHASTA2_ASSERT(AlignedExtendedBase(poly).value == poly.value);
        SHASTA2_ASSERT(not AlignedExtendedBase(plain).isGap());
        SHASTA2_ASSERT(AlignedExtendedBase(poly).isPoly());
        SHASTA2_ASSERT(ExtendedBase(AlignedExtendedBase(poly)) == poly);
    }

    // The gap lives in AlignedExtendedBase only. ExtendedBase cannot represent
    // one, which is the point: a sequence cannot contain a gap.
    {
        const AlignedExtendedBase gap = AlignedExtendedBase::gap();
        SHASTA2_ASSERT(gap.isGap());
        SHASTA2_ASSERT(not gap.isPoly());
        SHASTA2_ASSERT(gap.isValid());
        SHASTA2_ASSERT(gap.complement().isGap());
        SHASTA2_ASSERT(gap.character() == '-');

        // '-' is not a valid ExtendedBase.
        bool threw = false;
        try {
            ExtendedBase::fromCharacter('-');
        } catch(const std::exception&) {
            threw = true;
        }
        SHASTA2_ASSERT(threw);

        // The alphabet has exactly 8 symbols.
        for(uint64_t value=0; value<8; value++) {
            SHASTA2_ASSERT(ExtendedBase::fromInteger(value).isValid());
        }
        SHASTA2_ASSERT(not ExtendedBase::fromInteger(8UL).isValid());
    }

    // character() and fromCharacter() round trip over all symbols.
    for(uint64_t value=0; value<8; value++) {
        const ExtendedBase e = ExtendedBase::fromInteger(value);
        SHASTA2_ASSERT(ExtendedBase::fromCharacter(e.character()) == e);
    }
    for(uint64_t value=0; value<=AlignedExtendedBase::gapValue; value++) {
        const AlignedExtendedBase e = AlignedExtendedBase::fromInteger(value);
        SHASTA2_ASSERT(AlignedExtendedBase::fromCharacter(e.character()) == e);
    }

    // An invalid character throws.
    {
        bool threw = false;
        try {
            ExtendedBase::fromCharacter('X');
        } catch(const std::exception&) {
            threw = true;
        }
        SHASTA2_ASSERT(threw);
    }

    // toString and the fromString helpers round trip.
    // This matters because it is how the encoding reaches theseus and comes
    // back: msa1AlignExtended passes toString(encoded) in and parses the MSA
    // rows with vectorOfAlignedExtendedBasesFromString. If either mapping
    // folded 'a' into 'A', as AlignedBase::fromCharacter does, every poly
    // symbol would silently become a plain one.
    {
        SHASTA2_ASSERT(toString(vectorOfExtendedBasesFromString("ACGTacgt")) == "ACGTacgt");
        SHASTA2_ASSERT(toString(vectorOfAlignedExtendedBasesFromString("ACGTacgt-")) == "ACGTacgt-");
        SHASTA2_ASSERT(ExtendedBase::fromCharacter('a') != ExtendedBase::fromCharacter('A'));
        SHASTA2_ASSERT(ExtendedBase::fromCharacter('a').isPoly());
        SHASTA2_ASSERT(not ExtendedBase::fromCharacter('A').isPoly());
        SHASTA2_ASSERT(ExtendedBase::fromCharacter('a').base() ==
            ExtendedBase::fromCharacter('A').base());

        // AlignedBase, by contrast, would lose it. This documents exactly why
        // msa1AlignExtended must not use vectorOfAlignedBasesFromString.
        SHASTA2_ASSERT(AlignedBase::fromCharacter('a') == AlignedBase::fromCharacter('A'));
    }



    // Round trip all 19 sequences, at a range of thresholds, and check that they
    // all encode to the SAME symbol string. That is the property that makes the
    // homopolymer ambiguity structurally impossible: if every read has the same
    // run structure there is nothing to align, and no scoring choice can
    // misplace the G between the two A runs.
    for(uint64_t threshold=1; threshold<=12; threshold++) {
        ExtendedSequence firstEncoded;
        bool uniform = true;

        for(uint64_t i=0; i<msa1TestSequences.size(); i++) {
            const vector<Base> sequence = vectorOfBasesFromString(msa1TestSequences[i]);

            ExtendedSequence encoded;
            encodeExtended(sequence, threshold, encoded);

            // The encoding is lossless.
            vector<Base> decoded;
            decodeExtended(encoded, decoded);
            SHASTA2_ASSERT(decoded == sequence);

            // A plain symbol always stands for exactly one base, and a poly
            // symbol for more than the threshold.
            for(const auto& [e, runLength]: encoded) {
                if(e.isPoly()) {
                    SHASTA2_ASSERT(runLength > threshold);
                } else {
                    SHASTA2_ASSERT(runLength == 1);
                }
            }

            // A poly symbol stands for a whole run on its own, so it is never
            // adjacent to a symbol with the same base.
            for(uint64_t j=1; j<encoded.size(); j++) {
                if(encoded[j].first.isPoly() or encoded[j-1].first.isPoly()) {
                    SHASTA2_ASSERT(encoded[j].first.base() != encoded[j-1].first.base());
                }
            }

            if(i == 0) {
                firstEncoded = encoded;
            } else if(toString(encoded) != toString(firstEncoded)) {
                uniform = false;
            } else {
                SHASTA2_ASSERT(encoded.size() == firstEncoded.size());
            }
        }

        cout << "Threshold " << threshold << ": " <<
            (uniform ? "all 19 sequences encode to the same " : "encodings DIFFER, first has ") <<
            firstEncoded.size() << " symbols." << endl;

        // The shortest A run in this set is 6 bases and the others at the same
        // two positions are 7 or longer, so a threshold below 6 collapses all of
        // them and every sequence encodes identically. At a threshold of 6 the
        // run of 6 is no longer above the threshold, stays plain, and the
        // guarantee is lost. Asserted in both directions.
        if(threshold < 6) {
            SHASTA2_ASSERT(uniform);
        }
        if(threshold == 6) {
            SHASTA2_ASSERT(not uniform);
        }

        // The default must keep the encoding uniform here. This is the property
        // the whole alphabet exists to provide, so a change to the default that
        // loses it should not pass silently.
        if(threshold == defaultHomopolymerThreshold) {
            SHASTA2_ASSERT(uniform);
        }
    }



    // A sequence with no run longer than the threshold encodes to itself, one
    // plain symbol per base, so behavior on sequences without long homopolymers
    // is unchanged by this alphabet.
    {
        const vector<Base> sequence = vectorOfBasesFromString("ACGTACGTAACCGGTT");
        ExtendedSequence encoded;
        encodeExtended(sequence, defaultHomopolymerThreshold, encoded);
        SHASTA2_ASSERT(encoded.size() == sequence.size());
        for(uint64_t i=0; i<encoded.size(); i++) {
            SHASTA2_ASSERT(not encoded[i].first.isPoly());
            SHASTA2_ASSERT(encoded[i].first.base() == sequence[i]);
            SHASTA2_ASSERT(encoded[i].second == 1);
        }
    }



    // Behavior of a single run as its length crosses the threshold.
    //
    // A run of length threshold encodes to threshold plain symbols, and a run of
    // length threshold+1 to a single poly symbol. So the encoding is
    // discontinuous at the threshold: one extra base changes the symbol string
    // by threshold-1 symbols, and to an aligner a run that straddles the
    // threshold looks like a large indel rather than a small one. That is the
    // price of collapsing a whole run to one symbol, and it is why the threshold
    // has to be kept below the shortest run whose length is unreliable.
    {
        const uint64_t threshold = defaultHomopolymerThreshold;
        for(uint64_t runLength=1; runLength<=threshold+3; runLength++) {
            const vector<Base> sequence(runLength, Base::fromCharacter('A'));
            ExtendedSequence encoded;
            encodeExtended(sequence, threshold, encoded);

            if(runLength > threshold) {
                SHASTA2_ASSERT(encoded.size() == 1);
                SHASTA2_ASSERT(encoded[0].first.isPoly());
                SHASTA2_ASSERT(encoded[0].second == runLength);
            } else {
                SHASTA2_ASSERT(encoded.size() == runLength);
                for(const auto& [e, ignore]: encoded) {
                    SHASTA2_ASSERT(not e.isPoly());
                }
            }

            vector<Base> decoded;
            decodeExtended(encoded, decoded);
            SHASTA2_ASSERT(decoded == sequence);
        }

        // All runs longer than the threshold encode to the same single SYMBOL,
        // regardless of length. This is the property that removes long run
        // lengths from the alignment entirely: the symbol strings are equal, and
        // only the lengths carried alongside differ.
        ExtendedSequence shortRun;
        ExtendedSequence longRun;
        encodeExtended(vector<Base>(threshold + 1, Base::fromCharacter('A')),
            threshold, shortRun);
        encodeExtended(vector<Base>(1000, Base::fromCharacter('A')),
            threshold, longRun);
        SHASTA2_ASSERT(toString(shortRun) == toString(longRun));
        SHASTA2_ASSERT(shortRun.size() == 1);
        SHASTA2_ASSERT(shortRun[0].second != longRun[0].second);

        // The discontinuity at the threshold.
        ExtendedSequence atThreshold;
        encodeExtended(vector<Base>(threshold, Base::fromCharacter('A')),
            threshold, atThreshold);
        SHASTA2_ASSERT(atThreshold.size() == threshold);
    }



    // A threshold of 0 is legal and makes every run a poly symbol, so the
    // symbol string becomes the run length encoding of the sequence. It stays
    // lossless: the lengths move into the pairs.
    {
        const vector<Base> sequence = vectorOfBasesFromString("AACCCGTTTT");
        ExtendedSequence encoded;
        encodeExtended(sequence, 0, encoded);
        SHASTA2_ASSERT(toString(encoded) == "acgt");
        for(const auto& [symbol, runLength]: encoded) {
            SHASTA2_ASSERT(symbol.isPoly());
        }
        vector<Base> decoded;
        decodeExtended(encoded, decoded);
        SHASTA2_ASSERT(decoded == sequence);

        // Two sequences differing only in run lengths encode identically at 0,
        // which is what lets them stack with no aligner, but not at 1, where a
        // run of one stays plain and a run of two does not.
        ExtendedSequence a0, b0, a1, b1;
        const vector<Base> x = vectorOfBasesFromString("GGGGGTACAAGT");
        const vector<Base> y = vectorOfBasesFromString("GGGGGTACAGTT");
        encodeExtended(x, 0, a0);
        encodeExtended(y, 0, b0);
        encodeExtended(x, 1, a1);
        encodeExtended(y, 1, b1);
        SHASTA2_ASSERT(toString(a0) == toString(b0));
        SHASTA2_ASSERT(toString(a1) != toString(b1));
    }



    // Degenerate inputs.
    {
        ExtendedSequence encoded;
        vector<Base> decoded;

        // Empty sequence.
        encodeExtended(vector<Base>(), defaultHomopolymerThreshold, encoded);
        SHASTA2_ASSERT(encoded.empty());
        decodeExtended(encoded, decoded);
        SHASTA2_ASSERT(decoded.empty());

        // Single base.
        const vector<Base> single = vectorOfBasesFromString("G");
        encodeExtended(single, defaultHomopolymerThreshold, encoded);
        SHASTA2_ASSERT(encoded.size() == 1);
        SHASTA2_ASSERT(not encoded[0].first.isPoly());
        decodeExtended(encoded, decoded);
        SHASTA2_ASSERT(decoded == single);

        // One run spanning the whole sequence: a single poly symbol.
        const vector<Base> oneRun(100, Base::fromCharacter('T'));
        encodeExtended(oneRun, defaultHomopolymerThreshold, encoded);
        SHASTA2_ASSERT(encoded.size() == 1);
        SHASTA2_ASSERT(encoded[0].first.isPoly());
        SHASTA2_ASSERT(encoded[0].second == oneRun.size());
        decodeExtended(encoded, decoded);
        SHASTA2_ASSERT(decoded == oneRun);
    }



    // Round trip random sequences over a range of thresholds.
    // A simple deterministic generator, so a failure is reproducible.
    {
        uint64_t random = 38495721;
        const auto next = [&random]() {
            random = random * 6364136223846793005ULL + 1442695040888963407ULL;
            return uint64_t(random >> 33);
        };

        for(uint64_t threshold=1; threshold<=8; threshold++) {
            for(uint64_t trial=0; trial<200; trial++) {
                vector<Base> sequence;
                const uint64_t runCount = next() % 40;
                Base previousBase = Base::fromInteger(uint64_t(next() % 4));
                for(uint64_t i=0; i<runCount; i++) {
                    // Consecutive runs must have different bases, so the run
                    // structure is what we intend it to be.
                    Base base = previousBase;
                    if(i > 0) {
                        base = Base::fromInteger(uint64_t((previousBase.value + 1 + next() % 3) % 4));
                        SHASTA2_ASSERT(base != previousBase);
                    }
                    const uint64_t runLength = 1 + next() % 20;
                    for(uint64_t j=0; j<runLength; j++) {
                        sequence.push_back(base);
                    }
                    previousBase = base;
                }

                ExtendedSequence encoded;
                encodeExtended(sequence, threshold, encoded);

                vector<Base> decoded;
                decodeExtended(encoded, decoded);
                SHASTA2_ASSERT(decoded == sequence);

                // The encoding never grows the sequence.
                SHASTA2_ASSERT(encoded.size() <= max(sequence.size(), 1UL));
            }
        }
    }

    cout << "testMsa1ExtendedBase passed." << endl;
}


// Test the consensus computed from an alignment over the extended alphabet.
//
// Theseus is not exercised here: extendedConsensus, attachRunLengths and
// expandExtendedAlignment are pure functions, so the alignment is supplied
// directly. That also lets the real, defective abpoa alignment be used as input.
void shasta2::testMsaRepairConsensus()
{
    // Confirm the input really is defective. If this stops being true the test
    // has lost its point.
    {
        const uint64_t impureColumnCount = msaRepairImpureColumnCount(msa1BadAlignmentRows);
        cout << "The input alignment has " << impureColumnCount <<
            " columns containing more than one base." << endl;
        SHASTA2_ASSERT(impureColumnCount == 4);
    }

    // Recover the sequences from that alignment and encode them. They all encode
    // to the same symbol string, so they line up without any further work; this
    // stands in for what theseus would return. 11 distinct base sequences
    // collapse to 1 distinct encoding.
    const uint64_t threshold = defaultHomopolymerThreshold;
    vector<AlignedExtendedSequence> alignment;
    vector<uint64_t> weights;
    for(const string& row: msa1BadAlignmentRows) {
        string ungapped;
        for(const char c: row) {
            if(c != '-') {
                ungapped.push_back(c);
            }
        }
        ExtendedSequence encoded;
        encodeExtended(vectorOfBasesFromString(ungapped), threshold, encoded);

        // No gaps needed: every row has the same symbol string.
        AlignedExtendedSequence alignedRow;
        for(const auto& [e, runLength]: encoded) {
            alignedRow.push_back(make_pair(AlignedExtendedBase(e), runLength));
        }
        alignment.push_back(alignedRow);
        weights.push_back(1);
    }
    for(const AlignedExtendedSequence& row: alignment) {
        SHASTA2_ASSERT(toString(row) == toString(alignment.front()));
    }
    cout << "All " << alignment.size() << " sequences encode to the same " <<
        alignment.front().size() << " symbols." << endl;



    // Vote, using the default estimator.
    vector< pair<Base, uint64_t> > consensus;
    AlignedExtendedSequence alignedConsensus;
    extendedConsensus(alignment, weights, RunLengthEstimator::Mode,
        consensus, alignedConsensus);
    const string consensusString = msa1ToString(consensus);

    // The expected consensus. The two variable A runs come out at 12 and 11.
    //
    // These are the KNOWN TRUE run lengths at this locus, so this is a check
    // against ground truth and not merely against what some other tool produces.
    // The median gives 11 and 11 here and is wrong: 10 of the 19 reads report the
    // first run as shorter than 12, so the median lands one base low even though
    // 12 is the single most common length. See RunLengthEstimator.
    const string expectedConsensus =
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTT"
        "AAACTATAAAGTAAATTCCTCCCATAGTT"
        "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAGTT"
        "AAACTATAAAGTAAATTCCTCCCATAGTT";
    cout << "Consensus length " << consensusString.size() << "." << endl;
    SHASTA2_ASSERT(consensusString == expectedConsensus);

    // The consensus is the aligned consensus expanded.
    {
        string fromAligned;
        for(const auto& [e, runLength]: alignedConsensus) {
            if(e.isGap()) {
                SHASTA2_ASSERT(runLength == 0);
                continue;
            }
            for(uint64_t k=0; k<runLength; k++) {
                fromAligned.push_back(e.base().character());
            }
        }
        SHASTA2_ASSERT(fromAligned == consensusString);
    }

    // Coverage is the weight supporting the chosen length, and can never exceed
    // the total weight.
    {
        uint64_t totalWeight = 0;
        for(const uint64_t weight: weights) {
            totalWeight += weight;
        }
        for(const auto& [base, coverage]: consensus) {
            SHASTA2_ASSERT(coverage > 0);
            SHASTA2_ASSERT(coverage <= totalWeight);
        }
    }



    // The median gives a different, and here incorrect, answer at the first run:
    // 11 rather than the true 12, in each of the two tandem copies. Asserted so
    // the two estimators cannot silently converge, and so the direction of the
    // difference is recorded.
    {
        vector< pair<Base, uint64_t> > medianConsensus;
        AlignedExtendedSequence ignore;
        extendedConsensus(alignment, weights, RunLengthEstimator::Median,
            medianConsensus, ignore);
        const string medianString = msa1ToString(medianConsensus);
        cout << "Mode gives length " << consensusString.size() <<
            " (correct), median gives length " << medianString.size() << "." << endl;
        SHASTA2_ASSERT(medianString != consensusString);
        SHASTA2_ASSERT(medianString.size() == consensusString.size() - 2);
    }



    // Expanding must preserve every read and give every row the same length as
    // the aligned consensus.
    {
        vector< vector<AlignedBase> > expandedAlignment;
        vector<AlignedBase> expandedAlignedConsensus;
        expandExtendedAlignment(alignment, alignedConsensus,
            expandedAlignment, expandedAlignedConsensus);

        SHASTA2_ASSERT(expandedAlignment.size() == msa1BadAlignmentRows.size());
        vector<string> expandedRows;
        for(uint64_t i=0; i<expandedAlignment.size(); i++) {
            SHASTA2_ASSERT(expandedAlignment[i].size() == expandedAlignedConsensus.size());

            string ungapped;
            string full;
            for(const AlignedBase b: expandedAlignment[i]) {
                full.push_back(b.character());
                if(not b.isGap()) {
                    ungapped.push_back(b.character());
                }
            }
            expandedRows.push_back(full);

            string expected;
            for(const char c: msa1BadAlignmentRows[i]) {
                if(c != '-') {
                    expected.push_back(c);
                }
            }
            SHASTA2_ASSERT(ungapped == expected);
        }

        // No column of the expanded alignment holds two different bases. The
        // input had 4. This is structural: a poly symbol and the base next to it
        // are different symbols and get different columns, so they cannot be
        // confused.
        const uint64_t impureColumnCount = msaRepairImpureColumnCount(expandedRows);
        cout << "The expanded alignment has " << impureColumnCount <<
            " columns containing more than one base, in " <<
            expandedAlignedConsensus.size() << " columns." << endl;
        SHASTA2_ASSERT(impureColumnCount == 0);
    }



    // Weights must matter, for both estimators. Give the two reads whose first
    // run is 10 and whose second is 12 overwhelming weight, and the consensus
    // must follow them rather than the other seventeen reads.
    {
        vector<uint64_t> skewedWeights = weights;
        skewedWeights[9] = 100;
        skewedWeights[10] = 100;

        for(const RunLengthEstimator estimator:
            {RunLengthEstimator::Median, RunLengthEstimator::Mode}) {
            vector< pair<Base, uint64_t> > skewedConsensus;
            AlignedExtendedSequence ignore;
            extendedConsensus(alignment, skewedWeights, estimator,
                skewedConsensus, ignore);
            const string expected =
                "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC" + string(10, 'A') + "G" +
                string(12, 'A') + "GTTAAACTATAAAGTAAATTCCTCCCATAGTT" +
                "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTC" + string(10, 'A') + "G" +
                string(12, 'A') + "GTTAAACTATAAAGTAAATTCCTCCCATAGTT";
            SHASTA2_ASSERT(msa1ToString(skewedConsensus) == expected);
        }
    }



    // A real substitution must survive: the majority allele is kept, not
    // absorbed or dropped, and the consensus follows when the split moves.
    {
        vector<AlignedExtendedSequence> snpAlignment = alignment;

        uint64_t column = invalid<uint64_t>;
        for(uint64_t j=0; j<snpAlignment.front().size(); j++) {
            if(not snpAlignment.front()[j].first.isPoly()) {
                column = j;
                break;
            }
        }
        SHASTA2_ASSERT(column != invalid<uint64_t>);
        const AlignedExtendedBase original = snpAlignment.front()[column].first;
        const AlignedExtendedBase mutated =
            AlignedExtendedBase(ExtendedBase::fromBase(original.base().complement()));
        SHASTA2_ASSERT(mutated != original);

        for(uint64_t i=0; i<8; i++) {
            snpAlignment[i][column].first = mutated;
        }
        vector< pair<Base, uint64_t> > snpConsensus;
        AlignedExtendedSequence snpAlignedConsensus;
        extendedConsensus(snpAlignment, weights, RunLengthEstimator::Mode,
            snpConsensus, snpAlignedConsensus);
        SHASTA2_ASSERT(snpAlignedConsensus[column].first == original);
        SHASTA2_ASSERT(snpAlignedConsensus[column].second == 1);

        uint64_t positionInConsensus = 0;
        for(uint64_t j=0; j<column; j++) {
            positionInConsensus += snpAlignedConsensus[j].second;
        }
        SHASTA2_ASSERT(snpConsensus[positionInConsensus].second == 11);
        cout << "With a substitution in 8 of 19 rows, the majority allele is kept "
            "with coverage " << snpConsensus[positionInConsensus].second << "." << endl;

        // Flip it in 12 rows and the consensus must follow.
        for(uint64_t i=0; i<12; i++) {
            snpAlignment[i][column].first = mutated;
        }
        extendedConsensus(snpAlignment, weights, RunLengthEstimator::Mode,
            snpConsensus, snpAlignedConsensus);
        SHASTA2_ASSERT(snpAlignedConsensus[column].first == mutated);
    }



    // A column that is a gap in a majority of rows becomes a gap and contributes
    // nothing to the consensus.
    //
    // The extra column goes in the middle, not at either end. At an end it would
    // instead test the padding rule: a row whose first symbol is the inserted C
    // is the only one covering that column, and the rows that stop short of it
    // abstain rather than voting for a gap, so the C would win. That rule is
    // tested separately below.
    {
        vector<AlignedExtendedSequence> gapAlignment = alignment;
        const uint64_t column = gapAlignment.front().size() / 2;
        for(AlignedExtendedSequence& row: gapAlignment) {
            row.insert(row.begin() + int64_t(column),
                make_pair(AlignedExtendedBase::gap(), 0UL));
        }
        gapAlignment[0][column] = make_pair(
            AlignedExtendedBase::fromCharacter('C'), 1UL);

        vector< pair<Base, uint64_t> > gapConsensus;
        AlignedExtendedSequence gapAlignedConsensus;
        extendedConsensus(gapAlignment, weights, RunLengthEstimator::Mode,
            gapConsensus, gapAlignedConsensus);
        SHASTA2_ASSERT(gapAlignedConsensus[column].first.isGap());
        SHASTA2_ASSERT(gapAlignedConsensus[column].second == 0);
        SHASTA2_ASSERT(msa1ToString(gapConsensus) == consensusString);
    }



    // Padding at the ends of a row is not a vote for a gap.
    //
    // Two rows are cut back to the second half of the alignment, as a read
    // constrained on the right only comes out of theseus. The part they do not
    // cover is padded with the same gap symbol used for a deletion, and if that
    // padding voted it would carry the column and delete the first half of the
    // consensus. Instead those rows abstain there and the consensus is unchanged.
    {
        vector<AlignedExtendedSequence> paddedAlignment = alignment;
        const uint64_t column = paddedAlignment.front().size() / 2;
        for(uint64_t i=0; i<paddedAlignment.size(); i++) {
            if(i == 0) {
                continue;
            }
            for(uint64_t j=0; j<column; j++) {
                paddedAlignment[i][j] = make_pair(AlignedExtendedBase::gap(), 0UL);
            }
        }

        vector< pair<Base, uint64_t> > paddedConsensus;
        AlignedExtendedSequence paddedAlignedConsensus;
        extendedConsensus(paddedAlignment, weights, RunLengthEstimator::Mode,
            paddedConsensus, paddedAlignedConsensus);
        SHASTA2_ASSERT(msa1ToString(paddedConsensus) == consensusString);

        // Told explicitly that every row is anchored on both sides, the same
        // gaps are deletions and the first half does go away. The two readings
        // are genuinely different, which is why the caller has to say which one
        // applies rather than have it guessed.
        const vector< pair<uint64_t, uint64_t> > spans(
            paddedAlignment.size(),
            make_pair(0UL, paddedAlignment.front().size()));
        extendedConsensus(paddedAlignment, weights, RunLengthEstimator::Mode,
            spans, paddedConsensus, paddedAlignedConsensus);
        for(uint64_t j=0; j<column; j++) {
            SHASTA2_ASSERT(paddedAlignedConsensus[j].first.isGap());
        }
    }



    // No regression when nothing is long enough to become a poly symbol: the
    // result must equal plain weighted column majority.
    {
        const vector<string> rows = {"ACGTACGTAC", "ACGTACGTAC", "ACGTTCGTAC"};
        vector<AlignedExtendedSequence> a;
        for(const string& row: rows) {
            AlignedExtendedSequence r;
            for(const char c: row) {
                r.push_back(make_pair(AlignedExtendedBase::fromCharacter(c), 1UL));
                SHASTA2_ASSERT(not r.back().first.isPoly());
            }
            a.push_back(r);
        }
        vector< pair<Base, uint64_t> > c;
        AlignedExtendedSequence ac;
        extendedConsensus(a, {1, 1, 1}, RunLengthEstimator::Mode, c, ac);
        SHASTA2_ASSERT(msa1ToString(c) == "ACGTACGTAC");
        SHASTA2_ASSERT(c[4].second == 2);
        SHASTA2_ASSERT(c[0].second == 3);
    }



    // attachRunLengths puts the lengths back on the columns of a row.
    {
        ExtendedSequence encoded;
        encodeExtended(vectorOfBasesFromString(
            "CTC" + string(12, 'A') + "G" + string(11, 'A')), 4, encoded);
        SHASTA2_ASSERT(toString(encoded) == "CTCaGa");

        const vector<AlignedExtendedBase> row =
            vectorOfAlignedExtendedBasesFromString("C-TC-aG-a");
        AlignedExtendedSequence attached;
        attachRunLengths(row, encoded, Anchoring::BothSides, attached);
        const vector<uint64_t> expected = {1, 0, 1, 1, 0, 12, 1, 0, 11};
        SHASTA2_ASSERT(attached.size() == expected.size());
        for(uint64_t i=0; i<expected.size(); i++) {
            SHASTA2_ASSERT(attached[i].second == expected[i]);
        }
        SHASTA2_ASSERT(toString(attached) == "C-TC-aG-a");
    }



    // A row can cover only part of its encoding, because theseus trims the
    // overhang of a sequence fixed on one side only.
    {
        ExtendedSequence encoded;
        encodeExtended(vectorOfBasesFromString(
            "CTC" + string(12, 'A') + "G" + string(11, 'A') + string(7, 'T')), 4, encoded);
        SHASTA2_ASSERT(toString(encoded) == "CTCaGat");

        // Left fixed: the trailing poly symbol is trimmed.
        AlignedExtendedSequence attached;
        attachRunLengths(vectorOfAlignedExtendedBasesFromString("CTC-aGa"),
            encoded, Anchoring::LeftOnly, attached);
        const vector<uint64_t> expected = {1, 1, 1, 0, 12, 1, 11};
        for(uint64_t i=0; i<expected.size(); i++) {
            SHASTA2_ASSERT(attached[i].second == expected[i]);
        }
    }
    {
        // Right fixed: the leading poly symbol is trimmed, so the lengths must
        // come from offset 1, not 0.
        ExtendedSequence encoded;
        encodeExtended(vectorOfBasesFromString(
            string(9, 'G') + "CTC" + string(12, 'A') + "G" + string(11, 'A')), 4, encoded);
        SHASTA2_ASSERT(toString(encoded) == "gCTCaGa");

        AlignedExtendedSequence attached;
        attachRunLengths(vectorOfAlignedExtendedBasesFromString("CTCaG-a"),
            encoded, Anchoring::RightOnly, attached);
        const vector<uint64_t> expected = {1, 1, 1, 12, 1, 0, 11};
        for(uint64_t i=0; i<expected.size(); i++) {
            SHASTA2_ASSERT(attached[i].second == expected[i]);
            // The trimmed poly symbol's length of 9 must not appear anywhere.
            SHASTA2_ASSERT(attached[i].second != 9);
        }
    }



    // The window must come from the anchoring, not from searching the encoding
    // for the row's symbols.
    //
    // Here the encoding is a tandem repeat, A a C A a C, whose two poly runs have
    // different lengths, 10 and 20. The sequence is fixed on the right, so
    // theseus trimmed the first copy and the row is the second one, whose run is
    // 20. Searching for the row's symbols finds the FIRST copy and would pair the
    // row with 10 instead. Tandem repeats are exactly what this code is for, and
    // the 19 read locus it was developed on is itself a 2x tandem duplication.
    {
        ExtendedSequence encoded;
        encodeExtended(vectorOfBasesFromString(
            "A" + string(10, 'G') + "C" + "A" + string(20, 'G') + "C"), 4, encoded);
        SHASTA2_ASSERT(toString(encoded) == "AgCAgC");
        SHASTA2_ASSERT(encoded[1].second == 10);
        SHASTA2_ASSERT(encoded[4].second == 20);

        const vector<AlignedExtendedBase> row =
            vectorOfAlignedExtendedBasesFromString("AgC");

        AlignedExtendedSequence attached;
        attachRunLengths(row, encoded, Anchoring::RightOnly, attached);
        SHASTA2_ASSERT(attached[1].second == 20);

        // And the mirror image: fixed on the left, so the row is the first copy.
        attachRunLengths(row, encoded, Anchoring::LeftOnly, attached);
        SHASTA2_ASSERT(attached[1].second == 10);
    }



    // Degenerate inputs.
    {
        vector< pair<Base, uint64_t> > c;
        AlignedExtendedSequence ac;

        extendedConsensus({}, {}, RunLengthEstimator::Mode, c, ac);
        SHASTA2_ASSERT(c.empty());
        SHASTA2_ASSERT(ac.empty());

        // A single sequence: the consensus is that sequence.
        ExtendedSequence encoded;
        encodeExtended(vectorOfBasesFromString(
            "CTC" + string(12, 'A') + "G" + string(11, 'A')), 4, encoded);
        AlignedExtendedSequence one;
        for(const auto& [e, runLength]: encoded) {
            one.push_back(make_pair(AlignedExtendedBase(e), runLength));
        }
        extendedConsensus({one}, {7}, RunLengthEstimator::Mode, c, ac);
        SHASTA2_ASSERT(msa1ToString(c) == "CTC" + string(12, 'A') + "G" + string(11, 'A'));
        for(const auto& [base, coverage]: c) {
            SHASTA2_ASSERT(coverage == 7);
        }

        // Every row a gap at every column.
        const AlignedExtendedSequence allGapRow = {
            make_pair(AlignedExtendedBase::gap(), 0UL),
            make_pair(AlignedExtendedBase::gap(), 0UL)};
        extendedConsensus({allGapRow, allGapRow}, {1, 1},
            RunLengthEstimator::Mode, c, ac);
        SHASTA2_ASSERT(c.empty());
        SHASTA2_ASSERT(ac.size() == 2);
        SHASTA2_ASSERT(ac[0].first.isGap() and ac[1].first.isGap());
    }

    cout << "testMsa1Consensus passed." << endl;
}


// Test the local repair.
//
// The property that matters most here is the first one: on input with no bad
// region, nothing is modified at all. The whole point of repairing locally is
// that sequence accuracy in the rest of the assembly is left alone, so a test
// that the repair is a no-op where it should be is more important than a test
// that it works where it should.
void shasta2::testMsaRepairRepair()
{
    const uint64_t threshold = defaultHomopolymerThreshold;

    // These tests were written against the pattern trigger and assert on the
    // regions it finds, so they pin that trigger explicitly rather than
    // following the default.
    const MsaRepairTrigger trigger = MsaRepairTrigger::PatternOnly;

    // Build the real, defective abpoa alignment.
    vector< vector<AlignedBase> > alignment;
    vector<string> reads;
    for(const string& row: msa1BadAlignmentRows) {
        alignment.push_back(vectorOfAlignedBasesFromString(row));
        string ungapped;
        for(const char c: row) {
            if(c != '-') {
                ungapped.push_back(c);
            }
        }
        reads.push_back(ungapped);
    }
    const vector<uint64_t> weights(alignment.size(), 1);

    // Its column-wise consensus, which is what abpoa reports.
    vector<AlignedBase> alignedConsensus;
    vector< pair<Base, uint64_t> > consensus;
    msa1ColumnConsensus(alignment, weights, alignedConsensus, consensus);



    // The prescreen must fire on these reads, which contain two long A runs
    // separated by a single G.
    {
        vector< vector<Base> > sequences;
        for(const string& r: reads) {
            sequences.push_back(vectorOfBasesFromString(r));
        }
        SHASTA2_ASSERT(msaRepairTriggerPresent(sequences, MsaRepairTrigger::PatternOnly, threshold));

        // And it must not fire on ordinary sequence with no long runs. This is
        // the common case in real reads, where 97% of runs are 4 bases or
        // shorter, and it is what keeps the repair from running at all.
        const vector< vector<Base> > ordinary = {
            vectorOfBasesFromString("ACGTACGTAACCGGTTACGTACGT"),
            vectorOfBasesFromString("ACGTACGTAACCGGTTACGTACGT")};
        SHASTA2_ASSERT(not msaRepairTriggerPresent(ordinary, MsaRepairTrigger::PatternOnly, threshold));

        // Nor on a single long run with no separating base.
        const vector< vector<Base> > oneRun = {
            vectorOfBasesFromString("ACGT" + string(20, 'A') + "CGTA")};
        SHASTA2_ASSERT(not msaRepairTriggerPresent(oneRun, MsaRepairTrigger::PatternOnly, threshold));

        // Nor when the two long runs are separated by more than one base.
        const vector< vector<Base> > twoBases = {
            vectorOfBasesFromString(string(12, 'A') + "GT" + string(11, 'A'))};
        SHASTA2_ASSERT(not msaRepairTriggerPresent(twoBases, MsaRepairTrigger::PatternOnly, threshold));

        // But yes when they are separated by exactly one.
        const vector< vector<Base> > oneBase = {
            vectorOfBasesFromString(string(12, 'A') + "G" + string(11, 'A'))};
        SHASTA2_ASSERT(msaRepairTriggerPresent(oneBase, MsaRepairTrigger::PatternOnly, threshold));

        // The exact boundaries of the pattern, on the single sequence primitive.
        // A run must be strictly longer than the threshold on BOTH sides, and
        // the separator must be exactly one base.
        const auto pattern = [&](uint64_t left, uint64_t middle, uint64_t right) {
            return msaRepairPatternPresent(vectorOfBasesFromString(
                string(left, 'A') + string(middle, 'G') + string(right, 'A')),
                threshold);
        };
        SHASTA2_ASSERT(pattern(threshold + 1, 1, threshold + 1));
        SHASTA2_ASSERT(not pattern(threshold, 1, threshold + 1));  // left too short
        SHASTA2_ASSERT(not pattern(threshold + 1, 1, threshold));  // right too short
        SHASTA2_ASSERT(not pattern(threshold + 1, 2, threshold + 1)); // separator too long
        SHASTA2_ASSERT(pattern(100, 1, 100));

        // The pattern is found wherever it sits, including at the very start and
        // the very end, which a scan that keeps only two previous run lengths
        // could get wrong.
        SHASTA2_ASSERT(msaRepairPatternPresent(vectorOfBasesFromString(
            string(12, 'A') + "G" + string(11, 'A') + "CGTACGT"), threshold));
        SHASTA2_ASSERT(msaRepairPatternPresent(vectorOfBasesFromString(
            "CGTACGT" + string(12, 'A') + "G" + string(11, 'A')), threshold));

        // Degenerate inputs.
        SHASTA2_ASSERT(not msaRepairPatternPresent(vector<Base>(), threshold));
        SHASTA2_ASSERT(not msaRepairPatternPresent(
            vectorOfBasesFromString("A"), threshold));
        SHASTA2_ASSERT(not msaRepairPatternPresent(
            vector<Base>(1000, Base::fromCharacter('A')), threshold));
    }



    // THE IMPORTANT TEST. On an alignment with no bad region, nothing changes.
    {
        const vector<string> cleanRows = {
            "ACGTACGTAACCGGTTACGTACGT",
            "ACGTACGTAACCGGTTACGTACGT",
            "ACGTACGTAACCTGTTACGTACGT"};
        vector< vector<AlignedBase> > a = msaRepairAlignmentFromStrings(cleanRows);
        vector<AlignedBase> ac;
        vector< pair<Base, uint64_t> > c;
        const vector<uint64_t> w(a.size(), 1);
        msa1ColumnConsensus(a, w, ac, c);

        // Note this alignment DOES have a discordant column, at the T for C
        // substitution, so the region finder has something to cluster. It is the
        // pattern confirmation that must reject it.
        const auto aBefore = a;
        const auto acBefore = ac;
        const auto cBefore = c;

        const uint64_t repaired = msaRepair(a, ac, c, w, {}, msa1TestOptions(trigger));
        SHASTA2_ASSERT(repaired == 0);
        SHASTA2_ASSERT(a == aBefore);
        SHASTA2_ASSERT(ac == acBefore);
        SHASTA2_ASSERT(c == cBefore);
        cout << "Repair left a clean alignment untouched." << endl;
    }



    // On the real defective alignment, the repair fires and fixes it.
    {
        vector< vector<AlignedBase> > a = alignment;
        vector<AlignedBase> ac = alignedConsensus;
        vector< pair<Base, uint64_t> > c = consensus;

        // Before: 4 columns hold more than one base.
        SHASTA2_ASSERT(msaRepairImpureColumnCount(a) == 4);

        vector<MsaRepairRegion> regions;
        msaRepairFindBadRegions(a, ac, trigger, threshold, 10, 20, {}, regions);
        cout << "Found " << regions.size() << " bad region(s) in the real alignment." << endl;
        SHASTA2_ASSERT(not regions.empty());

        const uint64_t repaired = msaRepair(a, ac, c, weights, {}, msa1TestOptions(trigger));
        cout << "Repaired " << repaired << " region(s)." << endl;
        SHASTA2_ASSERT(repaired > 0);

        // After: no column holds more than one base.
        {
            const uint64_t impure = msaRepairImpureColumnCount(a);
            cout << "After repair the alignment has " << impure <<
                " columns containing more than one base." << endl;
            SHASTA2_ASSERT(impure == 0);
        }

        // Every row still ungaps to its original read: the repair rearranges
        // columns, it does not alter the reads.
        for(uint64_t i=0; i<a.size(); i++) {
            SHASTA2_ASSERT(msa1Ungap(a[i]) == reads[i]);
        }

        // All rows still have the length of the aligned consensus.
        for(const auto& row: a) {
            SHASTA2_ASSERT(row.size() == ac.size());
        }

        // consensus is alignedConsensus with the gaps removed.
        SHASTA2_ASSERT(msa1Ungap(ac) == msa1ToString(c));

        // The known true run lengths at this locus are 12 and 11 (see the 19
        // sequences above). The production default, MedianMarginGated, gets
        // the first run exactly right (median 11, cumulative support 10/19 =
        // 52.6%, below the 60% gate, nudged up to the true 12) but overcorrects
        // the second (median 11, cumulative support 11/19 = 57.9%, also below
        // the 60% gate, nudged up to 12 though 11 was already correct). This is
        // the documented trade-off in RunLengthEstimator, not a bug:
        // MedianMarginGated is the default because it minimizes total error
        // over a whole assembly (see msaRepair.hpp), not because it is exact on
        // every single locus, and this is the one locus this file has always
        // checked by hand where that trade-off is visible.
        const string s = msa1ToString(c);
        const string expected =
            "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAAGTT"
            "AAACTATAAAGTAAATTCCTCCCATAGTT"
            "TCCAGCCTGGGTGACAGAGCGAGACCCCAACTCAAAAAAAAAAAAGAAAAAAAAAAAAGTT"
            "AAACTATAAAGTAAATTCCTCCCATAGTT";
        cout << "Consensus after repair has length " << s.size() << "." << endl;
        SHASTA2_ASSERT(s == expected);
    }



    // Only the bad regions change. Everything outside them is untouched, which
    // is the whole reason for repairing locally.
    {
        vector< vector<AlignedBase> > a = alignment;
        vector<AlignedBase> ac = alignedConsensus;
        vector< pair<Base, uint64_t> > c = consensus;

        vector<MsaRepairRegion> regions;
        msaRepairFindBadRegions(a, ac, trigger, threshold, 10, 20, {}, regions);
        SHASTA2_ASSERT(not regions.empty());

        // Record the columns before the first region and after the last one.
        const uint64_t firstBegin = regions.front().begin;
        const uint64_t lastEnd = regions.back().end;
        const uint64_t tailLength = ac.size() - lastEnd;
        vector<string> prefixBefore;
        vector<string> suffixBefore;
        for(const auto& row: a) {
            prefixBefore.push_back(msaRepairToString(row).substr(0, firstBegin));
            suffixBefore.push_back(msaRepairToString(row).substr(lastEnd));
        }

        msaRepair(a, ac, c, weights, {}, msa1TestOptions(trigger));

        // The prefix is at the same columns and unchanged.
        for(uint64_t i=0; i<a.size(); i++) {
            SHASTA2_ASSERT(msaRepairToString(a[i]).substr(0, firstBegin) == prefixBefore[i]);
        }
        // The suffix is unchanged, though it has moved if the repair changed the
        // number of columns.
        for(uint64_t i=0; i<a.size(); i++) {
            const string row = msaRepairToString(a[i]);
            SHASTA2_ASSERT(row.substr(row.size() - tailLength) == suffixBefore[i]);
        }
        cout << "Columns outside the repaired regions are unchanged." << endl;
    }



    // The consensus outside a repaired region is untouched, coverage included.
    //
    // This is the reason the consensus is spliced rather than recomputed from
    // the repaired alignment. Recomputing it rewrote the coverage of every
    // position in the sequence, because the coverage abpoa reports is computed
    // its own way and is not what a recount of the rows gives. That went
    // unnoticed at first because the test built the consensus with the same
    // formula the recomputation used, so the two agreed by construction. Here
    // the coverage is deliberately set to a value no recount could produce.
    {
        vector< vector<AlignedBase> > a = alignment;
        vector<AlignedBase> ac = alignedConsensus;
        vector< pair<Base, uint64_t> > c = consensus;
        const uint64_t sentinel = 777;
        for(auto& [base, coverage]: c) {
            coverage = sentinel;
        }
        const auto before = c;

        // Which consensus bases lie inside a region, worked out before the
        // repair moves anything.
        vector<MsaRepairRegion> regions;
        msaRepairFindBadRegions(a, ac, trigger, threshold, 10, 20, {}, regions);
        SHASTA2_ASSERT(not regions.empty());
        vector<bool> isInsideRegion(c.size(), false);
        {
            uint64_t base = 0;
            for(uint64_t j=0; j<ac.size(); j++) {
                if(ac[j].isGap()) {
                    continue;
                }
                for(const MsaRepairRegion& region: regions) {
                    if((j >= region.begin) and (j < region.end)) {
                        isInsideRegion[base] = true;
                    }
                }
                ++base;
            }
        }
        uint64_t insideCount = 0;
        for(const bool b: isInsideRegion) {
            if(b) {
                ++insideCount;
            }
        }
        SHASTA2_ASSERT(insideCount > 0);
        SHASTA2_ASSERT(insideCount < c.size());

        SHASTA2_ASSERT(msaRepair(a, ac, c, weights, {}, msa1TestOptions(trigger)) > 0);

        // The bases before the first repaired region and after the last must be
        // unchanged, coverage and all.
        uint64_t firstInside = 0;
        while((firstInside < isInsideRegion.size()) and (not isInsideRegion[firstInside])) {
            ++firstInside;
        }
        uint64_t lastInside = 0;
        for(uint64_t i=0; i<isInsideRegion.size(); i++) {
            if(isInsideRegion[i]) {
                lastInside = i;
            }
        }
        const uint64_t tailLength = before.size() - lastInside - 1;

        for(uint64_t i=0; i<firstInside; i++) {
            SHASTA2_ASSERT(c[i] == before[i]);
        }
        for(uint64_t i=0; i<tailLength; i++) {
            SHASTA2_ASSERT(c[c.size() - 1 - i] == before[before.size() - 1 - i]);
        }

        // And the coverage inside the regions really was replaced, or the test
        // above would pass for the wrong reason.
        //
        // Compared against before.size() - insideCount, both measured on the
        // PRE-repair consensus, not against the post-repair c.size(): a repair
        // can change how many consensus BASES a region covers (this is exactly
        // what happens when a run length changes), so insideCount and c.size()
        // taken from different sides of the repair are not the same base
        // count and are not interchangeable.
        uint64_t stillSentinel = 0;
        for(const auto& [base, coverage]: c) {
            if(coverage == sentinel) {
                ++stillSentinel;
            }
        }
        SHASTA2_ASSERT(stillSentinel == before.size() - insideCount);
        cout << "Consensus outside the repaired regions is unchanged: " <<
            stillSentinel << " of " << c.size() <<
            " positions keep their original coverage." << endl;
    }



    // A row with no base in a region votes for a deletion there, and is outvoted
    // when it is in the minority.
    //
    // Such a row has nothing to align, so it is left out of the realignment and
    // comes back as the gaps it already was. That is not the same as leaving it
    // out of the vote: it is fixed on both sides, so it covers this region and
    // says the region is deleted. Here two rows say otherwise and the sequence
    // stays.
    {
        const string sequence = "AAAAAAAAAAAAGAAAAAAAAAAA";
        const vector<string> rows = {sequence, sequence, string(24, '-')};
        vector< vector<AlignedBase> > a = msaRepairAlignmentFromStrings(rows);
        vector<AlignedBase> ac = vectorOfAlignedBasesFromString(sequence);
        vector< pair<Base, uint64_t> > c;
        for(const char ch: sequence) {
            c.push_back(make_pair(Base::fromCharacter(ch), 2UL));
        }
        const vector<uint64_t> w(3, 1);
        msaRepair(a, ac, c, w, {}, msa1TestOptions(trigger));
        SHASTA2_ASSERT(msa1ToString(c) == sequence);
        SHASTA2_ASSERT(msa1Ungap(a[2]).empty());
    }



    // The same row in the majority, and now the deletion is the answer.
    {
        const string sequence = "AAAAAAAAAAAAGAAAAAAAAAAA";
        const vector<string> rows = {sequence, string(24, '-'), string(24, '-')};
        vector< vector<AlignedBase> > a = msaRepairAlignmentFromStrings(rows);
        vector<AlignedBase> ac(24, AlignedBase::gap());
        vector< pair<Base, uint64_t> > c;
        const vector<uint64_t> w(3, 1);
        msaRepair(a, ac, c, w, {}, msa1TestOptions(trigger));
        SHASTA2_ASSERT(c.empty());
        SHASTA2_ASSERT(msa1Ungap(a[0]) == sequence);
    }



    // Adversarial cases, and then many random ones. What is checked in both is
    // the set of invariants that must hold whatever the input:
    //  - every row still has the length of the aligned consensus
    //  - no read was altered, only the columns it occupies
    //  - the consensus is still the aligned consensus with the gaps removed
    //  - the alignment was not made worse
    //  - if nothing was repaired, nothing at all changed
    {
        // Build an alignment from rows, give it a column majority consensus with
        // recognisable coverage, repair it, and check the invariants.
        const auto checkOne = [&](const vector<string>& rows, MsaRepairTrigger t) -> bool
        {
            vector< vector<AlignedBase> > a = msaRepairAlignmentFromStrings(rows);
            const vector<uint64_t> w(a.size(), 1);

            vector<AlignedBase> ac;
            vector< pair<Base, uint64_t> > c;
            msa1ColumnConsensus(a, w, ac, c);

            // Give every base a distinct coverage that no vote could produce.
            // The splice must carry these through untouched outside a repaired
            // region, and a sentinel catches a rebuild that recomputes them.
            for(uint64_t i=0; i<c.size(); i++) {
                c[i].second = 500 + i;
            }

            const auto aBefore = a;
            const auto cBefore = c;
            vector<string> readsBefore;
            for(const auto& row: a) {
                readsBefore.push_back(msa1Ungap(row));
            }
            const uint64_t impureBefore = msaRepairImpureColumnCount(a);

            const uint64_t repaired = msaRepair(a, ac, c, w, {}, msa1TestOptions(t));

            for(const auto& row: a) {
                SHASTA2_ASSERT(row.size() == ac.size());
            }
            for(uint64_t i=0; i<a.size(); i++) {
                SHASTA2_ASSERT(msa1Ungap(a[i]) == readsBefore[i]);
            }
            SHASTA2_ASSERT(msa1Ungap(ac) == msa1ToString(c));
            const uint64_t impureAfter = msaRepairImpureColumnCount(a);
            SHASTA2_ASSERT(impureAfter <= impureBefore);
            if(repaired == 0) {
                SHASTA2_ASSERT(a == aBefore);
                SHASTA2_ASSERT(c == cBefore);
            }
            return repaired > 0;
        };

        const string A12(12, 'A'), A11(11, 'A'), A10(10, 'A'), A6(6, 'A');
        const vector< vector<string> > cases = {
            // The pattern flush against each end of the alignment, and filling it.
            {A12 + "G" + A11 + "CGT", A10 + "--G" + A11 + "CGT", A12 + "G" + A11 + "CGT"},
            {"CGT" + A12 + "G" + A11, "CGT" + A10 + "--G" + A11, "CGT" + A12 + "G" + A11},
            {A12 + "G" + A11, A10 + "--G" + A11, A12 + "G" + A11},
            // Two patterns, close enough to merge and far enough to stay apart.
            {A12 + "G" + A11 + "C" + A12 + "G" + A11,
             A10 + "--G" + A11 + "C" + A10 + "--G" + A11,
             A12 + "G" + A11 + "C" + A12 + "G" + A11},
            {A12 + "G" + A11 + string(40, 'C') + A12 + "G" + A11,
             A10 + "--G" + A11 + string(40, 'C') + A10 + "--G" + A11,
             A12 + "G" + A11 + string(40, 'C') + A12 + "G" + A11},
            // Rows that contribute nothing, or much less than the others.
            {A12 + "G" + A11, A10 + "--G" + A11, string(24, '-')},
            {A12 + "G" + A11, A12 + "G" + A11, "-----" + A6 + "G" + A6 + "------"},
            // Very long runs.
            {string(100, 'A') + "G" + string(100, 'A'),
             string(98, 'A') + "--G" + string(100, 'A'),
             string(100, 'A') + "G" + string(100, 'A')},
            // A real substitution where the separating base is, which is not a
            // homopolymer problem and must be left alone.
            {A12 + "G" + A11, A12 + "T" + A11, A12 + "G" + A11},
            // Three long runs in a row, and two patterns sharing a run.
            {A12 + "G" + string(11, 'C') + "T" + string(12, 'G'),
             A10 + "--G" + string(11, 'C') + "T" + string(12, 'G'),
             A12 + "G" + string(11, 'C') + "T" + string(12, 'G')},
            {A12 + "G" + A11 + "G" + A12, A10 + "--G" + A11 + "G" + A12, A12 + "G" + A11 + "G" + A12},
            // Small and degenerate shapes.
            {A12 + "G" + A11, A10 + "--G" + A11},
            {A12 + "G" + A11},
            {A12 + "G" + A11, A12 + "G" + A11, A12 + "G" + A11},
            // Runs sitting on the threshold, a gap column inside a run, a long
            // gap stretch, and every row a different run length.
            {"CC" + string(4, 'A') + "G" + string(4, 'A') + "CC",
             "CC" + string(5, 'A') + "G" + string(4, 'A') + "C",
             "CC" + string(4, 'A') + "G" + string(4, 'A') + "CC"},
            {A6 + "-" + A6 + "G" + A11, A6 + "-" + A6 + "G" + A11, A6 + "A" + A6 + "G" + A11},
            {A12 + "G" + string(10, '-') + A11, A10 + "--G" + string(10, '-') + A11,
             A12 + "G" + string(10, '-') + A11},
            {string(9, 'A') + "---G" + A11, string(10, 'A') + "--G" + A11,
             string(11, 'A') + "-G" + A11, string(12, 'A') + "G" + A11},
            // Reads constrained on one side only. Theseus pads the part of the
            // alignment such a read does not cover with the same '-' it uses for
            // a deletion, so these rows look like rows full of deletions and are
            // not. They must not vote where they do not reach.
            {A12 + "G" + A11, A12 + "G" + A11, string(13, '-') + A11},
            {A12 + "G" + A11, string(13, '-') + A11, string(13, '-') + A11},
            {A12 + "G" + A11, A12 + string(12, '-'), A12 + string(12, '-')},
            {A12 + "G" + A11, A12 + string(12, '-'), string(13, '-') + A11},
            {A12 + "G" + A11, string(14, '-') + string(10, 'A'),
             string(14, '-') + string(10, 'A')}
        };
        for(const MsaRepairTrigger t: {MsaRepairTrigger::PatternOnly, MsaRepairTrigger::AnyLongRun}) {
            uint64_t repairedCases = 0;
            for(const vector<string>& rows: cases) {
                if(checkOne(rows, t)) {
                    ++repairedCases;
                }
            }
            cout << "Checked " << cases.size() << " adversarial alignments with the " <<
                ((t == MsaRepairTrigger::PatternOnly) ? "pattern" : "any long run") <<
                " trigger, " << repairedCases << " of them repaired." << endl;
            SHASTA2_ASSERT(repairedCases > 0);
        }



        // Reads constrained on one side only, checked on the consensus rather
        // than only on the invariants.
        //
        // Theseus pads the part of the alignment a read does not cover with the
        // same '-' it uses for a deletion, and nothing in the alignment tells
        // the two apart. Two padded rows against one that reaches will carry a
        // column if the padding is allowed to vote, and the consensus then loses
        // bases that every read covering them agrees on.
        {
            const string truth = A12 + "G" + A11;
            const auto check = [&](
                const vector<string>& rows,
                const vector<Anchoring>& anchoring,
                const string& expected) {
                vector< vector<AlignedBase> > a = msaRepairAlignmentFromStrings(rows);
                const vector<uint64_t> w(a.size(), 1);

                // Plain column majority, which is what an aligner that counts
                // padding as a gap produces.
                vector<AlignedBase> ac;
                vector< pair<Base, uint64_t> > c;
                msa1ColumnConsensus(a, w, ac, c);

                msaRepair(a, ac, c, w, anchoring);
                SHASTA2_ASSERT(msa1ToString(c) == expected);
            };
            const Anchoring both = Anchoring::BothSides;
            const Anchoring left = Anchoring::LeftOnly;
            const Anchoring right = Anchoring::RightOnly;

            // One padded row, then a padded majority, then the mirror image, then
            // both kinds at once. In every case the answer is the whole sequence.
            check({truth, truth, string(13, '-') + A11},
                {both, both, right}, truth);
            check({truth, string(13, '-') + A11, string(13, '-') + A11},
                {both, right, right}, truth);
            check({truth, A12 + string(12, '-'), A12 + string(12, '-')},
                {both, left, left}, truth);
            check({truth, A12 + string(12, '-'), string(13, '-') + A11},
                {both, left, right}, truth);

            // Two reads that start inside the second run, so all they carry of it
            // is the 10 bases after their own start. That says the run is at
            // least 10, not that it is 10, so it must not vote: counting it would
            // give mode{11, 10, 10} = 10 and shorten a run on the strength of
            // where the reads happen to begin.
            check({truth, string(14, '-') + string(10, 'A'),
                string(14, '-') + string(10, 'A')}, {both, right, right}, truth);

            // No read spans the window, so there is nothing to align the rest
            // against and the aligner's own placement has to stand.
            check({A12 + string(12, '-'), string(13, '-') + A11},
                {left, right}, A12 + A11);

            // The same rows, but every read really is fixed on both sides, as
            // abpoa always reports. Now the gaps are deletions and a majority of
            // deletions is a deletion. The two readings differ, which is why the
            // caller has to say which one applies.
            check({truth, string(13, '-') + A11, string(13, '-') + A11},
                {both, both, both}, A11);

            cout << "Reads constrained on one side only do not vote outside "
                "the part of the alignment they cover." << endl;
        }



        // Random alignments. The generator makes rows that differ only in the
        // lengths of their long runs, which is the case the repair targets, then
        // sprinkles gaps to make an alignment.
        uint64_t random = 1;
        const auto next = [&random]() {
            random = random * 6364136223846793005ULL + 1442695040888963407ULL;
            return uint64_t(random >> 33);
        };
        uint64_t randomRepaired = 0;
        const uint64_t trialCount = 3000;
        for(uint64_t trial=0; trial<trialCount; trial++) {
            random = 1000003ULL * (trial + 1);

            const uint64_t rowCount = 1 + next() % 6;
            const uint64_t blockCount = 1 + next() % 5;
            vector<string> ungapped(rowCount);
            uint64_t previousBase = 4;
            for(uint64_t b=0; b<blockCount; b++) {
                uint64_t baseIndex = next() % 4;
                if(baseIndex == previousBase) {
                    baseIndex = (baseIndex + 1) % 4;
                }
                previousBase = baseIndex;
                const char base = "ACGT"[baseIndex];
                const bool isLong = (next() % 2) == 0;
                for(uint64_t i=0; i<rowCount; i++) {
                    const uint64_t runLength = isLong ? (5 + next() % 10) : (1 + next() % 2);
                    for(uint64_t k=0; k<runLength; k++) {
                        ungapped[i].push_back(base);
                    }
                }
                if(b + 1 < blockCount) {
                    uint64_t separatorIndex = next() % 4;
                    if(separatorIndex == previousBase) {
                        separatorIndex = (separatorIndex + 1) % 4;
                    }
                    for(uint64_t i=0; i<rowCount; i++) {
                        ungapped[i].push_back("ACGT"[separatorIndex]);
                    }
                    previousBase = separatorIndex;
                }
            }

            uint64_t width = 0;
            for(const string& s: ungapped) {
                width = max(width, uint64_t(s.size()));
            }
            width += next() % 5;
            vector<string> rows;
            for(uint64_t i=0; i<rowCount; i++) {
                string padded = ungapped[i];
                while(padded.size() < width) {
                    const uint64_t p = next() % (padded.size() + 1);
                    padded.insert(padded.begin() + int64_t(p), '-');
                }
                rows.push_back(padded);
            }
            if(checkOne(rows, (trial % 2) ? MsaRepairTrigger::AnyLongRun : MsaRepairTrigger::PatternOnly)) {
                ++randomRepaired;
            }
        }
        cout << "Checked " << trialCount << " random alignments, " <<
            randomRepaired << " of them repaired." << endl;
        SHASTA2_ASSERT(randomRepaired > 0);
    }



    // Degenerate inputs.
    {
        vector< vector<AlignedBase> > a;
        vector<AlignedBase> ac;
        vector< pair<Base, uint64_t> > c;
        SHASTA2_ASSERT(msaRepair(a, ac, c, {}, {}, msa1TestOptions(trigger)) == 0);

        // A single row cannot disagree with the consensus, so nothing is found.
        a.push_back(vectorOfAlignedBasesFromString(
            string(12, 'A') + "G" + string(11, 'A')));
        ac = a.front();
        c.clear();
        for(const AlignedBase b: ac) {
            c.push_back(make_pair(Base(b), 1UL));
        }
        const auto aBefore = a;
        SHASTA2_ASSERT(msaRepair(a, ac, c, {}, {}, msa1TestOptions(trigger)) == 0);
        SHASTA2_ASSERT(a == aBefore);
    }

    cout << "testMsa1Repair passed." << endl;
}
