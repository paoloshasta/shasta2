#pragma once

// Shasta.
#include "Base.hpp"
#include "ReadId.hpp"

// Standard library.
#include "array.hpp"
#include "cstdint.hpp"
#include "string.hpp"
#include "vector.hpp"

namespace shasta2 {
    class HomopolymerModel;
    void testHomopolymerModel();
}



// An empirical model of homopolymer length errors: the probability
// P(m | n, base, strand, left, right) that a read reports a homopolymer run
// of length m when its true length is n. Used by msaRepair to choose the length of
// a long homopolymer run from the lengths observed in the reads that cover it.
//
// base is the base of the run as it appears in the assembly, and strand is
// the OrientedReadId strand of the read. A run of base X seen on strand 1 is
// a run of complement(X) as the read was sequenced, so the two strands of the
// same base have different error profiles. left and right are the bases
// next to the run, also as they appear in the assembly. The error profile
// depends on them too: at the same n, the mean length error differs by more
// than half a base between flanking contexts.
//
// The model is read from a csv file with this header:
//
//     base,strand,left,right,n,m,count,probability,logProbabilityDb
//
// one line per (base, strand, left, right, n, m):
// - left and right are A, C, G, T, or * for the model that ignores the
//   flanking bases. The * rows are required: they are used when a flanking
//   base is not known, or no row exists for it.
// - count is the number of training observations. It is only used for the
//   prior P(n), the fraction of training observations with true length n,
//   from the * rows, pooled over bases and strands, with one added to every
//   count.
// - probability is P(m | n, ...), for reading only.
// - logProbabilityDb is the same probability in decibels,
//   10 * log10(probability), which is what the model uses. For example a
//   probability of 0.1 is -10 dB. Working in logarithms means the evidence of
//   many reads is added instead of multiplied, which avoids underflow.
//
// n starts at 2: msaRepair only asks about runs of at least two bases. A value of
// n may be missing entirely for a given base and strand, typically because
// training saw no runs of that length, and it is then never chosen. An n that
// is present must have a line for every m from 0 to the largest m in the
// file. The largest n and the largest m each stand for "that length or more".
class shasta2::HomopolymerModel {
public:

    explicit HomopolymerModel(istream&);

    // The value of a flank argument when the flanking base is not known.
    // Defined in the cpp file.
    static const uint64_t unknownFlank;

    // A flank argument from its character: A, C, G, T, or * for
    // unknownFlank. Throws on anything else.
    static uint64_t flankFromString(const string&);

    // log P(m | n, base, strand, left, right) in dB. left and right are
    // Base values, or unknownFlank. n must be present in the model for this
    // base and strand. An m larger than the largest m in the model is
    // treated as that largest m.
    double logProbabilityDb(
        Base base, uint64_t strand, uint64_t left, uint64_t right,
        uint64_t n, uint64_t m) const;

    // Given the lengths m of a homopolymer run of this base, between these
    // flanking bases, observed in the reads on each strand, return the most
    // likely true length n: the one that maximizes
    // log P(n) + sum of log P(m | n, ...) over the observations.
    // Ties go to the shorter length.
    //
    // If logPosteriorDb is not null, it is also filled with the log
    // probability, in dB, of every n from 0 to maxN() given the observations,
    // normalized so that the probabilities sum to 1. An n absent from the
    // model gets -infinity.
    uint64_t mostLikelyLength(
        Base base,
        uint64_t left,
        uint64_t right,
        const array<vector<uint64_t>, 2>& observedLengths,
        vector<double>* logPosteriorDb = nullptr) const;

    uint64_t maxN() const
    {
        return nCount - 1;
    }
    uint64_t maxM() const
    {
        return mCount - 1;
    }

private:
    uint64_t nCount = 0;
    uint64_t mCount = 0;

    // log P(m | n, ...) in dB, indexed by [base][strand][left][right][n], one
    // entry per m. left and right are Base values, or unknownFlank for the *
    // rows. Empty for an n absent from the model.
    array<array<array<array<vector< vector<double> >, 5>, 5>, 2>, 4> table;

    // The row to use: the one for these flanks if present, otherwise the *
    // row. Returns null if n is absent from the model.
    const vector<double>* row(
        uint64_t base, uint64_t strand, uint64_t left, uint64_t right, uint64_t n) const;

    // log P(n) in dB, indexed by n.
    vector<double> logPriorDb;
};
