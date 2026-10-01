// Shasta.
#include "HomopolymerModel.hpp"
#include "invalid.hpp"
#include "msaRepair.hpp"
#include "SHASTA2_ASSERT.hpp"
#include "tmpDirectory.hpp"
using namespace shasta2;

// Standard library.
#include "algorithm.hpp"
#include "fstream.hpp"
#include "iostream.hpp"
#include "memory.hpp"
#include "stdexcept.hpp"
#include <cmath>
#include <filesystem>
#include <limits>
#include <sstream>
#include <unistd.h>

const uint64_t HomopolymerModel::unknownFlank = 4;


// See HomopolymerModel.hpp for the file format.
HomopolymerModel::HomopolymerModel(istream& file)
{
    SHASTA2_ASSERT(file);

    const string expectedHeader = "base,strand,left,right,n,m,count,probability,logProbabilityDb";
    string line;
    std::getline(file, line);
    if(not line.empty() and line.back() == '\r') {
        line.pop_back();
    }
    if(line != expectedHeader) {
        throw runtime_error("Homopolymer model: invalid format. Missing expected header " + expectedHeader);
    }

    // Read all the lines first, because the size of the table is only known
    // at the end.
    class Line {
    public:
        uint64_t base;
        uint64_t strand;
        uint64_t left;
        uint64_t right;
        uint64_t n;
        uint64_t m;
        uint64_t count;
        double logProbabilityDb;
    };
    vector<Line> lines;
    uint64_t lineNumber = 1;
    while(std::getline(file, line)) {
        ++lineNumber;
        if(not line.empty() and line.back() == '\r') {
            line.pop_back();
        }
        if(line.empty()) {
            continue;
        }
        const string where = "line " + to_string(lineNumber);

        vector<string> tokens;
        std::istringstream s(line);
        string token;
        while(std::getline(s, token, ',')) {
            tokens.push_back(token);
        }
        if(tokens.size() != 9) {
            throw runtime_error("Expected 9 fields at " + where);
        }

        Line l;
        try {
            l.base = flankFromString(tokens[0]);
        } catch(const std::exception&) {
            throw runtime_error("Invalid base at " + where);
        }
        if(l.base == unknownFlank) {
            throw runtime_error("Invalid base at " + where);
        }
        try {
            l.left = flankFromString(tokens[2]);
            l.right = flankFromString(tokens[3]);
        } catch(const std::exception&) {
            throw runtime_error("Invalid flanking base at " + where);
        }
        try {
            l.strand = std::stoull(tokens[1]);
            l.n = std::stoull(tokens[4]);
            l.m = std::stoull(tokens[5]);
            l.count = std::stoull(tokens[6]);
            l.logProbabilityDb = std::stod(tokens[8]);
        } catch(const std::exception&) {
            throw runtime_error("Invalid number at " + where);
        }
        if(l.strand > 1) {
            throw runtime_error("Invalid strand at " + where);
        }
        if(l.n < 2) {
            throw runtime_error("n must be at least 2 at " + where);
        }
        if(not (l.logProbabilityDb <= 0.)) {
            throw runtime_error("logProbabilityDb must be finite and not positive at " + where);
        }
        nCount = max(nCount, l.n + 1);
        mCount = max(mCount, l.m + 1);
        lines.push_back(l);
    }
    if(lines.empty()) {
        throw runtime_error("Homopolymer model: invalid format. Definition is empty.");
    }

    // Fill the table. A NaN marks an entry not yet seen, so a missing or
    // duplicate line can be detected.
    const double missing = std::numeric_limits<double>::quiet_NaN();
    for(auto& strands: table) {
        for(auto& lefts: strands) {
            for(auto& rights: lefts) {
                for(auto& rows: rights) {
                    rows.resize(nCount);
                }
            }
        }
    }
    const string flankCharacters = "ACGT*";
    for(const Line& l: lines) {
        vector<double>& row = table[l.base][l.strand][l.left][l.right][l.n];
        if(row.empty()) {
            row.assign(mCount, missing);
        }
        if(not std::isnan(row[l.m])) {
            throw runtime_error("Homopolymer model format error. Duplicate line for " +
                string(1, flankCharacters[l.base]) + " strand " + to_string(l.strand) +
                " flanks " + flankCharacters[l.left] + flankCharacters[l.right] +
                " n " + to_string(l.n) + " m " + to_string(l.m));
        }
        row[l.m] = l.logProbabilityDb;
    }

    // Every row present must be complete. Every base and strand must have *
    // rows, and a row for given flanks is only allowed for an n that also has
    // a * row, because the * rows decide which n can be chosen.
    for(uint64_t base=0; base<4; base++) {
        for(uint64_t strand=0; strand<2; strand++) {
            const string what = string(1, flankCharacters[base]) + " strand " + to_string(strand);
            bool found = false;
            for(uint64_t left=0; left<5; left++) {
                for(uint64_t right=0; right<5; right++) {
                    for(uint64_t n=0; n<nCount; n++) {
                        const vector<double>& r = table[base][strand][left][right][n];
                        if(r.empty()) {
                            continue;
                        }
                        const bool isStar = (left == unknownFlank) and (right == unknownFlank);
                        found = found or isStar;
                        if(not isStar and table[base][strand][unknownFlank][unknownFlank][n].empty()) {
                            throw runtime_error("Homopolymer model format error."
                                    "Flank row without a * row for " + what +
                                " n " + to_string(n));
                        }
                        for(uint64_t m=0; m<mCount; m++) {
                            if(std::isnan(r[m])) {
                                throw runtime_error("Homopolymer model format error. Missing line for " + what +
                                    " flanks " + flankCharacters[left] + flankCharacters[right] +
                                    " n " + to_string(n) + " m " + to_string(m));
                            }
                        }
                    }
                }
            }
            if(not found) {
                throw runtime_error("Homopolymer model format error. No * lines for " + what);
            }
        }
    }

    // The prior: how often each n occurs in the training data, from the *
    // rows, pooled over bases and strands, with one added to the count of
    // every n present.
    vector<double> priorCount(nCount, 0.);
    for(const Line& l: lines) {
        if((l.left == unknownFlank) and (l.right == unknownFlank)) {
            priorCount[l.n] += double(l.count);
        }
    }
    double priorTotal = 0.;
    for(uint64_t n=0; n<nCount; n++) {
        bool present = false;
        for(uint64_t base=0; base<4; base++) {
            present = present or not table[base][0][unknownFlank][unknownFlank][n].empty();
        }
        if(present) {
            priorCount[n] += 1.;
        }
        priorTotal += priorCount[n];
    }
    logPriorDb.resize(nCount);
    for(uint64_t n=0; n<nCount; n++) {
        logPriorDb[n] = (priorCount[n] > 0.) ?
            10. * std::log10(priorCount[n] / priorTotal) :
            -std::numeric_limits<double>::infinity();
    }
}



uint64_t HomopolymerModel::flankFromString(const string& s)
{
    if(s == "*") {
        return unknownFlank;
    }
    if(s.size() != 1) {
        throw runtime_error("Invalid flanking base " + s);
    }
    return Base::fromCharacter(s[0]).value;
}



const vector<double>* HomopolymerModel::row(
    uint64_t base, uint64_t strand, uint64_t left, uint64_t right, uint64_t n) const
{
    SHASTA2_ASSERT(base < 4);
    SHASTA2_ASSERT(strand < 2);
    SHASTA2_ASSERT(left <= unknownFlank);
    SHASTA2_ASSERT(right <= unknownFlank);
    if(n >= nCount) {
        return nullptr;
    }
    const vector<double>& flankRow = table[base][strand][left][right][n];
    if(not flankRow.empty()) {
        return &flankRow;
    }
    const vector<double>& starRow = table[base][strand][unknownFlank][unknownFlank][n];
    return starRow.empty() ? nullptr : &starRow;
}



double HomopolymerModel::logProbabilityDb(
    Base base, uint64_t strand, uint64_t left, uint64_t right,
    uint64_t n, uint64_t m) const
{
    const vector<double>* r = row(base.value, strand, left, right, n);
    SHASTA2_ASSERT(r);
    return (*r)[min(m, maxM())];
}



uint64_t HomopolymerModel::mostLikelyLength(
    Base base,
    uint64_t left,
    uint64_t right,
    const array<vector<uint64_t>, 2>& observedLengths,
    vector<double>* logPosteriorDb) const
{
    const double minusInfinity = -std::numeric_limits<double>::infinity();

    // log P(n) + sum over observations of log P(m | n, ...), in dB.
    vector<double> logLikelihood(nCount, minusInfinity);
    uint64_t bestN = invalid<uint64_t>;
    for(uint64_t n=0; n<nCount; n++) {
        const array<const vector<double>*, 2> rows = {
            row(base.value, 0, left, right, n),
            row(base.value, 1, left, right, n)};
        if(not rows[0] or not rows[1]) {
            continue;
        }
        double sum = logPriorDb[n];
        for(uint64_t strand=0; strand<2; strand++) {
            for(const uint64_t m: observedLengths[strand]) {
                sum += (*rows[strand])[min(m, maxM())];
            }
        }
        logLikelihood[n] = sum;
        if((bestN == invalid<uint64_t>) or (sum > logLikelihood[bestN])) {
            bestN = n;
        }
    }
    SHASTA2_ASSERT(bestN != invalid<uint64_t>);

    if(logPosteriorDb) {
        // Normalize so the probabilities sum to 1. Subtracting the best
        // value first keeps every term of the sum at most 1.
        const double best = logLikelihood[bestN];
        double sum = 0.;
        for(const double x: logLikelihood) {
            if(x != minusInfinity) {
                sum += std::pow(10., (x - best) / 10.);
            }
        }
        const double logSum = best + 10. * std::log10(sum);
        logPosteriorDb->resize(nCount);
        for(uint64_t n=0; n<nCount; n++) {
            (*logPosteriorDb)[n] = (logLikelihood[n] == minusInfinity) ?
                minusInfinity : logLikelihood[n] - logSum;
        }
    }

    return bestN;
}



// A small hand-built model: true lengths 2 and 3, observed lengths 0 to 4,
// the same for every base. In the * rows, strand 0 reads a run correctly 90%
// of the time, and strand 1 reads it one base short 90% of the time. Both
// lengths have the same count, so the prior does not favor either. For base
// A between C and G only, strand 0 also reads one base short.
void shasta2::testHomopolymerModel()
{
    const string fileName = tmpDirectory() + "testHomopolymerModel-" + to_string(::getpid()) + ".csv";
    {
        ofstream csv(fileName);
        csv << "base,strand,left,right,n,m,count,probability,logProbabilityDb\n";
        const auto write = [&](char base, uint64_t strand, const string& left, const string& right,
            uint64_t n, uint64_t peak, uint64_t count) {
            for(uint64_t m=0; m<=4; m++) {
                const double p = (m == peak) ? 0.9 : 0.025;
                csv << base << "," << strand << "," << left << "," << right << "," << n << "," << m << "," <<
                    ((m == peak) ? count : 0) << "," << p << "," << 10. * std::log10(p) << "\n";
            }
        };
        for(const char base: string("ACGT")) {
            for(uint64_t strand=0; strand<2; strand++) {
                for(uint64_t n=2; n<=3; n++) {
                    write(base, strand, "*", "*", n, (strand == 0) ? n : n - 1, 100);
                }
            }
        }
        for(uint64_t n=2; n<=3; n++) {
            write('A', 0, "C", "G", n, n - 1, 0);
        }
    }
    ifstream file(fileName);
    const auto modelPointer = make_shared<const HomopolymerModel>(file);
    std::filesystem::remove(fileName);

    const Base A = Base::fromCharacter('A');
    const uint64_t C = Base::fromCharacter('C').value;
    const uint64_t G = Base::fromCharacter('G').value;
    const uint64_t T = Base::fromCharacter('T').value;
    const uint64_t any = HomopolymerModel::unknownFlank;
    SHASTA2_ASSERT(HomopolymerModel::flankFromString("C") == C);
    SHASTA2_ASSERT(HomopolymerModel::flankFromString("*") == any);
    SHASTA2_ASSERT(modelPointer->maxN() == 3);
    SHASTA2_ASSERT(modelPointer->maxM() == 4);
    SHASTA2_ASSERT(std::abs(modelPointer->logProbabilityDb(A, 0, any, any, 3, 3) - 10. * std::log10(0.9)) < 1e-5);

    // An m past the largest in the model is treated as the largest.
    SHASTA2_ASSERT(modelPointer->logProbabilityDb(A, 0, any, any, 3, 100) == modelPointer->logProbabilityDb(A, 0, any, any, 3, 4));

    // Strand is used. A strand 0 read of 3 and a strand 1 read of 2 both say 3.
    // Scored as if both were on strand 0 they would tie, and the tie would go to 2.
    SHASTA2_ASSERT(modelPointer->mostLikelyLength(A, any, any, {vector<uint64_t>{3}, vector<uint64_t>{2}}) == 3);
    SHASTA2_ASSERT(modelPointer->mostLikelyLength(A, any, any, {vector<uint64_t>{3, 2}, vector<uint64_t>{}}) == 2);
    SHASTA2_ASSERT(modelPointer->mostLikelyLength(A, any, any, {vector<uint64_t>{}, vector<uint64_t>{2, 2}}) == 3);

    // The flanks are used. Between C and G a strand 0 read of 2 says 3.
    // Between other flanks, which have no row of their own, the * row is used
    // and the same read says 2.
    SHASTA2_ASSERT(modelPointer->mostLikelyLength(A, C, G, {vector<uint64_t>{2, 2}, vector<uint64_t>{}}) == 3);
    SHASTA2_ASSERT(modelPointer->mostLikelyLength(A, G, C, {vector<uint64_t>{2, 2}, vector<uint64_t>{}}) == 2);
    SHASTA2_ASSERT(modelPointer->mostLikelyLength(A, T, G, {vector<uint64_t>{2, 2}, vector<uint64_t>{}}) == 2);

    // The posterior sums to 1, and lengths absent from the model are impossible.
    vector<double> logPosteriorDb;
    SHASTA2_ASSERT(modelPointer->mostLikelyLength(A, any, any, {vector<uint64_t>{3, 2}, vector<uint64_t>{2}}, &logPosteriorDb) == 3);
    SHASTA2_ASSERT(logPosteriorDb.size() == 4);
    SHASTA2_ASSERT(logPosteriorDb[0] == -std::numeric_limits<double>::infinity());
    SHASTA2_ASSERT(logPosteriorDb[1] == -std::numeric_limits<double>::infinity());
    SHASTA2_ASSERT(std::abs(std::pow(10., logPosteriorDb[2] / 10.) + std::pow(10., logPosteriorDb[3] / 10.) - 1.) < 1e-9);

    // Many reads do not underflow. The product of the probabilities for n = 2
    // would be 0.025^2000, far below the smallest double.
    const vector<uint64_t> manyReads(2000, 3);
    SHASTA2_ASSERT(modelPointer->mostLikelyLength(A, any, any, {manyReads, vector<uint64_t>{}}, &logPosteriorDb) == 3);
    SHASTA2_ASSERT(std::abs(logPosteriorDb[3]) < 1e-9);
    SHASTA2_ASSERT(std::isfinite(logPosteriorDb[2]) and logPosteriorDb[2] < -1000.);

    // The same through msa1's consensus, which must pass the strand of each
    // row. One poly A column, a strand 0 row of length 3 and a strand 1 row of
    // length 2. With the model this is 3. The Mode estimator sees a tie and
    // takes 2.
    const AlignedExtendedBase polyA(ExtendedBase::polyFromBase(A));
    const vector<uint64_t> weights = {1, 1};
    const vector< array<uint64_t, 2> > strandWeights = {{1, 0}, {0, 1}};
    {
        const vector<AlignedExtendedSequence> alignment = {
            AlignedExtendedSequence(1, make_pair(polyA, 3UL)),
            AlignedExtendedSequence(1, make_pair(polyA, 2UL))};
        const vector< pair<uint64_t, uint64_t> > spans(2, make_pair(0UL, 1UL));
        vector< pair<Base, uint64_t> > consensus;
        AlignedExtendedSequence alignedConsensus;
        extendedConsensus(alignment, weights, RunLengthEstimator::Mode, spans,
            consensus, alignedConsensus, modelPointer, strandWeights);
        SHASTA2_ASSERT(alignedConsensus.size() == 1 and alignedConsensus.front().second == 3);
        extendedConsensus(alignment, weights, RunLengthEstimator::Mode, spans,
            consensus, alignedConsensus);
        SHASTA2_ASSERT(alignedConsensus.size() == 1 and alignedConsensus.front().second == 2);
    }

    // msa1's consensus passes the flanking consensus bases, in the order they
    // appear in the assembly, skipping gap columns. Columns C, gap, poly A, G,
    // with two strand 0 rows of length 2: between C and G this is 3. The same
    // with the flanks swapped, G and C, has no row of its own and gives 2.
    const auto flankedLength = [&](char left, char right) {
        const AlignedExtendedBase l(ExtendedBase::fromCharacter(left));
        const AlignedExtendedBase r(ExtendedBase::fromCharacter(right));
        const AlignedExtendedSequence row = {
            make_pair(l, 1UL), make_pair(AlignedExtendedBase::gap(), 0UL),
            make_pair(polyA, 2UL), make_pair(r, 1UL)};
        const vector<AlignedExtendedSequence> alignment = {row, row};
        const vector< pair<uint64_t, uint64_t> > spans(2, make_pair(0UL, 4UL));
        vector< pair<Base, uint64_t> > consensus;
        AlignedExtendedSequence alignedConsensus;
        extendedConsensus(alignment, weights, RunLengthEstimator::Mode, spans,
            consensus, alignedConsensus, modelPointer, {{1, 0}, {1, 0}});
        return alignedConsensus[2].second;
    };
    SHASTA2_ASSERT(flankedLength('C', 'G') == 3);
    SHASTA2_ASSERT(flankedLength('G', 'C') == 2);

    // When a flank is not known the * row is used, even if the other flank
    // is known: C with an unknown right flank, or an unknown left flank with
    // G, is not the C to G context and gives 2.
    SHASTA2_ASSERT(modelPointer->mostLikelyLength(A, C, any, {vector<uint64_t>{2, 2}, vector<uint64_t>{}}) == 2);
    SHASTA2_ASSERT(modelPointer->mostLikelyLength(A, any, G, {vector<uint64_t>{2, 2}, vector<uint64_t>{}}) == 2);

    // In msa1's consensus, a run at the beginning or end of the alignment has
    // no flank on that side, and neither does a run whose only neighbors on
    // that side are gap columns. Each of these gives 2, not the 3 of C to G.
    const auto edgeLength = [&](const string& columns) {
        // One character per column: A is the poly A run, - is a gap, and
        // anything else is that plain base.
        AlignedExtendedSequence row;
        uint64_t runColumn = invalid<uint64_t>;
        for(const char c: columns) {
            if(c == 'A') {
                runColumn = row.size();
                row.push_back(make_pair(polyA, 2UL));
            } else if(c == '-') {
                row.push_back(make_pair(AlignedExtendedBase::gap(), 0UL));
            } else {
                row.push_back(make_pair(AlignedExtendedBase(ExtendedBase::fromCharacter(c)), 1UL));
            }
        }
        const vector<AlignedExtendedSequence> alignment = {row, row};
        const vector< pair<uint64_t, uint64_t> > spans(2, make_pair(0UL, uint64_t(row.size())));
        vector< pair<Base, uint64_t> > consensus;
        AlignedExtendedSequence alignedConsensus;
        extendedConsensus(alignment, weights, RunLengthEstimator::Mode, spans,
            consensus, alignedConsensus, modelPointer, {{1, 0}, {1, 0}});
        return alignedConsensus[runColumn].second;
    };
    SHASTA2_ASSERT(edgeLength("CAG") == 3);     // Both flanks: C to G.
    SHASTA2_ASSERT(edgeLength("AG") == 2);      // Run at the beginning.
    SHASTA2_ASSERT(edgeLength("CA") == 2);      // Run at the end.
    SHASTA2_ASSERT(edgeLength("--AG") == 2);    // Only gaps to the left.
    SHASTA2_ASSERT(edgeLength("CA--") == 2);    // Only gaps to the right.
    SHASTA2_ASSERT(edgeLength("A") == 2);       // No flank on either side.

    cout << "testHomopolymerModel passed." << endl;
}
