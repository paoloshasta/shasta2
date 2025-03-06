#pragma once

// Shasta.
#include <cstdint.hpp>
#include <string.hpp>
#include <vector.hpp>

// CLI11. Requires package libcli11-dev.
#include "CLI/CLI.hpp"

namespace shasta {
    class AssemblerOptions;
}


class shasta::AssemblerOptions : public CLI::App {
public:

    AssemblerOptions(int argumentCount, char** arguments);
    bool isHelp = false;

    string configName;
    vector <string> inputFileNames;
    string assemblyDirectory = "ShastaRun";
    string command = "assemble";
    string memoryMode = "anonymous";
    string memoryBacking= "4K";
    uint64_t threadCount = 0;
    string exploreAccess = "user";
    uint16_t port = 17100;

    uint64_t minReadLength = 0;

    uint64_t k = 60;
    double markerDensity = 0.05;

    uint64_t minAnchorCoverage = 10;
    uint64_t maxAnchorCoverage = 60;

    class LocalAssemblyOptions {
    public:

        // The estimated offset gets extended by this ratio to
        // decide how much to extend reads that only appear in edgeIdA or edgeIdB.
        double estimatedOffsetRatio = 1.1;

        // Vertex sampling rate, used to set minVertexCoverage.
        // Only used if minVertexCoverage is 0 on input to
        // mode3::LocalAssembly constructor.
        double vertexSamplingRate = 0.8;

        // Alignment parameters.
        int64_t matchScore = 6;
        int64_t mismatchScore = -1;
        int64_t gapScore = -1;

        // Number of bases (not markers) that can be skipped by an alignment.
        uint64_t maxSkipBases = 500;

        // The maximum tolerated length drift of each read.
        // Used to compute the band for banded alignments.
        double maxDrift = 0.005;

        // Minimum half band, in markers.
        uint64_t minHalfBand = 100;

        // Minimum ration of score to best possible score for
        // an alignment to be used.
        double minScoreRatio = 0.7;

        // The maximum length of an MSA alignment we are willing to compute.
        uint64_t maxMsaLength = 5000;

    };
    LocalAssemblyOptions localAssemblyOptions;
};


