#pragma once

/*******************************************************************************

A TangleMatrix is defined by two sets of AssemblyGraph edges (segments):
- The entrances.
- The exits.
Any two sets of edges can be used to define a TangleMatrix,
regardless of their location and connectivity in the AssemblyGraph.
The two sets are not required to be disjoint.

We use the second to last AnchorId of each entrance and the second
AnchorId of each exit to count common reads.

*******************************************************************************/

// Shasta.
#include "AssemblyGraph.hpp"

// Standard library.
#include <vector.hpp>

namespace shasta {
    class TangleMatrix;

}


class shasta::TangleMatrix {
public:
    using edge_descriptor = AssemblyGraph::edge_descriptor;

    TangleMatrix(
        const AssemblyGraph&,
        vector<edge_descriptor> entranceEdges,
        vector<edge_descriptor> exitEdges);

    TangleMatrix() {}

    void construct(
        const AssemblyGraph&,
        vector<edge_descriptor> entranceEdges,
        vector<edge_descriptor> exitEdges);

    void writeHtml(
        const AssemblyGraph&,
        ostream&) const;

    class Entrance {
        public:
        edge_descriptor e;

        // The second to last AnchorId of this AssemblyGraphEdge.
        AnchorId anchorId;

        // Total coverage for this Anchor.
        uint64_t coverage;

        // Common coverage for this Anchor (sum of tangle matrix row).
        uint64_t commonCoverage = 0;

        Entrance(
            edge_descriptor e,
            AnchorId anchorId,
            uint64_t coverage) :
            e(e), anchorId(anchorId), coverage(coverage) {}
    };
    vector<Entrance> entrances;



    class Exit {
        public:
        edge_descriptor e;

        // The second to AnchorId of this AssemblyGraphEdge.
        AnchorId anchorId;

        // Total coverage for this Anchor.
        uint64_t coverage;

        // Common coverage for this Anchor (sum of tangle matrix column).
        uint64_t commonCoverage = 0;

        Exit(
            edge_descriptor e,
            AnchorId anchorId,
            uint64_t coverage) :
            e(e), anchorId(anchorId), coverage(coverage) {}
    };
    vector<Exit> exits;

    // The tangle matrix is indexed by [iEntrance][iExit].
    vector < vector<uint64_t> > tangleMatrix;
};
