#include "AssemblyGraph.hpp"
#include "deduplicate.hpp"
#include "GTest.hpp"
#include "html.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
#include "Tangle.hpp"
#include "TangleMatrix.hpp"
#include "timestamp.hpp"
using namespace shasta2;



// This turns the AssemblyGraph into a "pseudo" Assembly graph
// containing Segments that are are as long as possible,
// but that generally contain haplotype switches and
// other assembly errors (so-called "pseudo-haplotypes").
// This should be called after the AssemblyGraph has
// already been made single-stranded.



void AssemblyGraph::makePseudo()
{
    AssemblyGraph& assemblyGraph = *this;
    ostream noOutput(0);
    performanceLog << timestamp << "AssemblyGraph::makePseudo begins." << endl;
    const bool debug = true;
    if(debug) {
        cout << "AssemblyGraph::makePseudo begins." << endl;
    }

    // Gather candidate Segments.
    vector<Segment> segmentsToBeProcessed;
    BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {
        const vertex_descriptor v0 = source(segment, assemblyGraph);
        if(out_degree(v0, assemblyGraph) != 1) {
            continue;
        }
        if(in_degree(v0, assemblyGraph) != 2) {
            continue;
        }
        const vertex_descriptor v1 = target(segment, assemblyGraph);
        if(in_degree(v1, assemblyGraph) != 1) {
            continue;
        }
        if(out_degree(v1, assemblyGraph) != 2) {
            continue;
        }
        segmentsToBeProcessed.push_back(segment);
    }
    if(debug) {
        cout << "Found " << segmentsToBeProcessed.size() << " segments to be processed." << endl;
    }


    // Process them one by one.
    for(const Segment segment: segmentsToBeProcessed) {
        if(debug) {
            cout << "Working on " << id(segment) << endl;
        }
        const vertex_descriptor v0 = source(segment, assemblyGraph);
        const vertex_descriptor v1 = target(segment, assemblyGraph);

        // Create a Tangle consisting of this edge.
        // If degenerate, skip it (this can happen if an exit is also an entrance).
        const vector<vertex_descriptor> tangleVertices = {v0, v1};
        const Tangle tangle(assemblyGraph, tangleVertices);
        if((tangle.entrances.size() != 2) or(tangle.exits.size() != 2)) {
            if(debug) {
                cout << "Skipping degenerate tangle." << endl;
            }
            continue;
        }

        // Compute the TangleMatrix and run the likelihood ratio test.
        const TangleMatrix tangleMatrix(assemblyGraph, tangle.entrances, tangle.exits, noOutput);
        const GTest gTest(tangleMatrix.tangleMatrix, assemblyGraph.options.detangleEpsilon, true, true);
        SHASTA2_ASSERT(gTest.success);
        SHASTA2_ASSERT(gTest.hypotheses.size() == 2);
        const auto& hypothesis0 = gTest.hypotheses[0];
        const auto& hypothesis1 = gTest.hypotheses[1];
        if(debug) {
            cout << "Hypothesis 0: " <<
                (hypothesis0.connectivityMatrix[0][0] ? "in-phase" : "out-of-phase") <<
                " G=" << hypothesis0.G << endl;
            cout << "Hypothesis 1: " <<
                (hypothesis1.connectivityMatrix[0][0] ? "in-phase" : "out-of-phase") <<
                " G=" << hypothesis1.G << endl;
        }

        // Decide if we connect the entrances to the exits in-phase or
        // out-of-phase.
        const bool inPhase = hypothesis0.connectivityMatrix[0][0];

        if(debug) {
            cout << (inPhase ? "In-phase" : "Out-of-phase") << endl;
        }

        // Now do the two connections.
        for(uint64_t i=0; i<2; i++) {
            const Segment entrance = tangle.entrances[i];
            const Segment exit = tangle.exits[inPhase ? i : (1-i)];
            if(debug) {
                cout << "Connecting " << id(entrance) << " with " << id(exit) << endl;
            }

            // Disconnect the entrance from its target vertex.
            const vertex_descriptor v0 = source(entrance, assemblyGraph);
            const vertex_descriptor v1Old = target(entrance, assemblyGraph);
            const AnchorId anchorId1 = assemblyGraph[v1Old].anchorId;
            const vertex_descriptor v1New = add_vertex(AssemblyGraphVertex(anchorId1, assemblyGraph.nextVertexId++), assemblyGraph);
            add_edge(v0, v1New, assemblyGraph[entrance], assemblyGraph);
            boost::remove_edge(entrance, assemblyGraph);

            // Disconnect the exit from its source vertex.
            const vertex_descriptor v2Old = source(exit, assemblyGraph);
            const vertex_descriptor v3 = target(exit, assemblyGraph);
            const AnchorId anchorId2 = assemblyGraph[v2Old].anchorId;
            const vertex_descriptor v2New = add_vertex(AssemblyGraphVertex(anchorId2, assemblyGraph.nextVertexId++), assemblyGraph);
            add_edge(v2New, v3, assemblyGraph[exit], assemblyGraph);
            boost::remove_edge(exit, assemblyGraph);

            // Make a copy of our tangle Segment.
            // The copy for the first haplotype keeps its id.
            auto[newSegment, _] = add_edge(v1New, v2New, assemblyGraph[segment], assemblyGraph);
            if(i == 1) {
                assemblyGraph[newSegment].id = nextEdgeId++;
            }
        }

        // Now we can remove our Segment.
        boost::remove_edge(segment, assemblyGraph);
    }



    performanceLog << timestamp << "AssemblyGraph::makePseudo ends." << endl;
    if(debug) {
        cout << "AssemblyGraph::makePseudo ends." << endl;
    }
}

