// Shasta2.
#include "AssemblyGraph.hpp"
#include "color.hpp"
#include "DisjointSets.hpp"
#include "GTest.hpp"
#include "html.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
#include "Tangle.hpp"
#include "TangleMatrix.hpp"
#include "timestamp.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



void AssemblyGraph::detangleVertices()
{
    performanceLog << timestamp << "AssemblyGraph::detangleVertices begins." << endl;

    AssemblyGraph& assemblyGraph = *this;
    const bool debug = false;
    ostream html(0);

    const vector< vector<bool> > inPhaseConnectivityMatrix    = { {true , false}, {false, true } };
    const vector< vector<bool> > outOfPhaseConnectivityMatrix = { {false, true }, {true , false} };

    // For now we only handle vertices with in-degree and out-degree equal to 2.
    // Gather the ones with id less that the id of their reverse complement.
    vector<vertex_descriptor> candidateVertices;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        if(in_degree(v, assemblyGraph) != 2) {
            continue;
        }
        if(out_degree(v, assemblyGraph) != 2) {
            continue;
        }
        const AssemblyGraphVertex& vertex = assemblyGraph[v];
        const vertex_descriptor vRc = vertex.vRc;
        SHASTA2_ASSERT(vRc != null_vertex());
        SHASTA2_ASSERT(vRc != v);

        if(id(v) < id(vRc)) {
            candidateVertices.push_back(v);
        }
    }
    cout << "Found " << candidateVertices.size() << " candidate vertex pairs for detangling." << endl;



    // Now loop over out detangling candidates.
    vector<edge_descriptor> inEdges;        // in-edges of v
    vector<edge_descriptor> outEdges;       // out-edges of v
    vector<edge_descriptor> inEdgesRc;      // in-edges of vRc
    vector<edge_descriptor> outEdgesRc;     // out-edges of vRc
    for(const vertex_descriptor v: candidateVertices) {
        const AssemblyGraphVertex& vertex = assemblyGraph[v];
        const vertex_descriptor vRc = vertex.vRc;
        SHASTA2_ASSERT(vRc != null_vertex());
        SHASTA2_ASSERT(vRc != v);

        if(debug) {
            cout << "Working on vertex pair " << id(v) << " " << id(vRc) << endl;
        }

        // Gather incoming/outgoing edges.
        inEdges.clear();
        outEdges.clear();
        inEdgesRc.clear();
        outEdgesRc.clear();
        BGL_FORALL_INEDGES(v, e, assemblyGraph, AssemblyGraph) {
            inEdges.push_back(e);
        }
        BGL_FORALL_OUTEDGES(v, e, assemblyGraph, AssemblyGraph) {
            outEdges.push_back(e);
        }
        BGL_FORALL_INEDGES(vRc, eRc, assemblyGraph, AssemblyGraph) {
            inEdgesRc.push_back(eRc);
        }
        BGL_FORALL_OUTEDGES(vRc, eRc, assemblyGraph, AssemblyGraph) {
            outEdgesRc.push_back(eRc);
        }
        std::ranges::sort(inEdges, orderById);
        std::ranges::sort(outEdges, orderById);
        std::ranges::sort(inEdgesRc, orderById);
        std::ranges::sort(outEdgesRc, orderById);

        if(debug) {
            cout << "inEdges:";
            for(const edge_descriptor e: inEdges) {
                cout << " " << id(e);
            }
            cout << endl;
            cout << "outEdges:";
            for(const edge_descriptor e: outEdges) {
                cout << " " << id(e);
            }
            cout << endl;
            cout << "inEdgesRc:";
            for(const edge_descriptor e: inEdgesRc) {
                cout << " " << id(e);
            }
            cout << endl;
            cout << "outEdgesRc:";
            for(const edge_descriptor e: outEdgesRc) {
                cout << " " << id(e);
            }
            cout << endl;
        }

        // Create the tangle matrix.
        TangleMatrix tangleMatrix(assemblyGraph, inEdges, outEdges, html);
        if(debug) {
            cout << "Tangle  matrix: ";
            for(const auto& row: tangleMatrix.tangleMatrix) {
                std::ranges::copy(row, ostream_iterator<uint64_t>(cout, " "));
            }
            cout << endl;
        }

        // Run the G-test on this tangle matrix.
        GTest gTest(tangleMatrix.tangleMatrix, options.detangleEpsilon, false, false);

        // If the G-test failed, don't detangle.
        if(not gTest.success) {
            if(debug) {
                cout << "Likelihood ratio test was not successful." << endl;
            }
            continue;
        }

        const auto& bestHypothesis = gTest.hypotheses.front();
        const double bestG = bestHypothesis.G;
        if(debug) {
            cout << "Best hypothesis:";
            for(const auto& row: bestHypothesis.connectivityMatrix) {
                for(const bool b: row) {
                    cout << (b ? " 1" : " 0");
                }
            }
            cout << endl;
            cout << "G = " << bestG;
            if(gTest.hypotheses.size() > 1) {
                const double secondBestG = gTest.hypotheses[1].G;
                cout << ", second best G = " << secondBestG;
            }
            cout << endl;
        }

        // Check if the best hypothesis satisfies our options.
        if(bestG > options.detangleMaxLogP) {
            if(debug) {
                cout << "Best hypothesis G is too high." << endl;
            }
            continue;
        }
        if(gTest.hypotheses.size() > 1) {
            const double secondBestG = gTest.hypotheses[1].G;
            if(secondBestG - bestG < assemblyGraph.options.detangleMinLogPDelta) {
                if(debug) {
                    cout << "Second best hypothesis G is too low." << endl;
                }
                continue;
            }
        }

        // Figure out if in phase or out of phase.
        bool isInPhase = (bestHypothesis.connectivityMatrix == inPhaseConnectivityMatrix);
        bool isOutOfPhase = (bestHypothesis.connectivityMatrix == outOfPhaseConnectivityMatrix);
        if(not (isInPhase or isOutOfPhase)) {
            if(debug) {
                cout << "Not detangling because the best hypothesis is not a permutation." << endl;
            }
            continue;
        }
        if(debug) {
            if(isInPhase) {
                cout << "This vertex and its reverse complement will be detangled in phase." << endl;
            } else {
                SHASTA2_ASSERT(isOutOfPhase);
                cout << "This vertex and its reverse complement will be detangled out of phase." << endl;
            }
        }



        // If getting here, this vertex and its reverse complement will be detangled.
        // To do this, we need to disconnect (at the beginning or end) the Segments involved.
        // To do this we use disconnectAtBeginning and disconnectAtEnd, which change
        // edge_descriptors but leave the ids unchanged.
        // Some Segments can appear in more than one of the inEdges, outEdges,
        // inEdgesRc, outEdgesRc vectors. To avoid working with invalidated
        // edge_descriptors, we work with Segment ids instead.
        std::map<uint64_t, Segment> segmentMap;
        vector<uint64_t> inEdgeIds;
        for(const Segment e: inEdges) {
            inEdgeIds.push_back(id(e));
            segmentMap.insert(make_pair(id(e), e));
        }
        vector<uint64_t> outEdgeIds;
        for(const Segment e: outEdges) {
            outEdgeIds.push_back(id(e));
            segmentMap.insert(make_pair(id(e), e));
        }
        for(const Segment e: inEdgesRc) {
            segmentMap.insert(make_pair(id(e), e));
        }
        for(const Segment e: outEdgesRc) {
            segmentMap.insert(make_pair(id(e), e));
        }
        for(const uint64_t inEdgeId: inEdgeIds) {
            const Segment eNew = disconnectAtEnd(segmentMap.at(inEdgeId));
            segmentMap[id(eNew)] = eNew;
            const Segment eNewRc = assemblyGraph[eNew].eRc;
            segmentMap[id(eNewRc)] = eNewRc;
        }
        for(const uint64_t outEdgeId: outEdgeIds) {
            const Segment eNew = disconnectAtBeginning(segmentMap.at(outEdgeId));
            segmentMap[id(eNew)] = eNew;
            const Segment eNewRc = assemblyGraph[eNew].eRc;
            segmentMap[id(eNewRc)] = eNewRc;
        }

        // Check that all is good.
        for(const auto&[segmentId, e]: segmentMap) {
            SHASTA2_ASSERT(segmentId == id(e));
        }



        // Now we can just connect in-phase or out-of-phase inEdge/outEdge pairs.
        // Loop over the entrance/exit pairs to be connected.
        for(uint64_t i=0; i<2; i++) {
            const uint64_t entranceId = inEdgeIds[i];
            const uint64_t exitId = (isInPhase ? outEdgeIds[i] : outEdgeIds[1-i]);
            const Segment entrance = segmentMap.at(entranceId);
            const Segment exit = segmentMap.at(exitId);
            const vertex_descriptor vEntrance = target(entrance, assemblyGraph);
            const vertex_descriptor vExit = source(exit, assemblyGraph);

            // Connect them with a null edge.This will later be removed during compression.
            const auto[e, wasAdded] = add_edge(vEntrance, vExit, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
            createReverseComplementEdge(e);

        }

        // Now we can remove the vertex we detangled and its reverse complement.
        SHASTA2_ASSERT(in_degree(v, assemblyGraph) == 0);
        SHASTA2_ASSERT(out_degree(v, assemblyGraph) == 0);
        SHASTA2_ASSERT(in_degree(vRc, assemblyGraph) == 0);
        SHASTA2_ASSERT(out_degree(vRc, assemblyGraph) == 0);
        boost::remove_vertex(v, assemblyGraph);
        boost::remove_vertex(vRc, assemblyGraph);
    }

}



void AssemblyGraph::detangleEdges([[maybe_unused]] const string& debugOutputBaseName)
{
    AssemblyGraph& assemblyGraph = *this;
    ostream html(0);

    // Edges that are candidates for detangling define the tangles.
    vector< vector<vertex_descriptor> > tangles;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);
        if(out_degree(v0, assemblyGraph) != 1) {
            continue;
        }
        if(in_degree(v1, assemblyGraph) != 1) {
            continue;
        }
        if(in_degree(v0, assemblyGraph) != 2) {
            continue;
        }
        if(out_degree(v1, assemblyGraph) != 2) {
            continue;
        }
        tangles.push_back({v0, v1});
    }

    // Create a map that gives the tangle each vertex belongs to, if any.
    std::map<vertex_descriptor, uint64_t> tangleMap;
    for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
        const vector<vertex_descriptor>& tangle = tangles[tangleId];
        for(const vertex_descriptor v: tangle) {
            tangleMap.insert(make_pair(v, tangleId));
        }
    }

    // Find the reverse complement of each tangle.
    vector<uint64_t> tangleRc(tangles.size());
    for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
        const vector<vertex_descriptor>& tangle = tangles[tangleId];
        const vertex_descriptor v = tangle.front();
        const vertex_descriptor vRc = assemblyGraph[v].vRc;
        tangleRc[tangleId] = tangleMap.at(vRc);
    }

#if 0
    // Detangling, no read following.
    const bool attemptReadFollowing = false;
    detangleAndReadFollowing(tangles, tangleRc, attemptReadFollowing, debugOutputBaseName);
#endif


    // Detangle each pair.
    for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
        if(tangleId < tangleRc[tangleId]) {
            const Tangle tangle(assemblyGraph, tangles[tangleId]);
            detangleTanglePairStrict(tangle, html);
        }
    }

}



bool AssemblyGraph::detangleAndReadFollowingSuperbubbles(const string& debugOutputBaseName)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t lengthThreshold = 100000;

    performanceLog << timestamp << "AssemblyGraph::detangleAndReadFollowingSuperbubbles begins: " <<
        debugOutputBaseName << endl;

    // Create the superbubbles as tangles consisting of short Segments.
    vector< vector<vertex_descriptor> > tangles;
    vector<uint64_t> tangleRc;
    createTanglesBySegmentLength(lengthThreshold, tangles, tangleRc);

    // Write a csv file that can be imported into Bandage to see
    // the tangles.
    writeTangles(tangles, tangleRc, debugOutputBaseName + "-Tangles-Bandage.csv");

    // Detangling and read following.
    const bool attemptReadFollowing = true;
    return detangleAndReadFollowing(tangles, tangleRc, attemptReadFollowing, debugOutputBaseName);

    performanceLog << timestamp << "AssemblyGraph::detangleAndReadFollowingSuperbubbles ends: " <<
        debugOutputBaseName << endl;
}



bool AssemblyGraph::detangleAndReadFollowing(
    const vector< vector<vertex_descriptor> >& tangles,
    const vector<uint64_t>& tangleRc,
    bool attemptReadFollowing,
    const string& debugOutputBaseName)
{
    AssemblyGraph& assemblyGraph = *this;
    const bool debug = false;

    // Detangle or read following on each pair of reverse complemented tangles.
    bool somethingWasDone = false;
    for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
        if(tangleId <= tangleRc[tangleId]) {
            const Tangle tangle(assemblyGraph, tangles[tangleId]);

            ofstream html;
            if(debug) {

                cout << "Working on tangle " << tangleId <<
                    " and its reverse complement tangle " << tangleRc[tangleId] << endl;
                cout << "This tangle has " <<
                    tangle.entrances.size() << " entrances, " <<
                    tangle.exits.size() << " exits, " <<
                    tangle.tangleEdges.size() << " segments, " <<
                    tangle.entrances.size() << " vertices." << endl;

            }

            if(tangle.entrances.empty() or tangle.exits.empty()) {
                if(debug) {
                    cout << "Skipping tangle " << tangleId << " with no entrances and/or no exits." << endl;
                }
                continue;
            }

            if(debug) {
                // Tangle details go to the html file.
                html.open(debugOutputBaseName + "-Tangle-" + to_string(tangleId) + ".html");
                writeHtmlBegin(html, "Tangle " + to_string(tangleId));
                html << "<h1>Tangle " << tangleId << "</h1>";
                tangle.writeHtml(html);
            }

            // Try detangling.
            // If detangling was successful, we are done with this tangle.
            if(detangleStrandSymmetric(tangle, html)) {
                if(debug) {
                    cout << "Detangling was successful on this tangle." << endl;
                }
                somethingWasDone = true;
                continue;
            }

            // Detangling was not successful for this tangle.
            // Try read following, if requested.
            if(debug) {
                cout << "Detangling was not successful on this tangle." << endl;
            }
            if(attemptReadFollowing) {
                if(debug) {
                    cout << "Attempting read following." << endl;
                }
                if(readFollowingStrandSymmetric(tangleId, tangle, html)) {
                    if(debug) {
                        cout << "Read following was successful on this tangle." << endl;
                    }
                    somethingWasDone = true;
                } else {
                    if(debug) {
                        cout << "Read following was not successful on this tangle." << endl;
                    }
                }
            }
        }
    }

    return somethingWasDone;
}



// Write a csv file that can be imported into Bandage to see
// the tangles.
void AssemblyGraph::writeTangles(
    const vector< vector<vertex_descriptor> >& tangles,
    const vector<uint64_t>& tangleRc,
    const string& fileName) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream csv(fileName);
    csv << "Segment,Tangle,TangleRc,Color\n";

    for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
        const string color = randomHslColor(tangleId, 0.75, 0.5);
        const vector<vertex_descriptor>& tangle = tangles[tangleId];

        for(const vertex_descriptor v0: tangle) {

            BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
                const vertex_descriptor v1 = target(e, assemblyGraph);
                if(binary_search(tangle.begin(), tangle.end(), v1, orderById)) {
                    csv << id(e) << ",";
                    csv << tangleId << ",";
                    csv << tangleRc[tangleId] << ",";
                    csv << color << "\n";
                }
            }
        }
    }
}



bool AssemblyGraph::detangleStrandSymmetric(const Tangle& tangle, ostream& html)
{
    if(
        (tangle.isSelfComplementary()) and
        (tangle.entrances.size() == 2) and
        (tangle.exits.size() == 2)) {
        return detangleSelfComplementaryTangle2By2(tangle, html);
    } else {
        return detangleTanglePair(tangle, html);
    }
}



// This detangles the tangle passed in as an argument.
// If the tangle is not self-complementary,
// it also detangles its reverse complement.
bool AssemblyGraph::detangleTanglePair(
    const Tangle& tangle, ostream& html)
{
    AssemblyGraph& assemblyGraph = *this;
    SHASTA2_ASSERT(not tangle.entrances.empty());
    SHASTA2_ASSERT(not tangle.exits.empty());

    if(html) {
        html << "<h2>Detangling</h2>";
    }

    // Compute the block structure of the TangleMatrix.
    vector<TangleMatrix> blocks;
    tangle.tangleMatrix().findBlockStructure(blocks);



    // Run a G-test for each block.
    // Gather entrance/exit pairs to be connected.
    // The bool of each pair says whether that pair should be
    // connected using a "deep" call to connect.
    vector< pair< pair<uint64_t, uint64_t>, bool> > connectPairs;
    vector<GTest> gTests;
    bool connectionFailure = false;
    vector<uint64_t> failedBlocks;
    for(uint64_t iBlock=0; iBlock<blocks.size(); iBlock++) {
        const TangleMatrix& block = blocks[iBlock];

        const bool onlyConsiderPermutation = (block.entrances.size() == block.exits.size());
        const bool onlyConsiderInjective = onlyConsiderPermutation;
        const GTest& gTest = gTests.emplace_back(block.tangleMatrix, assemblyGraph.options.detangleEpsilon,
        onlyConsiderInjective, onlyConsiderPermutation);

        if(gTest.isPositive(assemblyGraph.options.detangleMaxLogP, assemblyGraph.options.detangleMinLogPDelta)) {

            // The likelihood ratio test gave a reliable result.
            // Find the pairs to be connected.
             for(uint64_t i=0; i<block.entrances.size(); i++) {
                const Segment entrance = block.entrances[i];
                for(uint64_t j=0; j<block.exits.size(); j++) {
                    const Segment exit = block.exits[j];
                    if(gTest.hypotheses[0].connectivityMatrix[i][j]) {
                        if(not canConnect(entrance, exit)) {
                            connectionFailure = true;
                            if(html) {
                                html << "<br>Cannot connect " << id(entrance) <<
                                    " " << id(exit);
                            }
                        }
                        connectPairs.push_back({{id(entrance), id(exit)}, false});
                    }
                }
            }
        } else {
            failedBlocks.push_back(iBlock);
        }
    }



    if(html) {
        for(uint64_t iBlock=0; iBlock<blocks.size(); iBlock++) {
            const TangleMatrix& block = blocks[iBlock];
            const GTest gTest = gTests[iBlock];

            html << "<h3>Tangle matrix block " << iBlock << "</h3>";
            block.writeTotalTangleMatrix(html);

            if(gTest.success) {
                const auto& bestHypothesis = gTest.hypotheses.front();
                const double bestG = bestHypothesis.G;
                html <<
                    "<br>G-test best hypothesis:"
                    "<table><tr><th>";
                for(const AssemblyGraph::edge_descriptor exit: block.exits) {
                    html << "<th>" << id(exit);
                }
                for(uint64_t i=0; i<block.entrances.size(); i++) {
                    html << "<tr><th>" << id(block.entrances[i]);
                    for(uint64_t j=0; j<block.exits.size(); j++) {
                        html << "<td class=centered>" << int(bestHypothesis.connectivityMatrix[i][j]);
                    }
                }
                html << "</table><br>G = " << bestG;
                if(gTest.hypotheses.size() > 1) {
                    const double secondBestG = gTest.hypotheses[1].G;
                    html << "<br>Second best G = " << secondBestG;
                    html << "<br>&Delta;G = " << secondBestG - bestG;
                }

                if(gTest.isPositive(assemblyGraph.options.detangleMaxLogP, assemblyGraph.options.detangleMinLogPDelta)) {
                    html << "<br>The G-test indicates that the best hypothesis for this block is reliable.";
                } else {
                    html << "<br>The G-test indicates that the best hypothesis for this block is not reliable.";
               }

            } else {
                html << "<br>The G-test for this block failed.";
            }
        }
    }

    // Don't tolerate any connection failures.
    if(connectionFailure) {
        if(html) {
            html << "<br>Not detangling because of connection failures.";
        }
        return false;
    }



    // In general, we don't want any failed blocks.
    // But we make an exception to handle an important special case:
    // - Only two block failed.
    // - One of the two failed blocks consists of just one entrance.
    // - One of the two failed blocks consists of just one exit.
    // - That entrance/exit pair can be connected using a "deep"
    //   RestrictedAnchorGraph.
    // In that case, we will also connect this entrance/exit pair,
    // making sure to use a deep RestrictedAnchorGraph.
    bool isSpecialCase = false;
    Segment specialCaseEntrance = assemblyGraphNullEdge;
    Segment specialCaseExit = assemblyGraphNullEdge;
    if(failedBlocks.size() == 2) {
        const TangleMatrix& block0 = blocks[failedBlocks[0]];
        const TangleMatrix& block1 = blocks[failedBlocks[1]];

        if(
            (block0.entrances.size() == 1) and
            (block0.exits.empty()) and
            (block1.entrances.empty()) and
            (block1.exits.size() == 1))
        {
            specialCaseEntrance = block0.entrances.front();
            specialCaseExit = block1.exits.front();
        }

        else if(
            (block1.entrances.size() == 1) and
            (block1.exits.empty()) and
            (block0.entrances.empty()) and
            (block0.exits.size() == 1))
        {
            specialCaseEntrance = block1.entrances.front();
            specialCaseExit = block0.exits.front();
        }

        // We still need to check if this entrance/exit pair
        // can be connected using a deep RestrictedAnchorGraph.
        if(specialCaseEntrance != assemblyGraphNullEdge) {
            if(canConnect(specialCaseEntrance, specialCaseExit, true)) {
                isSpecialCase = true;
                if(html) {
                    html << "<br>Special case: " <<
                        id(specialCaseEntrance) << " and " << id(specialCaseExit) <<
                        " can be connected using a deep RestrictedAnchorGraph." << endl;
                }
            } else {
                if(html) {
                    html << "<br>Special case: " <<
                        id(specialCaseEntrance) << " and " << id(specialCaseExit) <<
                        " cannot be connected using a deep RestrictedAnchorGraph." << endl;
                }
            }
        }
    }
    uint64_t specialCaseEntranceId = invalid<uint64_t>;
    uint64_t specialCaseExitId = invalid<uint64_t>;
    if(isSpecialCase) {
        specialCaseEntranceId = id(specialCaseEntrance);
        specialCaseExitId = id(specialCaseExit);
    }




    // If there are failed blocks and we are not in that special case, don't detangle.
    if((not failedBlocks.empty()) and (not isSpecialCase)) {
        if(html) {
            html << "<br>Not detangling because " << failedBlocks.size() <<
                " blocks out of " << blocks.size() <<
                " did not give a sufficiently reliable result.";
        }
        return false;
    }

    // If getting here, we will detangle.
    if(html) {
        html << "<br>This tangle and its reverse complement will be detangled." << endl;
    }

    // For the special case we have to use a "deep" connect.
    if(isSpecialCase) {
        connectPairs.push_back({{specialCaseEntranceId, specialCaseExitId}, true});
    }

    // Make the connections.
    detangleMakeConnections(tangle, connectPairs);

#if 0
    // Remove all the Segments internal to this tangle
    // and their reverse complements.
    for(const Segment e: tangle.tangleEdges) {
        if(not isSelfComplementary) {
            boost::remove_edge(assemblyGraph[e].eRc, assemblyGraph);
        }
        boost::remove_edge(e, assemblyGraph);
    }



    // If getting here, this tangle and its reverse complement will be detangled.
    // To do this, we need to disconnect (at the beginning or end) the Segments involved.
    // To do this we use disconnectAtBeginning and disconnectAtEnd, which change
    // edge_descriptors but leave the ids unchanged.
    // Some Segments can at the same time be an entrance and/or exit
    // of this tangle or its reverse complement. To avoid working with invalidated
    // edge_descriptors, we work with Segment ids instead.
    std::map<uint64_t, Segment> segmentMap;
    tangle.createSegmentMap(segmentMap);
    vector<uint64_t> entranceIds;
    for(const Segment e: tangle.entrances) {
        entranceIds.push_back(id(e));
    }
    vector<uint64_t> exitIds;
    for(const Segment e: tangle.exits) {
        exitIds.push_back(id(e));
    }

    // Create disconnected versions of the entrances and exits,
    // while keeping the segmentMap up to date.
    for(const uint64_t entranceId: entranceIds) {
        const Segment eNew = disconnectAtEnd(segmentMap.at(entranceId));
        segmentMap[id(eNew)] = eNew;
        const Segment eNewRc = assemblyGraph[eNew].eRc;
        segmentMap[id(eNewRc)] = eNewRc;
    }
    if(not isSelfComplementary) {
        for(const uint64_t exitId: exitIds) {
            const Segment eNew = disconnectAtBeginning(segmentMap.at(exitId));
            segmentMap[id(eNew)] = eNew;
            const Segment eNewRc = assemblyGraph[eNew].eRc;
            segmentMap[id(eNewRc)] = eNewRc;
        }
    }
    tangle.checkSegmentMap(segmentMap);



    // Make the connections.
    // This also does the same for the reverse complement of this Tangle.
    std::set<uint64_t> connectedEntrances;
    for(const auto&[connectPair, deep]: connectPairs) {
        const uint64_t entranceId = connectPair.first;
        const uint64_t exitId = connectPair.second;
        if(isSelfComplementary and (connectedEntrances.contains(entranceId))) {
            continue;
        }
        const AssemblyGraph::edge_descriptor entrance = segmentMap.at(entranceId);
        const AssemblyGraph::edge_descriptor exit = segmentMap.at(exitId);
        const edge_descriptor eNew = connect(entrance, exit, deep);
        createReverseComplementEdge(eNew);
        if(isSelfComplementary) {
            connectedEntrances.insert(entranceId);
            connectedEntrances.insert(id(assemblyGraph[exit].eRc));
        }
    }



    // Remove all the Tangle vertices and their reverse complements
    // that are now isolated. These are the ones that were
    // not connected to any entrance or exit.
    for(const vertex_descriptor v: tangle.tangleVertices) {
        const vertex_descriptor vRc = assemblyGraph[v].vRc;

        if((in_degree(v, assemblyGraph) == 0) and (out_degree(v, assemblyGraph) == 0)) {
            boost::remove_vertex(v, assemblyGraph);
        }

        if(not isSelfComplementary) {
            if((in_degree(vRc, assemblyGraph) == 0) and (out_degree(vRc, assemblyGraph) == 0)) {
                boost::remove_vertex(vRc, assemblyGraph);
            }
        }
    }
#endif

    return true;
}


// Make connections for detangling a detangle.
// The connectPairs contain pairs of entrance/exit segmentIds,
// and the bool parameter says whether a "deep" connect should be used for that pair.
void AssemblyGraph::detangleMakeConnections(
    const Tangle& tangle,
    const vector< pair< pair<uint64_t, uint64_t>, bool> >& connectPairs
    )
{
    AssemblyGraph& assemblyGraph = *this;
    const bool isSelfComplementary = tangle.isSelfComplementary();

    // Remove all the Segments internal to this tangle
    // and their reverse complements.
    for(const Segment e: tangle.tangleEdges) {
        if(not isSelfComplementary) {
            boost::remove_edge(assemblyGraph[e].eRc, assemblyGraph);
        }
        boost::remove_edge(e, assemblyGraph);
    }



    // If getting here, this tangle and its reverse complement will be detangled.
    // To do this, we need to disconnect (at the beginning or end) the Segments involved.
    // To do this we use disconnectAtBeginning and disconnectAtEnd, which change
    // edge_descriptors but leave the ids unchanged.
    // Some Segments can at the same time be an entrance and/or exit
    // of this tangle or its reverse complement. To avoid working with invalidated
    // edge_descriptors, we work with Segment ids instead.
    std::map<uint64_t, Segment> segmentMap;
    tangle.createSegmentMap(segmentMap);
    vector<uint64_t> entranceIds;
    for(const Segment e: tangle.entrances) {
        entranceIds.push_back(id(e));
    }
    vector<uint64_t> exitIds;
    for(const Segment e: tangle.exits) {
        exitIds.push_back(id(e));
    }

    // Create disconnected versions of the entrances and exits,
    // while keeping the segmentMap up to date.
    for(const uint64_t entranceId: entranceIds) {
        const Segment eNew = disconnectAtEnd(segmentMap.at(entranceId));
        segmentMap[id(eNew)] = eNew;
        const Segment eNewRc = assemblyGraph[eNew].eRc;
        segmentMap[id(eNewRc)] = eNewRc;
    }
    if(not isSelfComplementary) {
        for(const uint64_t exitId: exitIds) {
            const Segment eNew = disconnectAtBeginning(segmentMap.at(exitId));
            segmentMap[id(eNew)] = eNew;
            const Segment eNewRc = assemblyGraph[eNew].eRc;
            segmentMap[id(eNewRc)] = eNewRc;
        }
    }
    tangle.checkSegmentMap(segmentMap);



    // Make the connections.
    // This also does the same for the reverse complement of this Tangle.
    std::set<uint64_t> connectedEntrances;
    for(const auto&[connectPair, deep]: connectPairs) {
        const uint64_t entranceId = connectPair.first;
        const uint64_t exitId = connectPair.second;
        if(isSelfComplementary and (connectedEntrances.contains(entranceId))) {
            continue;
        }
        const AssemblyGraph::edge_descriptor entrance = segmentMap.at(entranceId);
        const AssemblyGraph::edge_descriptor exit = segmentMap.at(exitId);
        const edge_descriptor eNew = connect(entrance, exit, deep);
        createReverseComplementEdge(eNew);
        if(isSelfComplementary) {
            connectedEntrances.insert(entranceId);
            connectedEntrances.insert(id(assemblyGraph[exit].eRc));
        }
    }



    // Remove all the Tangle vertices and their reverse complements
    // that are now isolated. These are the ones that were
    // not connected to any entrance or exit.
    for(const vertex_descriptor v: tangle.tangleVertices) {
        const vertex_descriptor vRc = assemblyGraph[v].vRc;

        if((in_degree(v, assemblyGraph) == 0) and (out_degree(v, assemblyGraph) == 0)) {
            boost::remove_vertex(v, assemblyGraph);
        }

        if(not isSelfComplementary) {
            if((in_degree(vRc, assemblyGraph) == 0) and (out_degree(vRc, assemblyGraph) == 0)) {
                boost::remove_vertex(vRc, assemblyGraph);
            }
        }
    }

}



// This detangles the tangle passed in as an argument.
// If the tangle is not self-complementary,
// it also detangles its reverse complement.
// This "strict" version:
// - Does not use the block structure of the tangle matrix.
// - It allows all hypotheses to be considered in the
//   likelihood ratio test, but requires the top hypothesis
//   to be a permutation.
// - It requires the tangle not to be self-complementary.
bool AssemblyGraph::detangleTanglePairStrict(
    const Tangle& tangle, ostream& html)
{
    AssemblyGraph& assemblyGraph = *this;
    SHASTA2_ASSERT(not tangle.entrances.empty());
    SHASTA2_ASSERT(not tangle.exits.empty());

    if(tangle.entrances.size() != tangle.exits.size()) {
        return false;
    }

    if(tangle.isSelfComplementary()) {
        return false;
    }

    if(html) {
        html << "<h2>Detangling</h2>";
    }

    // Run the likelihood ratio test on the entire matrix,
    // without using its block structure. Allow all hypotheses.
    const bool onlyConsiderPermutation = false;
    const bool onlyConsiderInjective = false;
    const GTest gTest(tangle.tangleMatrix().tangleMatrix, assemblyGraph.options.detangleEpsilon,
        onlyConsiderInjective, onlyConsiderPermutation);

    // If the likelihood ratio test failed, don't do anything.
    if(not gTest.success) {
        return false;
    }

    // If the top hypothesis is not a permutation, don't do anything.
    SHASTA2_ASSERT(not gTest.hypotheses.empty());
    const GTest::Hypothesis& topHypothesis = gTest.hypotheses.front();
    if(not GTest::isForwardInjective(topHypothesis.connectivityMatrix)) {
        return false;
    }
    if(not GTest::isBackwardInjective(topHypothesis.connectivityMatrix)) {
        return false;
    }

    // If the likelihood ratio test is not sufficiently reliable,
    // don't do anything.
    if(not gTest.isPositive(
        assemblyGraph.options.detangleMaxLogP,
        assemblyGraph.options.detangleMinLogPDelta)) {
        return false;
    }

    // Create the connect pairs.
    vector< pair<pair<uint64_t, uint64_t>, bool> > connectPairs;
    for(uint64_t i=0; i<tangle.entrances.size(); i++) {
        for(uint64_t j=0; j<tangle.exits.size(); j++) {
            if(topHypothesis.connectivityMatrix[i][j]) {
                if(not canConnect(tangle.entrances[i], tangle.exits[j])) {
                    return false;
                }
                connectPairs.push_back(
                    {{id(tangle.entrances[i]), id(tangle.exits[j])}, false});
            }
        }
    }

    // Make the connections.
    detangleMakeConnections(tangle, connectPairs);
    return true;
}



bool AssemblyGraph::detangleSelfComplementaryTangle2By2(
    const Tangle& tangle,
    ostream& html)
{
    AssemblyGraph& assemblyGraph = *this;
    if(html) {
        html << "<br>AssemblyGraph::detangleSelfComplementaryTangle2By2 begins" << endl;
    }

    const vector<Segment>& entrances = tangle.entrances;
    const vector<Segment>& exits = tangle.exits;

    SHASTA2_ASSERT(tangle.isSelfComplementary());
    SHASTA2_ASSERT(entrances.size() == 2);
    SHASTA2_ASSERT(exits.size() == 2);

    // In a self-complementary tangle, each exit is the reverse complement
    // of an entrance. We can't connect an entrance with its own reverse
    // complement, so each entrance can only possibly be connected to
    // the one exit with is not its reverse complement.
    // Therefore we don't need to compute a TangleMatrix and run a G-test.
    // We just need to check that that connection is possible.

    // Store the reverse complement of each entrance.
    vector<Segment> entrancesRc(2);
    for(uint64_t i=0; i<2; i++) {
        entrancesRc[i] = assemblyGraph[entrances[i]].eRc;
        SHASTA2_ASSERT(std::ranges::contains(exits, entrancesRc[i]));
    }

    // Define the entrance/exit pair that we will explictly connect.
    // The remaining entrance/exit pair will be connected
    // automatically.
    const Segment entrance = entrances[0];
    const Segment exit = entrancesRc[1];
    const uint64_t entranceId = id(entrance);
    const uint64_t exitId = id(exit);
    if(html) {
        html << "<br>Will attempt to connect " << entranceId << " with " << exitId << endl;
    }

    // Check if  we can connect this entrance/exit pair.
    if(not assemblyGraph.canConnect(entrance, exit)) {
        if(html) {
            html << "<br>Not detangling because can't connect " << entranceId <<
                " with " << exitId << endl;
        }
        return false;
    }

    if(html) {
        html << "<br>This tangle will be detangled." << endl;
    }

    // Remove all the Segments internal to this tangle.
    for(const Segment e: tangle.tangleEdges) {
        boost::remove_edge(e, assemblyGraph);
    }



    // If getting here, this tangle will be detangled.
    // To do this, we need to disconnect (at the beginning or end) the Segments involved.
    // To do this we use disconnectAtBeginning and disconnectAtEnd, which change
    // edge_descriptors but leave the ids unchanged.
    // Some Segments can at the same time be an entrance and/or exit
    // of this tangle or its reverse complement. To avoid working with invalidated
    // edge_descriptors, we work with Segment ids instead.
    std::map<uint64_t, Segment> segmentMap;
    tangle.createSegmentMap(segmentMap);
    vector<uint64_t> entranceIds;
    for(const Segment e: tangle.entrances) {
        entranceIds.push_back(id(e));
    }
    vector<uint64_t> exitIds;
    for(const Segment e: tangle.exits) {
        exitIds.push_back(id(e));
    }

    // Disconnect the entrances and exits,
    // while keeping the segmentMap up to date.
    for(const uint64_t entranceId: entranceIds) {
        const Segment eNew = disconnectAtEnd(segmentMap.at(entranceId));
        segmentMap[id(eNew)] = eNew;
        const Segment eNewRc = assemblyGraph[eNew].eRc;
        segmentMap[id(eNewRc)] = eNewRc;
    }
    for(const uint64_t exitId: exitIds) {
        const Segment eNew = disconnectAtBeginning(segmentMap.at(exitId));
        segmentMap[id(eNew)] = eNew;
        const Segment eNewRc = assemblyGraph[eNew].eRc;
        segmentMap[id(eNewRc)] = eNewRc;
    }
    tangle.checkSegmentMap(segmentMap);



    // Make the connection and its reverse complement.
    const edge_descriptor eNew = connect(segmentMap.at(entranceId), segmentMap.at(exitId));
    createReverseComplementEdge(eNew);

    // Remove all the Tangle vertices that are now isolated.
    // These are the ones that were not connected to any entrance or exit.
    for(const vertex_descriptor v: tangle.tangleVertices) {
        if((in_degree(v, assemblyGraph) == 0) and (out_degree(v, assemblyGraph) == 0)) {
            boost::remove_vertex(v, assemblyGraph);
        }
    }

    return true;
}



// Create tangles as connected components when only short Segments are considered.
// If i the index of a tangle (in the tangles vector), tangleRc[i] gives the
// index of its reverse complement.
// If tangleRc[i] == i, the tangle with index i is self-complementary.
void AssemblyGraph::createTanglesBySegmentLength(
    uint64_t maxLength,
    vector< vector<vertex_descriptor> >& tangles,
    vector<uint64_t>& tangleRc
    ) const
{
    const AssemblyGraph& assemblyGraph = *this;

    // Map the vertices to integers.
    // This is needed below to compute connected components.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    std::vector<vertex_descriptor> vertexTable;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexIndexMap.insert(make_pair(v, vertexTable.size()));
        vertexTable.push_back(v);
    }

    // Compute connected components using only short edges.
    // Each non-trivial connected component will become a tangle.
    DisjointSets disjointSets(vertexTable.size());
    BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {
        if(assemblyGraph[segment].length() <= maxLength) {
            const vertex_descriptor v0 = source(segment, assemblyGraph);
            const vertex_descriptor v1 = target(segment, assemblyGraph);
            disjointSets.unionSet(vertexIndexMap.at(v0), vertexIndexMap.at(v1));
        }
    }


    // Get the connected components.
    // These will be our tangles.
    vector< vector<uint64_t> > componentsVertexIndexes;
    disjointSets.gatherComponents(1, componentsVertexIndexes);

    // Convert vertex indexes to vertex descriptors.
    // Sort each component so we can do binary searches.
    tangles.clear();
    for(const auto& componentVertexIndexes: componentsVertexIndexes) {
        vector<vertex_descriptor>& tangle = tangles.emplace_back();
        for(const uint64_t vertexIndex: componentVertexIndexes) {
            tangle.push_back(vertexTable[vertexIndex]);
        }
        std::ranges::sort(tangle, orderById);
    }

    // Create a map that gives the tangle each vertex belongs to.
    std::map<vertex_descriptor, uint64_t> tangleMap;
    for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
        const vector<vertex_descriptor>& tangle = tangles[tangleId];
        for(const vertex_descriptor v: tangle) {
            tangleMap.insert(make_pair(v, tangleId));
        }
    }

    // Find the reverse complement of each tangle.
    tangleRc.resize(tangles.size());
    for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
        const vector<vertex_descriptor>& tangle = tangles[tangleId];
        const vertex_descriptor v = tangle.front();
        const vertex_descriptor vRc = assemblyGraph[v].vRc;
        tangleRc[tangleId] = tangleMap.at(vRc);
    }

}

