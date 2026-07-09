// Shasta2.
#include "AssemblyGraph.hpp"
#include "Anchor.hpp"
#include "areSimilarSequences.hpp"
#include "deduplicate.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>



// Find Bubbles.
// The edges of each Bubble are sorted by id.
void AssemblyGraph::findBubbles(vector<Bubble>& bubbles) const
{
    performanceLog << timestamp << "AssemblyGraph::findBubbles begins." << endl;

    const AssemblyGraph& assemblyGraph = *this;
    bubbles.clear();

    // Look at bubbles with source v0.
    std::map<vertex_descriptor, vector<edge_descriptor> > m;
    BGL_FORALL_VERTICES(v0, assemblyGraph, AssemblyGraph) {
        m.clear();
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            m[v1].push_back(e);
        }

        for(const auto& p: m) {
            const vertex_descriptor v1 = p.first;
            const vector<edge_descriptor>& edges = p.second;
            if(edges.size() > 1) {
                Bubble bubble;
                bubble.v0 = v0;
                bubble.v1 = v1;
                bubble.edges = edges;
                std::ranges::sort(bubble.edges, orderById);
                bubbles.push_back(bubble);
            }
        }
    }

#if 0
    cout << "Found " << bubbles.size() << " bubbles." << endl;

    // Count the bubbles for each ploidy.
    vector<uint64_t> histogram;
    for(const Bubble& bubble: bubbles) {
        const uint64_t ploidy = bubble.edges.size();
        if(ploidy >= histogram.size()) {
            histogram.resize(ploidy+1, 0);
        }
        ++histogram[ploidy];
    }
    for(uint64_t ploidy=0; ploidy<histogram.size(); ploidy++) {
        const uint64_t frequency = histogram[ploidy];
        if(frequency) {
            cout << "Ploidy " << ploidy << ": " << frequency << " bubbles." << endl;
        }
    }
#endif

    performanceLog << timestamp << "AssemblyGraph::findBubbles ends." << endl;
}



uint64_t AssemblyGraph::bubbleCleanup()
{
    uint64_t modifiedCount = 0;

    vector< pair<vertex_descriptor, vertex_descriptor> > excludeList;
    while(true) {
        const uint64_t modifiedCountThisIteration = bubbleCleanupIterationMultithreaded(excludeList, options.threadCount);
        if(modifiedCountThisIteration == 0) {
            break;
        }
        modifiedCount += modifiedCountThisIteration;
    }

    return modifiedCount;
}



uint64_t AssemblyGraph::bubbleCleanupIterationMultithreaded(
    vector< pair<vertex_descriptor, vertex_descriptor> >& excludeList,
    uint64_t threadCount)
{
    performanceLog << timestamp << "Bubble cleanup iteration begins." << endl;
    AssemblyGraph& assemblyGraph = *this;

    // Find all bubbles.
    vector<Bubble> allBubbles;
    findBubbles(allBubbles);
    cout << "Found " << allBubbles.size() << " bubbles." << endl;

    // Find candidate bubbles.
    // These are the ones that are not in the exclude list and
    // in which no branch has offset greater than
    // options.bubbleCleanupMaxBubbleLength.
    vector<Bubble>& candidateBubbles = bubbleCleanupIterationData.candidateBubbles;
    candidateBubbles.clear();
    for(const Bubble& bubble: allBubbles) {

        if(std::ranges::contains(excludeList, make_pair(bubble.v0, bubble.v1))) {
            continue;
        }

        bool hasLongBranch = false;
        for(const edge_descriptor e: bubble.edges) {
            if(assemblyGraph[e].offset() > options.bubbleCleanupMaxBubbleLength) {
                hasLongBranch = true;
                break;
            }
        }

        if(not hasLongBranch) {
            candidateBubbles.push_back(bubble);
        }
    }
    // cout << candidateBubbles.size() << " bubbles are candidate for clean up." << endl;

    // Assemble sequence for all the edges of these bubbles.
    edgesToBeAssembled.clear();
    for(const Bubble& bubble: candidateBubbles) {
        for(const edge_descriptor e: bubble.edges) {
            if(not assemblyGraph[e].wasAssembled) {
                edgesToBeAssembled.push_back(e);
            }
        }
    }
    assemble();


    // Process the bubbles in multithreaded code.
    // All changes to the AssemblyGraph are done under mutex.
    bubbleCleanupIterationData.modifiedCount = 0;
    const uint64_t batchSize = 10;
    setupLoadBalancing(candidateBubbles.size(), batchSize);
    runThreads(&AssemblyGraph::bubbleCleanupIterationThreadFunction, threadCount);
    // cout << "Bubble cleanup modified " << modifiedCount << " bubbles." << endl;

    // Update the excludeList.
    for(const Bubble& bubble: candidateBubbles) {
        excludeList.push_back(make_pair(bubble.v0, bubble.v1));
    }
    std::ranges::sort(excludeList);

    performanceLog << timestamp << "Bubble cleanup iteration ends." << endl;
    return bubbleCleanupIterationData.modifiedCount;

}



void AssemblyGraph::bubbleCleanupIterationThreadFunction([[maybe_unused]] uint64_t threadId)
{
    vector<Bubble>& candidateBubbles = bubbleCleanupIterationData.candidateBubbles;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over candidate bubbles in this batch.
        for(uint64_t i=begin; i!=end; ++i) {
            if(bubbleCleanup(candidateBubbles[i])) {
                __sync_fetch_and_add(&bubbleCleanupIterationData.modifiedCount, 1);
            }
        }
    }
}






bool AssemblyGraph::bubbleCleanup(const Bubble& bubble)
{
    // EXPOSE WHEN CODE STABILIZES.
    const vector<uint64_t> minRepeatCount = {0, 2, 2, 2, 2, 2, 2};

    AssemblyGraph& assemblyGraph = *this;
    const uint64_t ploidy = bubble.edges.size();

    const bool debug = false; // (assemblyGraph[bubble.edges.front()].id)== 130581;
    if(debug) {
        cout << "Attempting cleanup for bubble";
        for(uint64_t i=0; i<ploidy; i++) {
            const edge_descriptor e = bubble.edges[i];
            cout << " (" << i << " " << assemblyGraph[e].id << ")";
        }
        cout << endl;
    }

    // Find which pairs of branches in the bubble have similar sequence
    // as defined by the minRepeatCount vector.
    // See analyzeBubble for more information.
    vector< pair<uint64_t, uint64_t> > similarPairs;
    analyzeBubble(bubble, minRepeatCount, similarPairs);

    // If no similar pairs were found, leave this bubble alone.
    if(similarPairs.empty()) {
        if(debug) {
            cout << "No similar branch sequences were found. Nothing done for this bubble." << endl;
        }
        return false;
    }

    if(debug) {
        cout << "These pairs of branches have similar sequence:" << endl;
        for(const auto& p: similarPairs) {
            const uint64_t i0 = p.first;
            const uint64_t i1 = p.second;
            cout << i0 << " " << i1 << " " << assemblyGraph[bubble.edges[i0]].id <<
                " " << assemblyGraph[bubble.edges[i1]].id << endl;
        }
    }

    // Find the OrientedReadIds that appear in all the steps of each branch of the bubble.
    vector< vector<OrientedReadId> > allOrientedReadIds(ploidy);
    for(uint64_t i=0; i<ploidy; i++) {
        const AssemblyGraphEdge& edge = assemblyGraph[bubble.edges[i]];
        for(const auto& step: edge) {
            std::ranges::copy(step.anchorPair.orientedReadIds, back_inserter(allOrientedReadIds[i]));
        }
        deduplicate(allOrientedReadIds[i]);
        if(debug) {
            cout << "Branch " << i << " has " << allOrientedReadIds[i].size() <<
                " total oriented read ids." << endl;
        }
    }

#if 0
    // Gather all OrientedReadIds that appear in exactly one of the branches.
    vector<OrientedReadId> unambiguousOrientedReadIdsAllBranches;
    for(const auto& v: allOrientedReadIds) {
        std::ranges::copy(v, back_inserter(unambiguousOrientedReadIdsAllBranches));
    }
    deduplicateAndCountAndKeepUnique(unambiguousOrientedReadIdsAllBranches);

    // The unambiguous OrientedReadIds for each branch are the intersection
    // of unambiguousOrientedReadIdsAllBranches with allOrientedReadIds for that branch.
    vector< vector<OrientedReadId> > unambiguousOrientedReadIds(bubble.edges.size());
    for(uint64_t i=0; i<ploidy; i++) {
        std::ranges::set_intersection(
            unambiguousOrientedReadIdsAllBranches,
            allOrientedReadIds[i],
            back_inserter(unambiguousOrientedReadIds[i]));
        if(debug) {
            cout << "Branch " << i << " has " << unambiguousOrientedReadIds[i].size() <<
                " unambiguous oriented read ids." << endl;
        }
    }
#endif

    // Create an AnchorPair between the source and target of this bubble.
    const AnchorId anchorId0 = assemblyGraph[bubble.v0].anchorId;
    const AnchorId anchorId1 = assemblyGraph[bubble.v1].anchorId;
    const AnchorPair bubbleAnchorPair(anchors, anchorId0, anchorId1, false);

    // The usable OrientedReadIds for each branch are the intersection of
    // allOrientedReads for the branch with the OrientedReadIds in the bubbleAnchorPair.
    vector< vector<OrientedReadId> > usableOrientedReadIds(ploidy);
    for(uint64_t i=0; i<ploidy; i++) {
        std::ranges::set_intersection(
            bubbleAnchorPair.orientedReadIds,
            allOrientedReadIds[i],
            back_inserter(usableOrientedReadIds[i]));
        if(debug) {
            cout << "Branch " << i << " has " << usableOrientedReadIds[i].size() <<
                " usable oriented read ids." << endl;
        }
    }



    // Compute groups of similar branches. Each group will generate a new branch.
    // This uses a disjoint sets data structure created from the similarPairs vector,
    // except for some common special cases.
    vector< vector<uint64_t> > branchGroups;
    if(ploidy == 2) {

        SHASTA2_ASSERT(similarPairs.size() == 1);    // We already checked that is is not empty.
        // Create a single branch group that includes both branches.
        branchGroups.push_back({0, 1});

    } else if(similarPairs.size() == (ploidy * (ploidy - 1)) / 2) {

        // Create a single branch group that includes all branches.
        branchGroups.resize(1);
        for(uint64_t i=0; i<ploidy; i++) {
            branchGroups.front().push_back(i);
        }

    } else {

        // Initialize the disjoint sets data structure.
        vector<uint64_t> rank(ploidy);
        vector<uint64_t> parent(ploidy);
        boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
        for(uint64_t i=0; i<ploidy; i++) {
            disjointSets.make_set(i);
        }

        // Connect the similar pairs.
        for(const auto& p: similarPairs) {
            disjointSets.union_set(p.first, p.second);
        }

        // Gather the branches in each group.
        vector< vector<uint64_t> > groups(ploidy);
        for(uint64_t i=0; i<ploidy; i++) {
            groups[disjointSets.find_set(i)].push_back(i);
        }

        // Each of the non-empty ones generates a branch group.
        for(const auto& group: groups) {
            if(not group.empty()) {
                branchGroups.push_back(group);
            }
        }
    }

    if(debug) {
        cout << "Found the following similar groups:" << endl;
        for(const auto& branchGroup: branchGroups) {
            for(const uint64_t i: branchGroup) {
                cout << i << " ";
            }
            cout << endl;
        }
    }



    // Each branch group generates a new branch, if it has enough coverage.
    for(const auto& branchGroup: branchGroups) {

        // If this branch group consists of a single branch, don't do anything.
        if(branchGroup.size() == 1) {
            continue;
        }

        if(debug) {
            cout << "Working on branch group";
            for(const uint64_t i: branchGroup) {
                cout << " " << i;
            }
            cout << endl;
        }

        // Create an AnchorPair using all of the usableOrientedReadIds
        // for the branches in this group.
        AnchorPair newAnchorPair;
        newAnchorPair.anchorIdA = anchorId0;
        newAnchorPair.anchorIdB = anchorId1;
        for(const uint64_t i: branchGroup) {
            std::ranges::copy(usableOrientedReadIds[i], back_inserter(newAnchorPair.orientedReadIds));
        }
        deduplicate(newAnchorPair.orientedReadIds);
        if(debug) {
            cout << "The new branch has coverage " << newAnchorPair.size() << endl;
        }

        if(newAnchorPair.size() < options.bubbleCleanupMinCommonCount) {
            if(debug) {
                cout << "Coverage for this branch group is too low, no new branch generated." << endl;
            }
            continue;
        }

        // Here we are making changes to the AssemblyGraph, so we need to grab the mutex.
        std::lock_guard<std::mutex> lock(mutex);
        edge_descriptor e;
        tie(e, ignore) = add_edge(bubble.v0, bubble.v1, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
        AssemblyGraphEdge& edge = assemblyGraph[e];
        edge.emplace_back(newAnchorPair, newAnchorPair.getAverageOffset(anchors));
        if(debug) {
            cout << "Generated a new branch " << edge.id << " by combining the branches in this group." << endl;
        }

        // Now we can delete the old branches in this group.
        for(const uint64_t i: branchGroup) {
            if(debug) {
                cout << "Removing " << assemblyGraph[bubble.edges[i]].id << endl;
            }
            boost::remove_edge(bubble.edges[i], assemblyGraph);
        }
    }


    return true;
}



// Analyze a Bubble and finds pairs of "similar" branches.
// See analyzeSimilarSequences.h for more information.
bool AssemblyGraph::analyzeBubble(
    const Bubble& bubble,
    const vector<uint64_t> minRepeatCount,
    vector< pair<uint64_t, uint64_t> >& similarPairs
    ) const
{
    const AssemblyGraph& assemblyGraph = *this;
    using shasta2::Base;

    const bool debug = false; // (assemblyGraph[bubble.edges.front()].id == 130581);
    if(debug) {
        cout << "Analyzing bubble";
        for (const edge_descriptor e: bubble.edges) {
            cout << " " << assemblyGraph[e].id;
        }
        cout << endl;
    }

    SHASTA2_ASSERT(bubble.edges.size() > 1);

    // Gather the sequences of all the sides of this bubble
    vector< vector<Base> > sequences;
    for(const edge_descriptor e: bubble.edges) {
        sequences.emplace_back();
        SHASTA2_ASSERT(assemblyGraph[e].wasAssembled);
        assemblyGraph[e].getSequence(sequences.back());
    }

    if(debug) {
        cout << "Bubble sequences:" << endl;
        for(uint64_t i=0; i<bubble.edges.size(); i++) {
            const edge_descriptor e = bubble.edges[i];
            const AssemblyGraphEdge& edge = assemblyGraph[e];
            cout << ">" << edge.id <<  " " << sequences[i].size() << endl;
            copy(sequences[i], ostream_iterator<Base>(cout));
            cout << endl;
        }
    }

    // Loop over pairs of bubble edges.
    ostream html(0);
    similarPairs.clear();
    for(uint64_t i0=0; i0<bubble.edges.size()-1; i0++) {
        const vector<Base>& sequence0 = sequences[i0];
        for(uint64_t i1=i0+1; i1<bubble.edges.size(); i1++) {
            const vector<Base>& sequence1 = sequences[i1];
             if(debug) {
                cout << "Checking " << assemblyGraph[bubble.edges[i0]].id <<
                    " against " << assemblyGraph[bubble.edges[i1]].id << endl;
            }

            if(areSimilarSequences(sequence0, sequence1, minRepeatCount, html)) {
                similarPairs.emplace_back(i0, i1);
            }
        }

    }

    return false;
}






// Find pairs of reverse complemented bubbles.
// The edges of the first bubble in each pair are sorted by id.
// The edges of the second bubble in each pair are sorted
// consistently with the ones in the first pair,
// that is, the reverse complement of p.first.edges[i] is p.second.edges[i].
void AssemblyGraph::findBubblePairs(vector<BubblePair>& bubblePairs) const
{
    performanceLog << timestamp << "AssemblyGraph::findBubblePairs begins." << endl;

    const AssemblyGraph& assemblyGraph = *this;
    bubblePairs.clear();

    // Each iteration of the loop generates a bubble and its reverse
    // complement. To generate each bubble only once we need to
    // keep track of the bubbles we already found.
    std::set< pair<vertex_descriptor, vertex_descriptor> > bubblesFound;

    // Look at bubbles with source v0.
    std::map<vertex_descriptor, vector<edge_descriptor> > m;
    BGL_FORALL_VERTICES(v0, assemblyGraph, AssemblyGraph) {

        // Gather edges, grouped by their target vertex.
        m.clear();
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            m[v1].push_back(e);
        }


        // Loop over groups with the same target vertex.
        // Each group with more than one edge generates
        // a pair of reverse complemented bubbles.
        for(auto& p: m) {
            const vertex_descriptor v1 = p.first;
            vector<edge_descriptor>& edges = p.second;
            if(edges.size() > 1) {

                // If we aleady generated this bubble, don't add it again.
                const vertex_descriptor v0Rc = assemblyGraph[v0].vRc;
                const vertex_descriptor v1Rc = assemblyGraph[v1].vRc;
                if(bubblesFound.contains(make_pair(v0, v1))) {
                    SHASTA2_ASSERT(bubblesFound.contains(make_pair(v1Rc, v0Rc)));
                    continue;
                }
                if(bubblesFound.contains(make_pair(v1Rc, v0Rc))) {
                    SHASTA2_ASSERT(bubblesFound.contains(make_pair(v0, v1)));
                    continue;
                }

                // Ok, if getting here we are going to generate a new pair of bubbles.
                bubblesFound.insert(make_pair(v0, v1));
                bubblesFound.insert(make_pair(v1Rc, v0Rc));

                std::ranges::sort(edges, orderById);
                auto& [bubble, bubbleRc] = bubblePairs.emplace_back();
                bubble.v0 = v0;
                bubble.v1 = v1;
                bubble.edges = edges;
                bubbleRc.v0 = v1Rc;
                bubbleRc.v1 = v0Rc;
                for(const edge_descriptor e: edges) {
                    bubbleRc.edges.push_back(assemblyGraph[e].eRc);
                }
            }
        }
    }

    performanceLog << timestamp << "AssemblyGraph::findBubblePairs ends." << endl;

}



uint64_t AssemblyGraph::bubblePairCleanup()
{
    uint64_t modifiedCount = 0;

    vector< pair<vertex_descriptor, vertex_descriptor> > excludeList;
    while(true) {
        const uint64_t modifiedCountThisIteration = bubblePairCleanupIterationMultithreaded(excludeList);
        if(modifiedCountThisIteration == 0) {
            break;
        }
        modifiedCount += modifiedCountThisIteration;
    }

    return modifiedCount;
}



uint64_t AssemblyGraph::bubblePairCleanupIterationMultithreaded(
    vector< pair<vertex_descriptor, vertex_descriptor> >& excludeList)
{

    performanceLog << timestamp << "Bubble cleanup iteration begins." << endl;
    AssemblyGraph& assemblyGraph = *this;

    // Find all the bubble pairs.
    vector<BubblePair> allBubblePairs;
    findBubblePairs(allBubblePairs);
    cout << "Found " << 2 * allBubblePairs.size() << " bubbles." << endl;
    performanceLog << timestamp << "Found " << 2 * allBubblePairs.size() << " bubbles." << endl;



    // Find candidate bubble pairs.
    // These are the ones that are not in the exclude list and
    // in which no branch has offset greater than
    // options.bubbleCleanupMaxBubbleLength.
    vector<BubblePair>& candidateBubblePairs = bubblePairCleanupIterationData.candidateBubblePairs;
    candidateBubblePairs.clear();
    for(const auto& [bubbleA, bubbleB]: allBubblePairs) {


        if(std::ranges::binary_search(excludeList, make_pair(bubbleA.v0, bubbleA.v1))) {
            SHASTA2_ASSERT(std::ranges::contains(excludeList, make_pair(bubbleB.v0, bubbleB.v1)));
            continue;
        }

        bool hasLongBranch = false;
        for(const edge_descriptor e: bubbleA.edges) {
            if(assemblyGraph[e].offset() > options.bubbleCleanupMaxBubbleLength) {
                hasLongBranch = true;
                break;
            }
        }

        if(not hasLongBranch) {
            candidateBubblePairs.push_back(BubblePair(bubbleA, bubbleB));
        }
    }
    cout << 2 * candidateBubblePairs.size() << " bubbles are candidate for clean up." << endl;
    performanceLog << timestamp << 2 * candidateBubblePairs.size() << " bubbles are candidate for clean up." << endl;



    // Assemble sequence for all the edges of the first bubble of each BubblePair.
    edgesToBeAssembled.clear();
    for(const auto& [bubbleA, bubbleB]: allBubblePairs) {
        for(const edge_descriptor e: bubbleA.edges) {
            if(not assemblyGraph[e].wasAssembled) {
                edgesToBeAssembled.push_back(e);
            }
        }
    }
    assemble();



    // Process the bubble pairs in multithreaded code.
    // All changes to the AssemblyGraph are done under mutex.
    bubblePairCleanupIterationData.modifiedCount = 0;
    const uint64_t batchSize = 10;
    performanceLog << timestamp << "Begin processing bubble pairs." << endl;
    setupLoadBalancing(candidateBubblePairs.size(), batchSize);
    runThreads(&AssemblyGraph::bubblePairCleanupIterationThreadFunction, options.actualThreadCount());
    cout << "Bubble cleanup modified " << bubblePairCleanupIterationData.modifiedCount << " bubbles." << endl;

    // Update the excludeList.
    for(const auto&[bubbleA, bubbleB]: candidateBubblePairs) {
        excludeList.push_back(make_pair(bubbleA.v0, bubbleA.v1));
        excludeList.push_back(make_pair(bubbleB.v0, bubbleB.v1));
    }
    std::ranges::sort(excludeList);

    performanceLog << timestamp << "Bubble cleanup iteration ends." << endl;


    return bubblePairCleanupIterationData.modifiedCount;
}



void AssemblyGraph::bubblePairCleanupIterationThreadFunction([[maybe_unused]] uint64_t threadId)
{
    vector<BubblePair>& candidateBubblePairs = bubblePairCleanupIterationData.candidateBubblePairs;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over candidate bubble pairs in this batch.
        for(uint64_t i=begin; i!=end; ++i) {
            if(bubblePairCleanup(candidateBubblePairs[i])) {
                __sync_fetch_and_add(&bubblePairCleanupIterationData.modifiedCount, 2);
            }
        }
    }
}



bool AssemblyGraph::bubblePairCleanup(const BubblePair& bubblePair)
{
    AssemblyGraph& assemblyGraph = *this;
    const bool debug = false;

    const Bubble& bubbleA = bubblePair.first;
    const Bubble& bubbleB = bubblePair.second;

    if(debug) {
        cout << "Bubble pair cleanup begins for:" << endl;
        cout << "Bubble A " <<
            anchorIdToString(assemblyGraph[bubbleA.v0].anchorId) << " ... " <<
            anchorIdToString(assemblyGraph[bubbleA.v1].anchorId) << endl;
        cout << "Bubble A edges before cleanup:";
        for(const edge_descriptor e: bubbleA.edges) {
            cout << " " << assemblyGraph[e].id;
        }
        cout << endl;
        cout << "Bubble B " <<
            anchorIdToString(assemblyGraph[bubbleB.v0].anchorId) << "..." <<
            anchorIdToString(assemblyGraph[bubbleB.v1].anchorId) << endl;
        cout << "Bubble B edges before cleanup:";
        for(const edge_descriptor e: bubbleB.edges) {
            cout << " " << assemblyGraph[e].id;
        }
        cout << endl;

    }

    if(bubbleCleanup(bubbleA)) {

        if(debug) {
            cout << "Edges of bubble A after cleanup:" << endl;
            BGL_FORALL_OUTEDGES(bubbleA.v0, eA, assemblyGraph, AssemblyGraph) {
                if(target(eA, assemblyGraph) != bubbleA.v1) {
                    continue;
                }
                const auto& edgeA = assemblyGraph[eA];
                cout << "    Edge " << edgeA.id << " " <<
                    anchorIdToString(assemblyGraph[source(eA, assemblyGraph)].anchorId) << " " <<
                    anchorIdToString(assemblyGraph[target(eA, assemblyGraph)].anchorId) << " with " <<
                    edgeA.size() << " steps." << endl;
                for(const AssemblyGraphEdgeStep& step: edgeA) {
                    cout << "        Step " << anchorIdToString(step.anchorPair.anchorIdA) << " " <<
                        anchorIdToString(step.anchorPair.anchorIdB) << endl;
                }
            }
        }

        // bubbleA was changed.
        // To maintain strand symmetry, we need to make
        // corresponding changes in bubbleB.

        // Here we are making changes to the AssemblyGraph, so we need to grab the mutex.
        std::lock_guard<std::mutex> lock(mutex);

        // First, remove the existing edges of bubbleB.
        for(const edge_descriptor e: bubbleB.edges) {
            boost::remove_edge(e, assemblyGraph);
        }

        // Now add edges to bubbleB, copying them one by one
        // from bubbleA, with reverse complement.
        // This also updates the eRc fields of eA and the newly aded edges.
        BGL_FORALL_OUTEDGES(bubbleA.v0, eA, assemblyGraph, AssemblyGraph) {
            if(target(eA, assemblyGraph) != bubbleA.v1) {
                continue;
            }
            const edge_descriptor eB = createReverseComplementEdge(eA);

            // Sanity checks.
            const AssemblyGraphEdge& edgeB = assemblyGraph[eB];
            if(debug) {
                cout << "Newly created edge " << edgeB.id <<
                    " is the reverse complement of  edge " << assemblyGraph[eA].id << endl;
                cout << "Steps of newly created edge:" << endl;
                for(const AssemblyGraphEdgeStep& step: edgeB) {
                    cout << anchorIdToString(step.anchorPair.anchorIdA) << " " <<
                        anchorIdToString(step.anchorPair.anchorIdB) << endl;
                }
            }
            SHASTA2_ASSERT(edgeB.front().anchorPair.anchorIdA == assemblyGraph[bubbleB.v0].anchorId);
            SHASTA2_ASSERT(edgeB.back().anchorPair.anchorIdB == assemblyGraph[bubbleB.v1].anchorId);
        }

        return true;

    } else {

        // bubbleA did not change, so we also leave bubbleB unchanged.
        return false;
    }
}



// This creates a vertex vB, identical to the reverse complement
// of vertex vA. It sets the vRc fields in vA and vB
// and returns vB.
AssemblyGraph::vertex_descriptor AssemblyGraph::createReverseComplementVertex(vertex_descriptor vA)
{
    const bool debug = false;

    AssemblyGraph& assemblyGraph = *this;
    AssemblyGraphVertex& vertexA = assemblyGraph[vA];

    // Get the AnchorIds.
    const AnchorId anchorIdA = vertexA.anchorId;
    const AnchorId anchorIdB = reverseComplementAnchorId(anchorIdA);

    // Create the new vertex.
    const vertex_descriptor vB = add_vertex(AssemblyGraphVertex(anchorIdB, nextVertexId++), assemblyGraph);
    AssemblyGraphVertex& vertexB = assemblyGraph[vB];

    // Set the vRc fields.
    vertexA.vRc = vB;
    vertexB.vRc = vA;

    if(debug) {
        cout << "Created a reverse complemented copy of a vertex: " <<
            vertexA.id << " " << anchorIdToString(vertexA.anchorId) <<
            vertexB.id << " " << anchorIdToString(vertexB.anchorId) << endl;

    }

    return vB;
}



// This creates an edge eB, identical to the reverse complement
// of edge eA. It sets the eRc fields in eA and eB
// and returns eB.
// The vertices of eB must already exist.
AssemblyGraph::edge_descriptor AssemblyGraph::createReverseComplementEdge(edge_descriptor eA)
{
    AssemblyGraph& assemblyGraph = *this;
    const bool debug = false;

    // Access the eA edge and its vertices.
    AssemblyGraphEdge& edgeA = assemblyGraph[eA];
    const vertex_descriptor vA0 = source(eA, assemblyGraph);
    const vertex_descriptor vA1 = target(eA, assemblyGraph);
    const AssemblyGraphVertex& vertexA0 = assemblyGraph[vA0];
    const AssemblyGraphVertex& vertexA1 = assemblyGraph[vA1];
    // cout << "MMM " << vertexA0.id << " " << vertexA1.id << endl;

    // Create the new edge eB.
    const vertex_descriptor vB0 = vertexA1.vRc;
    const vertex_descriptor vB1 = vertexA0.vRc;
    // cout << "NNN " << assemblyGraph[vB0].id << " " << assemblyGraph[vB1].id << endl;
    auto[eB, wasAdded] = add_edge(vB0, vB1, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
    // cout << "OOO" << endl;
    SHASTA2_ASSERT(wasAdded);
    AssemblyGraphEdge& edgeB = assemblyGraph[eB];

    // Update the eRc fields of both edges.
    edgeA.eRc = eB;
    edgeB.eRc = eA;

    // cout << "PPP " << edgeA.id << " " << edgeB.id << endl;
    // cout << "QQQ " << edgeA.size() << " " << edgeB.size() << endl;



    // Now add the steps of eA to eB, with reverse complement.
    for(const AssemblyGraphEdgeStep& stepA: edgeA) {
        AssemblyGraphEdgeStep& stepB = edgeB.emplace_back();
        stepB.offset = stepA.offset;
        stepB.anchorPair = stepA.anchorPair;

        // Swap and reverse complement the AnchorIds.
        std::swap(stepB.anchorPair.anchorIdA, stepB.anchorPair.anchorIdB);
        stepB.anchorPair.anchorIdA = reverseComplementAnchorId(stepB.anchorPair.anchorIdA);
        stepB.anchorPair.anchorIdB = reverseComplementAnchorId(stepB.anchorPair.anchorIdB);

        // Reverse complement the OrientedReadIds.
        for(OrientedReadId& orientedReadId: stepB.anchorPair.orientedReadIds) {
            orientedReadId.flipStrand();
        }
    }
    // Finally, reverse the steps of eB.
    std::ranges::reverse(edgeB);



    if(debug) {
        cout << "addReverseComplementEdge details:" << endl;

        cout << "edgeA: " << edgeA.id << " " <<
            anchorIdToString(assemblyGraph[source(eA, assemblyGraph)].anchorId) << " " <<
            anchorIdToString(assemblyGraph[target(eA, assemblyGraph)].anchorId) << endl;
        cout << "edgeA steps:" << endl;
        for(const auto& step: edgeA) {
            cout << anchorIdToString(step.anchorPair.anchorIdA) << " " <<
                anchorIdToString(step.anchorPair.anchorIdB) << endl;
        }

        cout << "edgeB: " << edgeB.id << " " <<
            anchorIdToString(assemblyGraph[source(eB, assemblyGraph)].anchorId) << " " <<
            anchorIdToString(assemblyGraph[target(eB, assemblyGraph)].anchorId) << endl;
        cout << "edgeB steps:" << endl;
        for(const auto& step: edgeB) {
            cout << anchorIdToString(step.anchorPair.anchorIdA) << " " <<
                anchorIdToString(step.anchorPair.anchorIdB) << endl;
        }
    }

    // Sanity check.
    if(not edgeB.empty()) {
        SHASTA2_ASSERT(edgeB.front().anchorPair.anchorIdA == assemblyGraph[vB0].anchorId);
        SHASTA2_ASSERT(edgeB.back().anchorPair.anchorIdB == assemblyGraph[vB1].anchorId);
    }

    return eB;
}
