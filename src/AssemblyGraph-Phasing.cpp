#include "AssemblyGraph.hpp"
#include "color.hpp"
#include "deduplicate.hpp"
#include "findConvergingVertex.hpp"
#include "findReachableVertices.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
#include "SuperbubbleChain.hpp"
#include "timestamp.hpp"
using namespace shasta2;



void AssemblyGraph::phaseSuperbubbleChains()
{
    performanceLog << timestamp << "AssemblyGraph::phaseSuperbubbleChains begins." << endl;

    PhaseSuperbubbleChainsData& data = phaseSuperbubbleChainsData;
    vector<SuperbubbleChain>& superbubbleChains = data.superbubbleChains;

    // Find superbubbles.
    vector<Superbubble> superbubbles;
    findSuperbubbles(superbubbles);
    writeSuperbubbles(superbubbles, "Superbubbles-WithOverlaps.csv");
    removeContainedSuperbubbles(superbubbles);
    cout << "Found " << superbubbles.size() << " non-overlapping superbubbles." << endl;
    writeSuperbubbles(superbubbles, "Superbubbles.csv");
    writeSuperbubblesForBandage(superbubbles, "Superbubbles-Bandage.csv");

    // Find Superbubble chains.
    findSuperbubbleChains(superbubbles, superbubbleChains);
    cout << "Found " << superbubbleChains.size() << " superbubble chains." << endl;
    writeSuperbubbleChains(superbubbleChains, "SuperbubbleChains.csv");
    writeSuperbubbleChainsForBandage(superbubbleChains, "SuperbubbleChains-Bandage.csv");

    // Phase them.
    setupLoadBalancing(superbubbleChains.size(), 1);
    runThreads(&AssemblyGraph::phaseSuperbubbleChainsThreadFunction, options.threadCount);
    superbubbleChains.clear();
    superbubbleChains.shrink_to_fit();

    compress();
    performanceLog << timestamp << "AssemblyGraph::phaseSuperbubbleChains ends." << endl;
}



void AssemblyGraph::phaseSuperbubbleChainsThreadFunction([[maybe_unused]] uint64_t threadId)
{
    PhaseSuperbubbleChainsData& data = phaseSuperbubbleChainsData;
    vector<SuperbubbleChain>& superbubbleChains = data.superbubbleChains;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all superbubble chains assigned to this batch.
        for(uint64_t superbubbleChainId=begin; superbubbleChainId<end; superbubbleChainId++) {
            SuperbubbleChain& superbubbleChain = superbubbleChains[superbubbleChainId];
            superbubbleChain.phase1(
                *this,
                superbubbleChainId);
        }
    }
}



void AssemblyGraph::findSuperbubbleChains(
    const vector<Superbubble>& superbubbles,
    vector<SuperbubbleChain>& superbubbleChains
    ) const
{
    // Index the superbubbles by their source and target vertex.
    // FOR REPRODUCIBILITY WE SHOULD CHANGE THIS TO INDEX BY VERTEX ID INSTEAD.
    std::map<vertex_descriptor, vector<uint64_t> > mapBySource;
    std::map<vertex_descriptor, vector<uint64_t> > mapByTarget;
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        const Superbubble& superbubble = superbubbles[superbubbleId];
        mapBySource[superbubble.sourceVertex].push_back(superbubbleId);
        mapByTarget[superbubble.targetVertex].push_back(superbubbleId);
    }

    // Sanity check: a vertex can only be a source or target of a single Superbubble.
    for(const auto& p: mapBySource) {
        SHASTA2_ASSERT(p.second.size() == 1);
    }
    for(const auto& p: mapByTarget) {
        SHASTA2_ASSERT(p.second.size() == 1);
    }

    // A vector to keep track of the Superbubbles we already added to a chain.
    vector<bool> wasUsed(superbubbles.size(), false);

    // Work vectors used below to construct chains.
    vector<uint64_t> forward;
    vector<uint64_t> backward;

    // Generate Superbubble chains.
    superbubbleChains.clear();
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        if(wasUsed[superbubbleId]) {
            continue;
        }
        // cout << "Starting a new chain at Superbubble " << superbubbleId << endl;

        // This Superbubble has not yet been added to any chain.
        // We will start a new chain here.

        // Create the forward portion of this chain.
        forward.clear();
        vertex_descriptor v = superbubbles[superbubbleId].targetVertex;
        while(true) {
            const auto it = mapBySource.find(v);
            if(it == mapBySource.end()) {
                break;
            }
            const vector<uint64_t>& nextVector = it->second;
            SHASTA2_ASSERT(nextVector.size() == 1);
            const uint64_t nextSuperbubbleId = nextVector.front();
            forward.push_back(nextSuperbubbleId);
            // cout << "Forward: " << nextSuperbubbleId << endl;

            v = superbubbles[nextSuperbubbleId].targetVertex;
        }

        // Create the backward portion of this chain.
        backward.clear();
        v = superbubbles[superbubbleId].sourceVertex;
        while(true) {
            const auto it = mapByTarget.find(v);
            if(it == mapByTarget.end()) {
                break;
            }
            const vector<uint64_t>& previousVector = it->second;
            SHASTA2_ASSERT(previousVector.size() == 1);
            const uint64_t previousSuperbubbleId = previousVector.front();
            backward.push_back(previousSuperbubbleId);
            // cout << "Backward: " << previousSuperbubbleId << endl;

            v = superbubbles[previousSuperbubbleId].sourceVertex;
        }

        // Now we can create the new SuperbubbleChain.
        superbubbleChains.emplace_back();
        SuperbubbleChain& superbubbleChain = superbubbleChains.back();
        std::reverse(backward.begin(), backward.end());
        for(const uint64_t id: backward) {
            wasUsed[id] = true;
            superbubbleChain.push_back(superbubbles[id]);
        }
        superbubbleChain.push_back(superbubbles[superbubbleId]);
        wasUsed[superbubbleId] = true;
        for(const uint64_t id: forward) {
            wasUsed[id] = true;
            superbubbleChain.push_back(superbubbles[id]);
        }
    }
}



void AssemblyGraph::strandSymmetricPhaseSuperbubbleChains()
{
    performanceLog << timestamp << "AssemblyGraph::strandSymmetricPhaseSuperbubbleChains begins." << endl;
    AssemblyGraph& assemblyGraph = *this;

    PhaseSuperbubbleChainsData& data = phaseSuperbubbleChainsData;
    vector<SuperbubbleChain>& superbubbleChains = data.superbubbleChains;

    // Find Superbubbles.
    vector<Superbubble> superbubbles;
    findSuperbubbles(superbubbles);
    writeSuperbubbles(superbubbles, "Superbubbles-WithOverlaps.csv");
    removeContainedSuperbubbles(superbubbles);
    cout << "Found " << superbubbles.size() << " non-overlapping superbubbles." << endl;
    writeSuperbubbles(superbubbles, "Superbubbles.csv");
    writeSuperbubblesForBandage(superbubbles, "Superbubbles-Bandage.csv");

    // Find SuperbubbleChain.
    findSuperbubbleChains(superbubbles, superbubbleChains);
    cout << "Found " << superbubbleChains.size() << " superbubble chains." << endl;
    writeSuperbubbleChains(superbubbleChains, "SuperbubbleChains.csv");
    writeSuperbubbleChainsForBandage(superbubbleChains, "SuperbubbleChains-Bandage.csv");



    // To find pairs of reverse complemented SuperbubbleChains, index the edges by superbubble chain.
    std::map<edge_descriptor, uint64_t> superbubbleChainMap;
    for(uint64_t superbubbleChainId=0; superbubbleChainId<superbubbleChains.size(); superbubbleChainId++) {
        const SuperbubbleChain& superbubbleChain = superbubbleChains[superbubbleChainId];
        for(const Superbubble& superbubble: superbubbleChain) {
            for(const edge_descriptor e: superbubble.internalEdges) {
                SHASTA2_ASSERT(not superbubbleChainMap.contains(e));
                superbubbleChainMap.insert(make_pair(e, superbubbleChainId));
            }
        }
    }
    cout << "Out of " << num_edges(assemblyGraph) << " assembly graph edges, " <<
        superbubbleChainMap.size() << " are in superbubble chains." << endl;



    // Find the reverse complement of each SuperbubbleChain.
    // This is very strict and can be simplified.
    vector<uint64_t> v;
    vector<uint64_t> superchainTable(superbubbleChains.size(), invalid<uint64_t>);
    for(uint64_t superbubbleChainId=0; superbubbleChainId<superbubbleChains.size(); superbubbleChainId++) {
        const SuperbubbleChain& superbubbleChain = superbubbleChains[superbubbleChainId];
        v.clear();
        for(const Superbubble& superbubble: superbubbleChain) {
            for(const edge_descriptor e: superbubble.internalEdges) {
                const edge_descriptor eRc = assemblyGraph[e].eRc;
                SHASTA2_ASSERT(eRc != assemblyGraphNullEdge);
                if(not superbubbleChainMap.contains(eRc)) {
                    cout << "Edge " << assemblyGraph[e].id <<
                        " is in a superbubble chain but its reverse complement " <<
                        assemblyGraph[eRc].id << " is not." << endl;
                }
                v.push_back(superbubbleChainMap.at(eRc));
            }
        }
        deduplicate(v);
        SHASTA2_ASSERT(v.size() == 1);
        superchainTable[superbubbleChainId] = v.front();
    }

    // Sanity check.
    for(uint64_t superbubbleChainId=0; superbubbleChainId<superbubbleChains.size(); superbubbleChainId++) {
        const uint64_t superbubbleChainIdRc = superchainTable[superbubbleChainId];
        SHASTA2_ASSERT(superbubbleChainIdRc != superbubbleChainId);
        SHASTA2_ASSERT(superchainTable[superbubbleChainIdRc] == superbubbleChainId);
    }

    // Fill in the pairs of SuperbubbleChains to be phased.
    data.superbubbleChainPairs.clear();
    for(uint64_t superbubbleChainId=0; superbubbleChainId<superbubbleChains.size(); superbubbleChainId++) {
        const uint64_t superbubbleChainIdRc = superchainTable[superbubbleChainId];
        if(superbubbleChainId < superbubbleChainIdRc) {
            data.superbubbleChainPairs.push_back(make_pair(superbubbleChainId, superbubbleChainIdRc));
        }
    }


    // Now, for each pair of reverse complement superbubble chains,
    // phase the first one, then make a reverse complemented copy to replace the second one.
    const uint64_t batchSize = 1;
    setupLoadBalancing(data.superbubbleChainPairs.size(), batchSize);
    runThreads(&AssemblyGraph::strandSymmetricPhaseSuperbubbleChainsThreadFunction, options.actualThreadCount());
    superbubbleChains.clear();
    superbubbleChains.shrink_to_fit();
    data.superbubbleChainPairs.clear();
    data.superbubbleChainPairs.shrink_to_fit();


    strandSymmetricCompress();
    performanceLog << timestamp << "AssemblyGraph::strandSymmetricPhaseSuperbubbleChains ends." << endl;
}



void AssemblyGraph::strandSymmetricPhaseSuperbubbleChainsThreadFunction([[maybe_unused]] uint64_t threadId)
{
    const PhaseSuperbubbleChainsData& data = phaseSuperbubbleChainsData;
    const vector<SuperbubbleChain>& superbubbleChains = data.superbubbleChains;
    const vector< pair<uint64_t, uint64_t> >& superbubbleChainPairs = data.superbubbleChainPairs;

    // Loop over SuperbubbleChain pairs batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over SuperbubbleChain pairs in this batch.
        for(uint64_t i=begin; i!=end; i++) {
            const auto&[superbubbleChainId, superbubbleChainIdRc] = superbubbleChainPairs[i];
            const SuperbubbleChain& superbubbleChain = superbubbleChains[superbubbleChainId];
            const SuperbubbleChain& superbubbleChainRc = superbubbleChains[superbubbleChainIdRc];
            strandSymmetricPhase(
                superbubbleChain, superbubbleChainId,
                superbubbleChainRc, superbubbleChainIdRc);
        }
    }
}



// This phases the first SuperbubbleChain of a reverse complemented pair,
// then copies it with reverse complement to replace the second SuperbubbleChain.
void AssemblyGraph::strandSymmetricPhase(
    const SuperbubbleChain& superbubbleChain,
    uint64_t superbubbleChainId,
    const SuperbubbleChain& superbubbleChainRc,
    uint64_t superbubbleChainIdRc)
{
    AssemblyGraph& assemblyGraph = *this;

    const bool debug = false;
    if(debug) {
        cout << "AssemblyGraph::strandSymmetricPhase begins for superbubble chains " <<
            superbubbleChainId << " " << superbubbleChainIdRc << endl;
    }

    // Phase the first bubble.
    superbubbleChain.phase1(assemblyGraph, superbubbleChainId);

    // The code below makes changes to the AssemblyGraph, so we have to acquire the mutex.
    std::lock_guard<std::mutex> lock(assemblyGraph.mutex);

    // Remove all internal vertices and edges of the second bubble.
    // Loop over all of its Superbubbles.
    if(debug) {
        cout << "Removing internal vertices and edges of superbubble chain " <<
            superbubbleChainIdRc << endl;
    }
    for(uint64_t i=0; i<superbubbleChainRc.size(); i++) {
        const Superbubble& superbubbleRc = superbubbleChainRc[i];

        // Remove the internal edges of this Superbubble.
        for(const edge_descriptor e: superbubbleRc.internalEdges) {
            boost::remove_edge(e, assemblyGraph);
        }

        // Remove the internal vertices of this Superbubble.
        for(const vertex_descriptor v: superbubbleRc.internalVertices) {
            SHASTA2_ASSERT(in_degree(v, assemblyGraph) == 0);
            SHASTA2_ASSERT(out_degree(v, assemblyGraph) == 0);
            boost::remove_vertex(v, assemblyGraph);
        }
    }

    // We still have to remove the sourceVertex of each Superbubble,
    // except the first one.
    for(uint64_t i=1; i<superbubbleChainRc.size(); i++) {
        const Superbubble& superbubble = superbubbleChainRc[i];
        const vertex_descriptor v = superbubble.sourceVertex;
        SHASTA2_ASSERT(in_degree(v, assemblyGraph) == 0);
        SHASTA2_ASSERT(out_degree(v, assemblyGraph) == 0);
        boost::remove_vertex(v, assemblyGraph);
    }



    // We need to replace the second bubble with a reverse complemented copy
    // of the first bubble, after phasing.
    // For this we need all the vertices and edges of the
    // phased version of the first bubble.
    if(debug) {
        cout << "Finding internal vertices and edges of superbubble chain " <<
            superbubbleChainId << endl;
    }
    std::set<vertex_descriptor> reachableVertices;
    std::set<edge_descriptor> reachableEdges;
    const vertex_descriptor v0 = superbubbleChain.front().sourceVertex;
    const vertex_descriptor v1 = superbubbleChain.back().targetVertex;
    findReachableWithStop(assemblyGraph, v0, v1, reachableVertices, reachableEdges);

    // Make reverse complemented copies of the vertices.
    if(debug) {
        cout << "Making reverse complemented copies of the vertices." << endl;
    }
    for(const vertex_descriptor v: reachableVertices) {
        createReverseComplementVertex(v);
    }

    if(debug) {
        cout << "Making reverse complemented copies of the edges." << endl;
    }
    // Make reverse complemented copies of the edges.
    for(const edge_descriptor e: reachableEdges) {
        if(debug) {
            cout << "Making a reverse complemented copy of " << assemblyGraph[e].id << endl;
        }
        createReverseComplementEdge(e);
    }

    if(debug) {
        cout << "AssemblyGraph::strandSymmetricPhase ends for superbubble chains " <<
            superbubbleChainId << " " << superbubbleChainIdRc << endl;
    }
}



// This finds all Superbubbles seen using the specified maxDistance.
// Some pairs of Superbubble can intersect (that is, they can have common edges).
void AssemblyGraph::findSuperbubbles(
    vector<Superbubble>& superbubbles) const
{
    const AssemblyGraph& assemblyGraph = *this;

    vector< pair<vertex_descriptor, vertex_descriptor> > forwardPairs;
    BGL_FORALL_VERTICES(vA, assemblyGraph, AssemblyGraph) {
        if(hasSelfEdge(vA)) {
            continue;
        }
        const vertex_descriptor vB = findConvergingVertexGeneral(assemblyGraph, vA, options.findSuperbubblesMaxDistance);
        if(vB != null_vertex()) {
            if(not hasSelfEdge(vB)) {
                forwardPairs.emplace_back(vA, vB);
                // cout << assemblyGraph[vA].id << "..." << assemblyGraph[vB].id << endl;
            }
        }
    }
    sort(forwardPairs.begin(), forwardPairs.end());
    // cout << "Found " << forwardPairs.size() << " forward pairs." << endl;

    const boost::reverse_graph<AssemblyGraph> reverseAssemblyGraph(assemblyGraph);
    vector< pair<vertex_descriptor, vertex_descriptor> > backwardPairs;
    BGL_FORALL_VERTICES(vA, assemblyGraph, AssemblyGraph) {
        if(hasSelfEdge(vA)) {
            continue;
        }
        const vertex_descriptor vB = findConvergingVertexGeneral(reverseAssemblyGraph, vA, options.findSuperbubblesMaxDistance);
        if(vB != null_vertex()) {
            if(not hasSelfEdge(vB)) {
                backwardPairs.emplace_back(vB, vA);
                // cout << assemblyGraph[vA].id << "..." << assemblyGraph[vB].id << endl;
            }
        }
    }
    sort(backwardPairs.begin(), backwardPairs.end());
    // cout << "Found " << backwardPairs.size() << " backward pairs." << endl;

    // Find pairs that appeared in both directions.
    vector< pair<vertex_descriptor, vertex_descriptor> > bidirectionalPairs;
    std::set_intersection(
        forwardPairs.begin(), forwardPairs.end(),
        backwardPairs.begin(), backwardPairs.end(),
        back_inserter(bidirectionalPairs));

    // Each of these pairs generates a Superbubble.
    superbubbles.clear();
    for(const auto& p: bidirectionalPairs) {
        superbubbles.emplace_back(assemblyGraph, p.first, p.second);

    }

}



// Remove Superbubbles that are entirely contained in a larger superbubble.
void AssemblyGraph::removeContainedSuperbubbles(vector<Superbubble>& superbubbles) const
{
    const AssemblyGraph& assemblyGraph = *this;

    // Find pairs of intersecting superbubbles.
    // Two superbubbles intersect if they have one or more internal edges in common.
    std::map<edge_descriptor, vector<uint64_t> > m;
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        const Superbubble& superbubble = superbubbles[superbubbleId];
        for(const edge_descriptor e: superbubble.internalEdges) {
            m[e].push_back(superbubbleId);
        }
    }
    std::set< pair<uint64_t, uint64_t> > intersectingPairs;
    for(const auto& p: m) {
        const vector<uint64_t>& edgeSuperbubbles = p.second;
        for(uint64_t i0=0; i0<edgeSuperbubbles.size()-1; i0++) {
            const uint64_t superbubbleId0 = edgeSuperbubbles[i0];
            for(uint64_t i1=i0+1; i1<edgeSuperbubbles.size(); i1++) {
                const uint64_t superbubbleId1 = edgeSuperbubbles[i1];
                intersectingPairs.insert({
                    min(superbubbleId0, superbubbleId1),
                    max(superbubbleId0, superbubbleId1)});
            }
        }
    }
    // cout << "Found " << intersectingPairs.size() << " intersecting superbubble pairs:" << endl;

    vector<edge_descriptor> commonEdges;
    vector<uint64_t> superbubblesToBeRemoved;
    for(const auto& p: intersectingPairs) {
        const uint64_t superbubbleId0 = p.first;
        const uint64_t superbubbleId1 = p.second;
        const Superbubble& superbubble0 = superbubbles[superbubbleId0];
        const Superbubble& superbubble1 = superbubbles[superbubbleId1];
        const auto& internalEdges0 = superbubble0.internalEdges;
        const auto& internalEdges1 = superbubble1.internalEdges;

        // Find the common edges.
        commonEdges.clear();
        std::set_intersection(
            internalEdges0.begin(), internalEdges0.end(),
            internalEdges1.begin(), internalEdges1.end(),
            back_inserter(commonEdges),
            orderById);

        if(commonEdges.size() == internalEdges0.size()) {
            // cout << "Superbubble " << superbubbleId0 << " is contained in superbubble " << superbubbleId1 << endl;
            superbubblesToBeRemoved.push_back(superbubbleId0);
        } else if(commonEdges.size() == internalEdges1.size()) {
            // cout << "Superbubble " << superbubbleId1 << " is contained in superbubble " << superbubbleId0 << endl;
            superbubblesToBeRemoved.push_back(superbubbleId1);
        } else {
            cout << "Superbubbles " << superbubbleId0 << " and " << superbubbleId1 << " intersect." << endl;
            cout << "Superbubbles " << superbubbleId0 << " has " << internalEdges0.size() << " internal edges:";
            for(const edge_descriptor e: internalEdges0) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << endl;
            cout << "Superbubbles " << superbubbleId1 << " has " << internalEdges1.size() << " internal edges:";
            for(const edge_descriptor e: internalEdges1) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << endl;
            cout << "Found " << commonEdges.size() << " common edges." << endl;
            // This should never happen, but just in case we remove both of them.
            SHASTA2_ASSERT(0);
            // superbubblesToBeRemoved.push_back(superbubbleId0);
            // superbubblesToBeRemoved.push_back(superbubbleId1);
        }
    }
    deduplicate(superbubblesToBeRemoved);
    // cout << "Removing " << superbubblesToBeRemoved.size() << " superbubbles that "
    //     "are contained in another superbubble." << endl;

    vector<Superbubble> newSuperbubbles;
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        if(not binary_search(superbubblesToBeRemoved.begin(), superbubblesToBeRemoved.end(), superbubbleId)) {
            newSuperbubbles.push_back(superbubbles[superbubbleId]);
        }
    }
    newSuperbubbles.swap(superbubbles);
}



// This creates a csv file with one line of information for each superbubble.
void AssemblyGraph::writeSuperbubbles(
    const vector<Superbubble>& superbubbles,
    const string& fileName) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream csv(fileName);
    csv << "Superbubble id,Type,Source ploidy,Target ploidy,Internal edges\n";

    for(uint64_t id=0; id<superbubbles.size(); id++) {
        const Superbubble& superbubble = superbubbles[id];

#if 0
        for(const auto& e: superbubble.sourceEdges) {
            cout << "Source edge " << assemblyGraph[e].id << endl;
        }
        for(const auto& e: superbubble.targetEdges) {
            cout << "Target edge " << assemblyGraph[e].id << endl;
        }
        for(const auto& e: superbubble.internalEdges) {
            cout << "Internal edge " << assemblyGraph[e].id << endl;
        }
#endif

        csv << id << ",";

        if(superbubble.isBubble()) {
            const uint64_t ploidy = superbubble.ploidy();
            if(ploidy == 1) {
                csv << "Edge,";
            } else {
                csv << "Bubble,";
            }
        } else {
            csv << "Superbubble,";
        }

        csv << superbubble.sourcePloidy() << ",";
        csv << superbubble.targetPloidy() << ",";

        for(const edge_descriptor e: superbubble.internalEdges) {
            csv << assemblyGraph[e].id << ",";
        }
        csv << "\n";
    }
}



// This creates a csv file that can be loaded in Bandage to see the Superbubbles.
void AssemblyGraph::writeSuperbubblesForBandage(
    const vector<Superbubble>& superbubbles,
    const string& fileName) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream csv(fileName);
    csv << "Segment,Color\n";

    for(uint64_t id=0; id<superbubbles.size(); id++) {
        const Superbubble& superbubble = superbubbles[id];
        const string color = randomHslColor(id, 0.75, 0.5);

        for(const edge_descriptor e: superbubble.internalEdges) {
            csv << assemblyGraph[e].id << "," << color << "\n";
        }
    }

}



void AssemblyGraph::writeSuperbubbleChains(
    const vector<SuperbubbleChain>& superbubbleChains,
    const string& fileName) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream csv(fileName);
    csv << "ChainId,Position,Type,Source ploidy,Target ploidy,"
        "Internal vertices count,Internal edges count,Internal edges\n";

    for(uint64_t chainId=0; chainId<superbubbleChains.size(); chainId++) {
        const  SuperbubbleChain& superbubbleChain = superbubbleChains[chainId];
        for(uint64_t position=0; position<superbubbleChain.size(); position++) {
            const Superbubble& superbubble = superbubbleChain[position];
            csv << chainId << ",";
            csv << position << ",";

            if(superbubble.isBubble()) {
                const uint64_t ploidy = superbubble.ploidy();
                if(ploidy == 1) {
                    csv << "Edge,";
                } else {
                    csv << "Bubble,";
                }
            } else {
                csv << "Superbubble,";
            }

            csv << superbubble.sourcePloidy() << ",";
            csv << superbubble.targetPloidy() << ",";
            csv << superbubble.internalVertices.size() << ",";
            csv << superbubble.internalEdges.size() << ",";
            for(const edge_descriptor e: superbubble.internalEdges) {
                csv << assemblyGraph[e].id << ",";
            }
            csv << "\n";
        }
    }

}



void AssemblyGraph::writeSuperbubbleChainsForBandage(
    const vector<SuperbubbleChain>& superbubbleChains,
    const string& fileName) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream csv(fileName);
    csv << "Segment,Color,Chain,Position\n";

    for(uint64_t chainId=0; chainId<superbubbleChains.size(); chainId++) {
        const SuperbubbleChain& superbubbleChain = superbubbleChains[chainId];
        const string color = randomHslColor(chainId, 0.75, 0.5);

        for(uint64_t position=0; position<superbubbleChain.size(); position++) {
            const Superbubble& superbubble = superbubbleChain[position];
            for(const edge_descriptor e: superbubble.internalEdges) {
                csv <<
                    assemblyGraph[e].id << "," <<
                    color << "," <<
                    chainId << "," <<
                    position << "\n";
            }
        }
    }
}
