// Must include this first due to some Boost include file problem.
#include <boost/pending/disjoint_sets.hpp>

// Shasta.
#include "AssemblyGraph.hpp"
#include "Anchor.hpp"
#include "AnchorGraph.hpp"
#include "areSimilarSequences.hpp"
#include "color.hpp"
#include "deduplicate.hpp"
#include "DisjointSets.hpp"
#include "findLinearChains.hpp"
#include "findConvergingVertex.hpp"
#include "Journeys.hpp"
#include "LocalAssembly3.hpp"
#include "LocalAssembly4.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
#include "RestrictedAnchorGraph.hpp"
#include "Superbubble.hpp"
#include "SuperbubbleChain.hpp"
#include "TangleMatrix1.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/filtered_graph.hpp>

// Standard library.
#include "chrono.hpp"



// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<AssemblyGraph>;



// Initial construction from the AnchorGraph.
AssemblyGraph::AssemblyGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    const AnchorGraph& anchorGraph,
    const Options& options) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<AssemblyGraph>(*this),
    anchors(anchors),
    journeys(journeys),
    options(options),
    orderById(*this)
{
    AssemblyGraph& assemblyGraph = *this;

    // Create a filtered AnchorGraph containing only the edges marked as "useForAssembly".
    class EdgePredicate {
    public:
        bool operator()(const AnchorGraph::edge_descriptor& e) const
        {
            return (*anchorGraph)[e].useForAssembly;
        }
        EdgePredicate(const AnchorGraph& anchorGraph) : anchorGraph(&anchorGraph) {}
        EdgePredicate() : anchorGraph(0) {}
        const AnchorGraph* anchorGraph;
    };
    using FilteredAnchorGraph = boost::filtered_graph<AnchorGraph, EdgePredicate>;
    FilteredAnchorGraph filteredAnchorGraph(anchorGraph, EdgePredicate(anchorGraph));



    // Find linear chains of edges in the FilteredAnchorGraph.
    vector< std::list<AnchorGraph::edge_descriptor> > chains;
    findLinearChains(filteredAnchorGraph, 1, chains);

    // Generate vertices.
    // At this stage there is a vertex for each AnchorGraph vertex
    // that is at the beginning or end of a linear chain,
    // so there is only one vertex for a given AnchorId.
    // However, after detangling there can be more than one vertex
    // with a given AnchorId. So the vertexMap is only used in this constructor.
    std::map<AnchorId, vertex_descriptor> vertexMap;
    for(const auto& chain: chains) {
        const AnchorId anchorId0 = anchorGraph[chain.front()].anchorPair.anchorIdA;
        const AnchorId anchorId1 = anchorGraph[chain.back()].anchorPair.anchorIdB;

        if(not vertexMap.contains(anchorId0)) {
            const vertex_descriptor v0 = add_vertex(AssemblyGraphVertex(anchorId0, nextVertexId++), assemblyGraph);
            vertexMap.insert(make_pair(anchorId0, v0));
        }

        if(not vertexMap.contains(anchorId1)) {
            const vertex_descriptor v1 = add_vertex(AssemblyGraphVertex(anchorId1, nextVertexId++), assemblyGraph);
            vertexMap.insert(make_pair(anchorId1, v1));
        }
    }
    SHASTA2_ASSERT(vertexMap.size() == num_vertices(assemblyGraph));



    // Generate the edges. There is an edge for each linear chain.
    for(const auto& chain: chains) {
        const AnchorId anchorId0 = anchorGraph[chain.front()].anchorPair.anchorIdA;
        const AnchorId anchorId1 = anchorGraph[chain.back()].anchorPair.anchorIdB;

        const vertex_descriptor v0 = vertexMap.at(anchorId0);
        const vertex_descriptor v1 = vertexMap.at(anchorId1);

        edge_descriptor e;
        bool edgeWasAdded;
        tie(e, edgeWasAdded) = add_edge(v0, v1, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
        AssemblyGraphEdge& edge = assemblyGraph[e];

        // Each AnchorGraph edge in the chain contributes a step to this AssemblyGraph edge.
        for(const AnchorGraph::edge_descriptor eA: chain) {
            const AnchorGraphEdge& edgeA = anchorGraph[eA];
            edge.emplace_back(edgeA.anchorPair, edgeA.offset);
        }
    }


    check();
}



// Deserialize constructor.
AssemblyGraph::AssemblyGraph(
    const Anchors& anchors,
    const Journeys& journeys,
    const Options& options,
    const string& stage) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<AssemblyGraph>(*this),
    anchors(anchors),
    journeys(journeys),
    options(options),
    orderById(*this)
{
    load(stage);
}



// Top level assembly function.
void AssemblyGraph::simplifyAndAssemble()
{
    writePerformanceStatistics("AssemblyGraph::simplifyAndAssemble begins");

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minComponentN50 = 100000;

    // Initial output.
    writeIntermediateStageIfRequested("A");



    // The iteration is needed because phasing can simplify some bubbles and create
    // more opportunities for additional bubble cleanup and phasing.
    for(uint64_t iteration=0; iteration<options.simplifyMaxIterationCount; iteration++) {

        uint64_t changeCount = 0;

        // Remove or simplify bubbles likely caused by errors.
        changeCount += bubbleCleanup();
        changeCount += compress();
        writeIntermediateStageIfRequested("B" + to_string(iteration));

        // Phase SuperbubbleChains, considering all hypotheses.
        changeCount += phaseSuperbubbleChains();
        writeIntermediateStageIfRequested("C" + to_string(iteration));

        // Read following.
        changeCount += findAndConnectAssemblyPaths();
        writeIntermediateStageIfRequested("D" + to_string(iteration));

        // Remove isolated vertices and connected components with small N50.
        removeIsolatedVertices();
        removeLowN50Components(minComponentN50);
        writeIntermediateStageIfRequested("E" + to_string(iteration));

        if(changeCount == 0) {
            break;
        }
    }




    // Sequence assembly.
    assembleAll();
    write("Final");
    writeFasta("Final");

    writePerformanceStatistics("AssemblyGraph::simplifyAndAssemble ends");
}



void AssemblyGraph::check() const
{
    const AssemblyGraph& assemblyGraph = *this;

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        SHASTA2_ASSERT(not edge.empty());



        // Check that the first/last AnchorIds of this edge are consistent
        // with the ones in the source/target vertices.
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);

        const AnchorId anchorId0 = assemblyGraph[v0].anchorId;
        const AnchorId anchorId1 = assemblyGraph[v1].anchorId;

        SHASTA2_ASSERT(edge.front().anchorPair.anchorIdA == anchorId0);
        SHASTA2_ASSERT(edge.back().anchorPair.anchorIdB == anchorId1);



        // Check that AnchorPairs in this edge are adjacent to each other.
        for(uint64_t i1=1; i1<edge.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            SHASTA2_ASSERT(edge[i0].anchorPair.anchorIdB == edge[i1].anchorPair.anchorIdA);
        }    }

}



uint64_t AssemblyGraphEdge::offset() const
{
    uint64_t sum = 0;
    for(const auto& step: *this) {
        sum += step.offset;
    }
    return sum;
}



void AssemblyGraph::write(const string& stage)
{
    cout << "Stage " << stage << ": " <<
        num_vertices(*this) << " vertices, " <<
        num_edges(*this) << " edges. Next edge id is " << nextEdgeId << "." << endl;

    save(stage);
    writeGfa("Assembly-" + stage + ".gfa");
    writeGraphviz("Assembly-" + stage + ".dot");
    writeCsv("Assembly-" + stage + ".csv");
}



void AssemblyGraph::writeIntermediateStageIfRequested(const string& name)
{
    if(options.writeIntermediateAssemblyStages) {
        write(name);
    }
}



void AssemblyGraph::writeFasta(const string& stage) const
{
    const AssemblyGraph& assemblyGraph = *this;;

    ofstream fasta("Assembly-" + stage + ".fasta");

    vector<shasta2::Base> sequence;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        edge.getSequence(sequence);

        fasta << ">" << edge.id << "\n";
        copy(sequence.begin(), sequence.end(), ostream_iterator<shasta2::Base>(fasta));
        fasta << "\n";
    }
}



void AssemblyGraph::writeGfa(const string& fileName) const
{
    ofstream gfa(fileName);
    writeGfa(gfa);
}



void AssemblyGraph::writeGfa(ostream& gfa) const
{
    const AssemblyGraph& assemblyGraph = *this;

    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";

    // Each edge generates a gfa segment.
    vector<shasta2::Base> sequence;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        const double coverage = edge.lengthWeightedAverageCoverage();

        // Record type.
        gfa << "S\t";

        // Name.
        gfa << edge.id << "\t";

        // Sequence.
        if(edge.wasAssembled) {
            edge.getSequence(sequence);
            copy(sequence.begin(), sequence.end(), ostream_iterator<shasta2::Base>(gfa));
            const uint64_t length = sequence.size();
            gfa << "\tLN:i:" << length;
            gfa << "\tRC:i:" << uint64_t(std::round(coverage * double(length)));
            gfa << "\n";

        } else {
            if(edge.empty()) {
                gfa << "*\tLN:i:0\n";

            } else {
                const uint64_t offset = edge.offset();
                gfa << "*\tLN:i:" << offset;
                gfa << "\tRC:i:" << uint64_t(std::round(coverage * double(offset)));
                gfa << "\n";
            }
        }
    }



    // For each vertex, generate a link between each pair of
    // incoming/outgoing edges.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        BGL_FORALL_INEDGES(v, e0, assemblyGraph, AssemblyGraph) {
            const uint64_t id0 = assemblyGraph[e0].id;
            BGL_FORALL_OUTEDGES(v, e1, assemblyGraph, AssemblyGraph) {
                const uint64_t id1 = assemblyGraph[e1].id;

                gfa <<
                    "L\t" <<
                    id0 << "\t+\t" <<
                    id1 << "\t+\t*\n";
            }
        }
    }


}



void AssemblyGraph::writeGraphviz(const string& fileName) const
{
    ofstream dot(fileName);
    writeGraphviz(dot);
}



void AssemblyGraph::writeGraphviz(ostream& dot) const
{
    const AssemblyGraph& assemblyGraph = *this;

    dot << "digraph AssemblyGraph {\n";



    // Write the vertices.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
    	const AssemblyGraphVertex& vertex = assemblyGraph[v];
    	dot <<
    		vertex.id <<
    		" [label=\"" << anchorIdToString(vertex.anchorId) << "\\n" << vertex.id << "\"]"
    	    ";\n";
    }



    // Write the edges.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
    	const vertex_descriptor v0 = source(e, assemblyGraph);
    	const vertex_descriptor v1 = target(e, assemblyGraph);
    	const AssemblyGraphVertex& vertex0 = assemblyGraph[v0];
    	const AssemblyGraphVertex& vertex1 = assemblyGraph[v1];
    	dot <<
    	    vertex0.id << "->" <<
    	    vertex1.id <<
    	    " [label=\"" << edge.id << "\\n" <<
    	    (edge.wasAssembled ? edge.sequenceLength() : edge.offset()) <<
    	    "\\n" << edge.size() <<
    	    "\"]"
    	    ";\n";
    }

    dot << "}\n";
}



// Assemble sequence for all edges.
void AssemblyGraph::assembleAll()
{
    writePerformanceStatistics("AssemblyGraph::assembleAll begins");

    const AssemblyGraph& assemblyGraph = *this;

    edgesToBeAssembled.clear();
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        edgesToBeAssembled.push_back(e);
    }
    assemble();
    edgesToBeAssembled.clear();

    writePerformanceStatistics("AssemblyGraph::assembleAll ends");
}



// Assemble sequence for the specified edge.
void AssemblyGraph::assemble(edge_descriptor e)
{
    edgesToBeAssembled.clear();
    edgesToBeAssembled.push_back(e);
    assemble();
}



// Assemble sequence for step i of the specified edge.
// This is the lowest level sequence assembly function and is not multithreaded.
// It runs a LocalAssembly2 on the AnchorPair for that step.
void AssemblyGraph::assembleStep(edge_descriptor e, uint64_t i)
{
    AssemblyGraph& assemblyGraph = *this;
    AssemblyGraphEdge& edge = assemblyGraph[e];
    AssemblyGraphEdgeStep& step = edge[i];

    // cout << "Begin local assembly for edge " << edge.id << " " << " step " << i << endl;

    if(step.anchorPair.anchorIdA == step.anchorPair.anchorIdB) {
        step.sequence.clear();
        return;
    }



    // Let LocalAssembly3 use OrientedReadIds from the previous and next step.
    vector<OrientedReadId> additionalOrientedReadIds;
    if(i > 0) {
        std::ranges::copy(edge[i - 1].anchorPair.orientedReadIds, back_inserter(additionalOrientedReadIds));
    }
    if(i < edge.size() - 1) {
        std::ranges::copy(edge[i + 1].anchorPair.orientedReadIds, back_inserter(additionalOrientedReadIds));
    }
    deduplicate(additionalOrientedReadIds);

    ostream html(0);
    LocalAssembly4 localAssembly(
        anchors,
        options.abpoaMaxLength,
        html,
        false,
        edge[i].anchorPair,
        additionalOrientedReadIds);
    step.sequence = localAssembly.sequence;
}



// Assemble sequence for all edges in the edgesToBeAssembled vector.
// This fills in the stepsToBeAssembled with all steps of those edges,
// then assembles each of the steps in parallel.
void AssemblyGraph::assemble()
{
    performanceLog << timestamp << "Sequence assembly begins for " << edgesToBeAssembled.size() <<
        " assembly graph edges." << endl;
    AssemblyGraph& assemblyGraph = *this;

    stepsToBeAssembled.clear();
    for(const edge_descriptor e: edgesToBeAssembled) {
        AssemblyGraphEdge& edge = assemblyGraph[e];
        for(uint64_t i=0; i<edge.size(); i++) {
            stepsToBeAssembled.push_back(make_pair(e, i));
        }
    }

    const uint64_t batchCount = 1;
    setupLoadBalancing(stepsToBeAssembled.size(), batchCount);
    runThreads(&AssemblyGraph::assembleThreadFunction, options.threadCount);

    // Mark them as assembled.
    for(const edge_descriptor e: edgesToBeAssembled) {
        assemblyGraph[e].wasAssembled = true;
    }

    edgesToBeAssembled.clear();
    stepsToBeAssembled.clear();

    performanceLog << timestamp << "Sequence assembly ends." << endl;
}



void AssemblyGraph::assembleThreadFunction(uint64_t /* threadId */)
{
    AssemblyGraph& assemblyGraph = *this;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all assembly steps assigned to this batch.
        for(uint64_t j=begin; j!=end; j++) {

            if((j % 1000) == 0) {
                std::lock_guard<std::mutex> lock(mutex);
                performanceLog << timestamp << "Starting sequence assembly step " << j << " of " <<
                    stepsToBeAssembled.size() << endl;
            }

            const auto& p = stepsToBeAssembled[j];
            const edge_descriptor e = p.first;
            const uint64_t i = p.second;
            AssemblyGraphEdge& edge = assemblyGraph[e];
            SHASTA2_ASSERT(i < edge.size());

            try {
                const auto t0 = steady_clock::now();
                assembleStep(e, i);
                const auto t1 = steady_clock::now();
                const double t01 = seconds(t1-t0);
                if(t01 > 10.) {
                    std::lock_guard<std::mutex> lock(mutex);
                    performanceLog << "Slow assembly step: segment " <<
                        assemblyGraph[e].id << " step " << i <<
                        " " << t01 << " s." << endl;

                }
            } catch(std::exception&) {
                std::lock_guard<std::mutex> lock(mutex);
                cout << "Error occurred during local assembly step: segment " <<
                    assemblyGraph[e].id << " step " << i << endl;
                throw;
            }
         }
    }
}



// Clear sequence from all steps of all edges.
void AssemblyGraph::clearSequence()
{
    AssemblyGraph& assemblyGraph = *this;

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        AssemblyGraphEdge& edge = assemblyGraph[e];
        edge.wasAssembled= false;
        for(AssemblyGraphEdgeStep& step: edge) {
            step.sequence.clear();
            step.sequence.shrink_to_fit();
        }
    }
}



void AssemblyGraphEdge::getSequence(vector<Base>& sequence) const
{
    sequence.clear();
    for(const auto& step: *this) {
        copy(step.sequence.begin(), step.sequence.end(), back_inserter(sequence));
    }
}



uint64_t AssemblyGraphEdge::sequenceLength() const
{
    SHASTA2_ASSERT(wasAssembled);

    uint64_t length = 0;
    for(const auto& step: *this) {
        length += step.sequence.size();
    }
    return length;
}



uint64_t AssemblyGraphEdge::length() const
{
    if(wasAssembled) {
        return sequenceLength();
    } else {
        return offset();
    }
}



void AssemblyGraph::findBubbles(vector<Bubble>& bubbles) const
{
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
}



uint64_t AssemblyGraph::bubbleCleanup()
{
    uint64_t modifiedCount = 0;

    vector< pair<vertex_descriptor, vertex_descriptor> > excludeList;
    while(true) {
        const uint64_t modifiedCountThisIteration = bubbleCleanupIteration(excludeList);
        if(modifiedCountThisIteration == 0) {
            break;
        }
        modifiedCount += modifiedCountThisIteration;
    }

    return modifiedCount;
}



uint64_t AssemblyGraph::bubbleCleanupIteration(
    vector< pair<vertex_descriptor, vertex_descriptor> >& excludeList)
{
    performanceLog << timestamp << "Bubble cleanup iteration begins." << endl;
    AssemblyGraph& assemblyGraph = *this;

    // Find all bubbles.
    vector<Bubble> allBubbles;
    findBubbles(allBubbles);
    // cout << "Found " << allBubbles.size() << " bubbles." << endl;

    // Find candidate bubbles.
    // These are the ones that are not in the exclude list and
    // in which no branch has offset greater than
    // options.bubbleCleanupMaxBubbleLength.
    vector<Bubble> candidateBubbles;
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


    uint64_t modifiedCount = 0;
    for(const Bubble& bubble: candidateBubbles) {
        if(bubbleCleanup(bubble)) {
            ++modifiedCount;
        }
    }
    // cout << "Bubble cleanup modified " << modifiedCount << " bubbles." << endl;

    // Update the excludeList.
    for(const Bubble& bubble: candidateBubbles) {
        excludeList.push_back(make_pair(bubble.v0, bubble.v1));
    }
    std::ranges::sort(excludeList);

    performanceLog << timestamp << "Bubble cleanup iteration ends." << endl;
    return modifiedCount;
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

    SHASTA2_ASSERT(out_degree(bubble.v0, assemblyGraph) > 0);
    SHASTA2_ASSERT(in_degree(bubble.v1, assemblyGraph) > 0);

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



// Compress linear chains of edges into a single edge.
uint64_t AssemblyGraph::compress()
{
    AssemblyGraph& assemblyGraph = *this;
    uint64_t compressCount = 0;

    // Find linear chains of 2 or more edges.
    vector< std::list<edge_descriptor> > chains;
    findLinearChains(assemblyGraph, 2, chains);

    for(const auto& chain: chains) {
        SHASTA2_ASSERT(chain.size() > 1);

        // Get the first and last edge of this chain.
        const edge_descriptor e0 = chain.front();
        const edge_descriptor e1 = chain.back();

        // Get the first and last edge of this chain.
        const vertex_descriptor v0 = source(e0, assemblyGraph);
        const vertex_descriptor v1 = target(e1, assemblyGraph);

        // Add the new edge.
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(v0, v1, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
        AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];

        // Concatenate the steps of all the edges in the chain.
        for(const edge_descriptor e: chain) {
            const AssemblyGraphEdge& edge = assemblyGraph[e];
            copy(edge.begin(), edge.end(), back_inserter(edgeNew));
        }



        // Minimal debug output.
        if(compressDebugLevel >= 1) {
            cout << "Compress " << assemblyGraph[chain.front()].id << "..." <<
                assemblyGraph[chain.back()].id <<
                " into " << edgeNew.id << endl;
        }



        // Compact debug output.
        if(compressDebugLevel >= 2) {
            cout << "Compress";
            for(const edge_descriptor e: chain) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << " into " << edgeNew.id << endl;
        }



        // Detailed debug output.
        if(compressDebugLevel >= 3) {
            uint64_t stepCount = 0;
            for(const edge_descriptor e: chain) {
                const AssemblyGraphEdge& edge = assemblyGraph[e];
                const uint64_t stepBegin = stepCount;
                const uint64_t stepEnd = stepBegin + edge.size();
                cout << edge.id << " " << edge.size() << " steps become " <<
                    edgeNew.id << " steps " << stepBegin << "-" << stepEnd << endl;
                stepCount = stepEnd;
            }
        }



        // Now we can remove the edges of the chain and its internal vertices.
        bool isFirst = true;
        for(const edge_descriptor e: chain) {
            if(isFirst) {
                isFirst = false;
            } else {
                const vertex_descriptor v = source(e, assemblyGraph);
                boost::clear_vertex(v, assemblyGraph);
                boost::remove_vertex(v, assemblyGraph);
            }
        }

        ++compressCount;
    }

    return compressCount;
}



void AssemblyGraph::save(ostream& s) const
{
    boost::archive::binary_oarchive archive(s);
    archive << *this;
}



void AssemblyGraph::load(istream& s)
{
    boost::archive::binary_iarchive archive(s);
    archive >> *this;
}



void AssemblyGraph::save(const string& stage) const
{
    // If not using persistent binary data, do nothing.
    if(largeDataFileNamePrefix.empty()) {
        return;
    }

    // First save to a string.
    std::ostringstream s;
    save(s);
    const string dataString = s.str();

    // Now save the string to binary data.
    const string name = largeDataName("AssemblyGraph-" + stage);
    MemoryMapped::Vector<char> data;
    data.createNew(name, largeDataPageSize);
    data.resize(dataString.size());
    const char* begin = dataString.data();
    const char* end = begin + dataString.size();
    copy(begin, end, data.begin());
}



void AssemblyGraph::load(const string& assemblyStage)
{
    // Access the binary data.
    MemoryMapped::Vector<char> data;
    try {
        const string name = largeDataName("AssemblyGraph-" + assemblyStage);
        data.accessExistingReadOnly(name);
    } catch (std::exception&) {
        throw runtime_error("Assembly graph at stage " + assemblyStage +
            " is not available.");
    }
    const string dataString(data.begin(), data.size());

    // Load it from here.
    std::istringstream s(dataString);
    try {
        load(s);
    } catch(std::exception& e) {
        throw runtime_error("Error reading assembly graph at stage " + assemblyStage +
            ": " + e.what());
    }
}



// Write a csv file that can be loaded in Bandage.
void AssemblyGraph::writeCsv(const string& fileName) const
{
    ofstream csv(fileName);
    writeCsv(csv);
}



// Write a csv file that can be loaded in Bandage.
void AssemblyGraph::writeCsv(ostream& csv) const
{
    const AssemblyGraph& assemblyGraph = *this;

    csv << "Segment,Number of steps,Average coverage,Estimated length,Actual length\n";
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        const uint64_t coverage = uint64_t(std::round(edge.averageCoverage()));
        csv <<
            edge.id << "," <<
            edge.size() << "," <<
            coverage << "," <<
            edge.offset() << ",";
        if(edge.wasAssembled) {
            csv << edge.sequenceLength();
        }
        csv << ",\n";
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



uint64_t AssemblyGraph::phaseSuperbubbleChains()
{
    performanceLog << timestamp << "AssemblyGraph::phaseSuperbubbleChains begins." << endl;

    PhaseSuperbubbleChainsData& data = phaseSuperbubbleChainsData;
    data.superbubbleChains = make_shared< vector<SuperbubbleChain> >();
    vector<SuperbubbleChain>& superbubbleChains = *(data.superbubbleChains);
    data.totalChangeCount = 0;

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
    data.superbubbleChains = 0;
    uint64_t changeCount = data.totalChangeCount;

    changeCount += compress();
    performanceLog << timestamp << "AssemblyGraph::phaseSuperbubbleChains ends." << endl;

    return changeCount;
}



void AssemblyGraph::phaseSuperbubbleChainsThreadFunction([[maybe_unused]] uint64_t threadId)
{
    PhaseSuperbubbleChainsData& data = phaseSuperbubbleChainsData;
    vector<SuperbubbleChain>& superbubbleChains = *(data.superbubbleChains);

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all superbubble chains assigned to this batch.
        for(uint64_t superbubbleChainId=begin; superbubbleChainId<end; superbubbleChainId++) {
            SuperbubbleChain& superbubbleChain = superbubbleChains[superbubbleChainId];
            const uint64_t changeCount = superbubbleChain.phase1(
                *this,
                superbubbleChainId);
            __sync_fetch_and_add(&data.totalChangeCount, changeCount);
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



// Find the non-trivial strongly connected components.
// Each component is stored with vertices sorted to permit binary searches.
void AssemblyGraph::findStrongComponents(
    vector< vector<vertex_descriptor> >& strongComponents) const
{
    const AssemblyGraph& assemblyGraph = *this;

    // Map the vertices to integers.
    uint64_t vertexIndex = 0;
    std::map<vertex_descriptor, uint64_t> vertexMap;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexMap.insert({v, vertexIndex++});
    }

    // Compute strong components.
    std::map<vertex_descriptor, uint64_t> componentMap;
    boost::strong_components(
        assemblyGraph,
        boost::make_assoc_property_map(componentMap),
        boost::vertex_index_map(boost::make_assoc_property_map(vertexMap)));

    // Gather the vertices in each strong component.
    std::map<uint64_t, vector<vertex_descriptor> > componentVertices;
    for(const auto& p: componentMap) {
        componentVertices[p.second].push_back(p.first);
    }

    strongComponents.clear();
    for(const auto& p: componentVertices) {
        const vector<vertex_descriptor>& component = p.second;
        if(component.size() > 1) {
            strongComponents.push_back(component);
            sort(strongComponents.back().begin(), strongComponents.back().end());
        }
    }


}



// This creates a csv file that can be loaded in bandage to see
// the strongly connected components.
void AssemblyGraph::colorStrongComponents() const
{
    const AssemblyGraph& assemblyGraph = *this;

    vector< vector<vertex_descriptor> > strongComponents;
    findStrongComponents(strongComponents);

    ofstream csv("StrongComponents.csv");
    csv << "Id,Color,Component\n";

    // Loop over the non-trivial strongly connected components.
    for(uint64_t i=0; i<strongComponents.size(); i++) {
        const vector<vertex_descriptor>& strongComponent = strongComponents[i];

        for(const vertex_descriptor v0: strongComponent) {
            BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
                const vertex_descriptor v1 = target(e, assemblyGraph);
                if(std::binary_search(strongComponent.begin(), strongComponent.end(), v1)) {
                    csv << assemblyGraph[e].id << ",Green," << i << "\n";
                }
            }
        }

    }


}




double AssemblyGraphEdge::averageCoverage() const
{
    uint64_t sum = 0;
    for(const AssemblyGraphEdgeStep& step: *this) {
        sum += step.anchorPair.orientedReadIds.size();
    }

    return double(sum) / double(size());
}



double AssemblyGraphEdge::lengthWeightedAverageCoverage() const
{
    uint64_t sum0 = 0;
    uint64_t sum1 = 0;
    for(const AssemblyGraphEdgeStep& step: *this) {
        const uint64_t length = wasAssembled ? step.sequence.size() : step.offset;
        const uint64_t coverage = step.anchorPair.orientedReadIds.size();
        sum0 += length;
        sum1 += length * coverage;
    }

    return double(sum1) / double(sum0);
}



// Compute oriented read journeys in the AssemblyGraph.
void AssemblyGraph::computeJourneys()
{
    const uint64_t orientedReadCount = journeys.size();

    AssemblyGraph& assemblyGraph = *this;


    // Compute uncompressed journeys for all oriented reads.
    class AssemblyGraphJourneyEntry {
    public:
        edge_descriptor e;
        uint64_t stepId;

        // These are positions in the standard Journeys on the Anchors.
        uint64_t positionInJourneyA;
        uint64_t positionInJourneyB;

        bool operator<(const AssemblyGraphJourneyEntry& that) const {
            return positionInJourneyA + positionInJourneyB <
                that.positionInJourneyA + that.positionInJourneyB;
        }
    };
    vector< vector<AssemblyGraphJourneyEntry> > assemblyGraphJourneys(orientedReadCount);

    // Loop over AssemblyGraph edges.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        // Loop over steps of this edge.
        for(uint64_t stepId=0; stepId<edge.size(); stepId++) {
            const AssemblyGraphEdgeStep& step = edge[stepId];

            // Locate the anchors for this step.
            const AnchorPair& anchorPair = step.anchorPair;
            const AnchorId anchorIdA = anchorPair.anchorIdA;
            const AnchorId anchorIdB = anchorPair.anchorIdB;
            const Anchor anchorA = anchors[anchorIdA];
            const Anchor anchorB = anchors[anchorIdB];

            // Loop over OrientedReadIds of this step.
            auto itA = anchorA.begin();
            auto itB = anchorB.begin();
            for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {

                // Locate this OrientedReadId in the two anchors.
                for(; (itA != anchorA.end()) and (itA->orientedReadId != orientedReadId); ++itA) {}
                SHASTA2_ASSERT(itA != anchorA.end());
                SHASTA2_ASSERT(itA->orientedReadId == orientedReadId);
                const AnchorMarkerInfo& infoA = *itA;
                for(; (itB != anchorB.end()) and (itB->orientedReadId != orientedReadId); ++itB) {}
                SHASTA2_ASSERT(itB != anchorB.end());
                const AnchorMarkerInfo& infoB = *itB;
                SHASTA2_ASSERT(itB->orientedReadId == orientedReadId);

                const uint32_t positionInJourneyA = infoA.positionInJourney;
                const uint32_t positionInJourneyB = infoB.positionInJourney;

                AssemblyGraphJourneyEntry entry;
                entry.e = e;
                entry.stepId = stepId;
                entry.positionInJourneyA = positionInJourneyA;
                entry.positionInJourneyB = positionInJourneyB;

                assemblyGraphJourneys[orientedReadId.getValue()].push_back(entry);

            }

        }
    }


    // Sort the journeys.
    for(vector<AssemblyGraphJourneyEntry>& v: assemblyGraphJourneys) {
        sort(v.begin(), v.end());
    }



    // Create the compressed journeys.
    // Here, we only consider transitions between AssemblyGraph edges.
    // Compressed journeys consisting of only one edge are considered empty.
    compressedJourneys.clear();
    compressedJourneys.resize(orientedReadCount);
    for(ReadId orientedReadIdValue=0; orientedReadIdValue<orientedReadCount; orientedReadIdValue++) {
        const vector<AssemblyGraphJourneyEntry>& assemblyGraphJourney = assemblyGraphJourneys[orientedReadIdValue];
        vector<edge_descriptor>& compressedJourney = compressedJourneys[orientedReadIdValue];

        for(uint64_t i1=0; i1<assemblyGraphJourney.size(); i1++) {
            const AssemblyGraphJourneyEntry& entry1 = assemblyGraphJourney[i1];
            const edge_descriptor e1 = entry1.e;

            if(i1 == 0) {
                compressedJourney.push_back(e1);
            } else {
                const uint64_t i0 = i1 - 1;
                const AssemblyGraphJourneyEntry& entry0 = assemblyGraphJourney[i0];
                const edge_descriptor e0 = entry0.e;
                if(e1 != e0) {
                    compressedJourney.push_back(e1);
                }
            }
        }
        if(compressedJourney.size() == 1) {
            compressedJourney.clear();
        }
    }



    // Write out the compressed journeys.
    {
        ofstream csv("AssemblyGraphCompressedJourneys.csv");
        for(ReadId orientedReadIdValue=0; orientedReadIdValue<orientedReadCount; orientedReadIdValue++) {
            const OrientedReadId orientedReadId = OrientedReadId::fromValue(orientedReadIdValue);
            const vector<edge_descriptor>& compressedJourney = compressedJourneys[orientedReadIdValue];
            csv << orientedReadId << ",";
            for(const edge_descriptor e: compressedJourney) {
                csv << assemblyGraph[e].id << ",";
            }
            csv << "\n";
        }
    }

    {
        ofstream csv("AssemblyGraphJourneys.csv");
        csv << "OrientedReadId,Segment,Step,PositionInJourneyA,PositionInJourneyB\n";
        for(ReadId orientedReadIdValue=0; orientedReadIdValue<orientedReadCount; orientedReadIdValue++) {
            const OrientedReadId orientedReadId = OrientedReadId::fromValue(orientedReadIdValue);
            const vector<AssemblyGraphJourneyEntry>& assemblyGraphJourney = assemblyGraphJourneys[orientedReadIdValue];
            for(const AssemblyGraphJourneyEntry& entry: assemblyGraphJourney) {
                csv << orientedReadId << ",";
                csv << assemblyGraph[entry.e].id << ",";
                csv << entry.stepId << ",";
                csv << entry.positionInJourneyA << ",";
                csv << entry.positionInJourneyB << "\n";
            }
        }
    }
}



// The detangling process can generate empty edges (edges without steps).
// This removes them by collapsing the vertices they join.
void AssemblyGraph::removeEmptyEdges()
{
    AssemblyGraph& assemblyGraph = *this;
    const bool debug = false;

    // We need to find groups of vertices that need to be collapsed together.
    // Usually it will be just two vertices to be collapsed together,
    // but it could be larger groups if there are adjacent empty edges.
    // So we need to find the connected components generated by the empty edges.

    // Map vertices to integer.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    vector<vertex_descriptor> vertexTable;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexIndexMap.insert(make_pair(v, vertexIndex++));
        vertexTable.push_back(v);
    }

    // Initialize the disjoint sets data structure.
    const uint64_t n = vertexIndexMap.size();
    // Initialize the disjoint sets data structure.
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }

    // Loop over the empty edges.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        if(assemblyGraph[e].empty()) {
            const vertex_descriptor v0 = source(e, assemblyGraph);
            const vertex_descriptor v1 = target(e, assemblyGraph);
            const uint64_t i0 = vertexIndexMap[v0];
            const uint64_t i1 = vertexIndexMap[v1];
            disjointSets.union_set(i0, i1);
            edgesToBeRemoved.push_back(e);
        }
    }

    // Gather the vertices in each connected component.
    vector< vector<uint64_t> > groups(n);
    for(uint64_t i=0; i<n; i++) {
        groups[disjointSets.find_set(i)].push_back(i);
    }

    // Collapse each group into a single vertex.
    for(const vector<uint64_t>& group: groups) {
        if(group.size() < 2) {
            continue;
        }
        if(debug) {
            cout << "Found a group of " << group.size() << " vertices to be collapsed:";
            for(const uint64_t vertexIndex: group) {
                const vertex_descriptor v = vertexTable[vertexIndex];
                cout << " " << assemblyGraph[v].id << "/" << anchorIdToString(assemblyGraph[v].anchorId);
            }
            cout << endl;
        }

        // Create the collapsed vertex.
        const AnchorId anchorId = assemblyGraph[vertexTable[group.front()]].anchorId;
        const vertex_descriptor vNew = add_vertex(
             AssemblyGraphVertex(anchorId, assemblyGraph.nextVertexId++), assemblyGraph);



        // For each vertex in this group, reroute all incoming/outgoing non-empty edges
        // to/from the new vertex.
        for(const uint64_t vertexIndex: group) {
            const vertex_descriptor v = vertexTable[vertexIndex];

            // Loop over non-empty incoming edges.
            BGL_FORALL_INEDGES(v, eOld, assemblyGraph, AssemblyGraph) {
                AssemblyGraphEdge& oldEdge = assemblyGraph[eOld];
                if(oldEdge.empty()) {
                    continue;
                }
                const vertex_descriptor u = source(eOld, assemblyGraph);

                // Create the new edge, with target vNew.
                AssemblyGraph::edge_descriptor eNew;
                tie(eNew, ignore) = add_edge(u, vNew, assemblyGraph);
                AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
                newEdge.id = oldEdge.id;
                newEdge.swapSteps(oldEdge);

                edgesToBeRemoved.push_back(eOld);
            }


            // Loop over non-empty outgoing edges.
            BGL_FORALL_OUTEDGES(v, eOld, assemblyGraph, AssemblyGraph) {
                AssemblyGraphEdge& oldEdge = assemblyGraph[eOld];
                if(oldEdge.empty()) {
                    continue;
                }
                const vertex_descriptor u = target(eOld, assemblyGraph);

                // Create the new edge, with source vNew.
                AssemblyGraph::edge_descriptor eNew;
                tie(eNew, ignore) = add_edge(vNew, u, assemblyGraph);
                AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
                newEdge.id = oldEdge.id;
                newEdge.swapSteps(oldEdge);

                edgesToBeRemoved.push_back(eOld);
            }
        }
    }

    deduplicate(edgesToBeRemoved);
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, assemblyGraph);
    }

    // Remove any vertices that were left isolated.
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        if((in_degree(v, assemblyGraph) == 0) and (out_degree(v, assemblyGraph) == 0)) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        boost::remove_vertex(v, assemblyGraph);
    }
}



// Return true if the two specified edges can be connected for assembly.
// If necessary, this constructs the RestrictedAnchorGraph between the two edges
// and checks that all is good.
bool AssemblyGraph::canConnect(edge_descriptor e0, edge_descriptor e1) const
{
    const AssemblyGraph& assemblyGraph = *this;

    // If the last AnchorId of e0 is the same as the first AnchorId of v1,
    // They can be trivially connected.
    const AssemblyGraph::vertex_descriptor v0 = target(e0, assemblyGraph);
    const AssemblyGraph::vertex_descriptor v1 = source(e1, assemblyGraph);
    if(v0 == v1) {
        return true;
    }
    const AnchorId anchorId0 = assemblyGraph[v0].anchorId;
    const AnchorId anchorId1 = assemblyGraph[v1].anchorId;
    if(anchorId0 == anchorId1) {
        return true;
    }



    // In the general case, we need to construct the RestrictedAnchorGraph
    // between e0 and e1. The RestrictedAnchorGraph requires a TangleMatrix,
    // so we create a 1x1 tangle matrix with e0 as the entrance and e1 as the exit.
    ostream html(0);
    TangleMatrix1 tangleMatrix(assemblyGraph, {e0}, {e1}, html);



    // Create the RestrictedAnchorGraph. If this fails, return false.
    try {
        ostream html(0);

        // Create the RestrictedAnchorGraph.
        RestrictedAnchorGraph restrictedAnchorGraph(anchors, journeys, tangleMatrix, 0, 0, html);
        vector<RestrictedAnchorGraph::edge_descriptor> longestPath;

        // Check that we have a path.
        restrictedAnchorGraph.findOptimalPath(anchorId0, anchorId1, longestPath);

        // Check coverage.
        uint64_t minCoverage = std::numeric_limits<uint64_t>::max();
        for(const RestrictedAnchorGraph::edge_descriptor e: longestPath) {
            const RestrictedAnchorGraphEdge& edge = restrictedAnchorGraph[e];
            minCoverage = min(minCoverage, edge.anchorPair.size());
        }
        if(minCoverage == 0) {
            return false;
        }

    } catch (const std::exception& e) {
        return false;
    }



    // If we get here, all is good.
    return true;
}



void AssemblyGraph::removeIsolatedVertices()
{
    AssemblyGraph& assemblyGraph = *this;


    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        if((out_degree(v, assemblyGraph) == 0) and (in_degree(v, assemblyGraph) == 0)) {
            verticesToBeRemoved.push_back(v);
        }
    }

    for(const vertex_descriptor v: verticesToBeRemoved) {
        boost::remove_vertex(v, assemblyGraph);
    }

}


// Remove connected components with a low N50.
void AssemblyGraph::removeLowN50Components([[maybe_unused]] uint64_t minN50)
{
    AssemblyGraph& assemblyGraph = *this;
    const bool debug = false;

    // Map the vertices to integers.
    uint64_t vertexIndex = 0;
    vector<vertex_descriptor> vertexTable;                  // Map integers to vertices
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;   // Map vertices to integers
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexTable.push_back(v);
        vertexIndexMap.insert({v, vertexIndex++});
    }



    // Compute connected components.
    // We can't use boost::connected_components because that only
    // supports undirected graphs.
    DisjointSets disjointSets(vertexIndexMap.size());
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);
        disjointSets.unionSet(vertexIndexMap[v0], vertexIndexMap[v1]);
    }
    vector< vector<uint64_t> > components;
    disjointSets.gatherComponents(1, components);
    if(debug) {
        cout << "Found " << components.size() << " connected components of the AssemblyGraph." << endl;
    }



    // Loop over connected components.
    for(const vector<uint64_t>& component: components) {

        // Gather the lengths of all edges in this component.
        vector<uint64_t> lengths;
        for(uint64_t vIndex0: component) {
            const vertex_descriptor v0 = vertexTable[vIndex0];
            BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
                const vertex_descriptor v1 = target(e, assemblyGraph);
                const uint64_t vIndex1 = vertexIndexMap[v1];
                SHASTA2_ASSERT(disjointSets.findSet(vIndex0) == disjointSets.findSet(vIndex1));
                lengths.push_back(assemblyGraph[e].length());
            }
        }
        const uint64_t totalLength = std::accumulate(lengths.begin(), lengths.end(), 0);

        // Find the N50.
        std::ranges::sort(lengths, std::greater<uint64_t>());
        uint64_t cumulativeLength = 0;
        uint64_t n50 = 0;
        for(const uint64_t length: lengths) {
            cumulativeLength += length;
            if(2 * cumulativeLength >= totalLength) {
                n50 = length;
                break;
            }
        }


        // Decide if we keep this component.
        const bool keep = (n50 >= minN50);
        if(debug) {
            if(keep) {
                cout << "Keeping";
            } else {
                cout << "Discarding";
            }
            cout << " a connected component with " << lengths.size() <<
                " segments, total length " << totalLength <<
                ", N50 " << n50 << endl;
        }

        // If not keeping this component, remove its vertices and edges.
        if(not keep) {
            for(uint64_t vIndex: component) {
                const vertex_descriptor v = vertexTable[vIndex];
                boost::clear_vertex(v, assemblyGraph);
                boost::remove_vertex(v, assemblyGraph);
            }

        }

    }

}
