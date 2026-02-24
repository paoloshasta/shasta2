#pragma once

// In the AssemblyGraph, each edge is a linear chains of adjacent AnchorPairs.
// It is initially created from linear chains of edges in the AnchorGraph.
// Each AssemblyGraphEdge generates a gfa segment.

// Shasta.
#include "AnchorPair.hpp"
#include "Base.hpp"
#include "invalid.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/serialization/vector.hpp>

// Standard library.
#include "memory.hpp"
#include "tuple.hpp"


namespace shasta2 {

    class AssemblyGraph;
    class AssemblyGraphVertex;
    class AssemblyGraphEdge;
    class AssemblyGraphEdgeStep;

    using AssemblyGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        AssemblyGraphVertex,
        AssemblyGraphEdge>;

    class AnchorGraph;
    class Anchors;
    class Assembler;
    class Options;
    class Detangler;
    class Superbubble;
    class SuperbubbleChain;

}



// The AnchorId of a vertex is the last AnchorId of the last AnchorPair
// of all incoming edges and the first AnchorId of the first AnchorPair
// of each outgoing edge.
// When the AssemblyGraph is initially created from the AnchorGraph,
// there can be at most one vertex for each AnchorId.
// However that is no longer true after detangling.
class shasta2::AssemblyGraphVertex {
public:
    AnchorId anchorId = invalid<AnchorId>;
    uint64_t id = invalid<uint64_t>;

    AssemblyGraphVertex(AnchorId anchorId, uint64_t id) :
        anchorId(anchorId), id(id) {}
    AssemblyGraphVertex() {}

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorId;
        ar & id;
    }
};



class shasta2::AssemblyGraphEdgeStep {
public:
    AnchorPair anchorPair;
    uint64_t offset = invalid<uint64_t>;
    vector<Base> sequence;

    AssemblyGraphEdgeStep(
        const AnchorPair& anchorPair,
        uint64_t offset) :
        anchorPair(anchorPair),
        offset(offset)
    {}

    AssemblyGraphEdgeStep()
    {}

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorPair;
        ar & offset;
        ar & sequence;
    }
};



class shasta2::AssemblyGraphEdge : public vector<AssemblyGraphEdgeStep> {
public:
    uint64_t id = invalid<uint64_t>;
    bool wasAssembled = false;

    AssemblyGraphEdge(uint64_t id = invalid<uint64_t>) : id(id) {}

    void check(const Anchors&) const;

    uint64_t offset() const;
    uint64_t sequenceLength() const;
    void getSequence(vector<Base>&) const;

    // If wasAssembled is set, this returns sequenceLength().
    // Otherwise, it returns offset().
    uint64_t length() const;

    double averageCoverage() const;
    double lengthWeightedAverageCoverage() const;

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object< vector<AssemblyGraphEdgeStep> >(*this);
        ar & id;
        ar & wasAssembled;
    }

    void swapSteps(AssemblyGraphEdge& that)
    {
        vector<AssemblyGraphEdgeStep>::swap(that);
    }

    AnchorId firstAnchorId() const
    {
        return front().anchorPair.anchorIdA;
    }
    AnchorId lastAnchorId() const
    {
        return back().anchorPair.anchorIdB;
    }
};



class shasta2::AssemblyGraph :
    public AssemblyGraphBaseClass,
    public MappedMemoryOwner,
    public MultithreadedObject<AssemblyGraph> {
public:

    // Initial construction from the AnchorGraph.
    AssemblyGraph(
        const Anchors&,
        const Journeys&,
        const AnchorGraph&,
        const Options&);

    // Deserialize constructor.
    AssemblyGraph(
        const Anchors&,
        const Journeys&,
        const Options&,
        const string& stage);

    // Detangle, phase, assemble sequence, output.
    void simplifyAndAssemble();

    const Anchors& anchors;
    const Journeys& journeys;
    uint64_t nextVertexId = 0;
    uint64_t nextEdgeId = 0;
    const Options& options;

private:

    void check() const;

    // Bubble cleanup.
    // A bubble is a set of parallel edges in the AssemblyGraph.
    // The edges are sorted by id.
    class Bubble {
    public:
        vertex_descriptor v0;
        vertex_descriptor v1;
        vector<edge_descriptor> edges;
    };
    void findBubbles(vector<Bubble>&) const;
public:
    uint64_t bubbleCleanup();
    uint64_t bubbleCleanupIteration(vector< pair<vertex_descriptor, vertex_descriptor> >& excludeList);
private:
    bool bubbleCleanup(const Bubble&);


    // Analyze a Bubble and finds pairs of "similar" branches.
    // See analyzeSimilarSequences.h for more information.
    bool analyzeBubble(
        const Bubble&,
        const vector<uint64_t> minRepeatCount,
        vector< pair<uint64_t, uint64_t> >& similarPairs
        ) const;

public:

    // Compress linear chains of edges into a single edge.
    uint64_t compress();
    uint64_t compressDebugLevel = 0; // 1=minimal, 2=compact, 3=detailed.



    // Class to order vertices or edges by id.
    class OrderById {
    public:
        OrderById(const AssemblyGraph& assemblyGraph): assemblyGraph(assemblyGraph) {}
        const AssemblyGraph& assemblyGraph;
        bool operator()(vertex_descriptor x, vertex_descriptor y) const
        {
            return assemblyGraph[x].id < assemblyGraph[y].id;
        }
        bool operator()(edge_descriptor x, edge_descriptor y) const
        {
            return assemblyGraph[x].id < assemblyGraph[y].id;
        }

        // Also order pairs of edges.
        using EdgePair = pair<edge_descriptor, edge_descriptor>;
        bool operator()(const EdgePair& x,const EdgePair& y) const
        {
            if(assemblyGraph[x.first].id < assemblyGraph[y.first].id) {
                return true;
            }
            if(assemblyGraph[x.first].id > assemblyGraph[y.first].id) {
                return false;
            }
            return assemblyGraph[x.second].id < assemblyGraph[y.second].id;
        }
    };
    const OrderById orderById;



    bool hasSelfEdge(vertex_descriptor v) const
    {
        bool edgeExists = false;
        tie(ignore, edgeExists) = boost::edge(v, v, *this);
        return edgeExists;
    }



    // Return true if the two specified edges can be connected for assembly.
    // If necessary, this construct the RestrictedAnchorGraph between the two edges
    // and checks that all is good.
    bool canConnect(edge_descriptor, edge_descriptor) const;

    // The detangling process can generate empty edges (edges without steps).
    // This removes them by collapsing the vertices they join.
    void removeEmptyEdges();

    // Remove isolated vertices.
    void removeIsolatedVertices();

    // Remove connected components with a low N50.
    void removeLowN50Components(uint64_t minN50);



public:

    // Compute compressed journeys in the AssemblyGraph.
    // The compressed journey of an oriented read
    // is the sequence of assembly graph edges it visits.
    void computeJourneys();
    vector< vector<edge_descriptor> > compressedJourneys;



    // Superbubbles and SuperbubbleChains.
    // - A Superbubble is a Tangle in which all entrance edges are into
    //   a single vertex and all exit edges are from a single vertex.
    // - A SuperbubbleChain is a linear sequence of Superbubbles.
    //   In the sequence, the target vertex of a Superbubble is the same
    //   as the source vertex of the next Superbubble in the SuperbubbleChain.

    // This finds all Superbubbles seen using options.findSuperbubblesMaxDistance.
    // Some pairs of Superbubble can intersect (that is, they can have common edges).
    void findSuperbubbles(vector<Superbubble>&) const;

    // Remove Superbubbles that are entirely contained in a larger superbubble.
    void removeContainedSuperbubbles(vector<Superbubble>&) const;

    void writeSuperbubbles(const vector<Superbubble>&, const string& fileName) const;
    void writeSuperbubblesForBandage(const vector<Superbubble>&, const string& fileName) const;

    void findSuperbubbleChains(
        const vector<Superbubble>&,
        vector<SuperbubbleChain>&
        ) const;

    void writeSuperbubbleChains(const vector<SuperbubbleChain>&, const string& fileName) const;
    void writeSuperbubbleChainsForBandage(const vector<SuperbubbleChain>&, const string& fileName) const;



    // Phasing of SuperbubbleChains.
    uint64_t phaseSuperbubbleChains();
    void phaseSuperbubbleChainsThreadFunction(uint64_t threadId);
    class PhaseSuperbubbleChainsData {
    public:
        shared_ptr< vector<SuperbubbleChain> > superbubbleChains;
        uint64_t totalChangeCount = 0;
    };
    PhaseSuperbubbleChainsData phaseSuperbubbleChainsData;



    // Strongly connected components.

    // Find the non-trivial strongly connected components.
    // Each component is stored with vertices sorted to permit binary searches.
    void findStrongComponents(vector< vector<vertex_descriptor> >&) const;

    // This creates a csv file that can be loaded in bandage to see
    // the strongly connected components.
    void colorStrongComponents() const;



    // Read following.
    // Note assemlbyPaths are not necessarily paths in the AssemblyGraph.
    // There may be jumps, which are bridged using local assemblies.
    void findAssemblyPaths(vector< vector<edge_descriptor> >& assemblyPaths) const;
    void connectAssemblyPaths(const vector< vector<edge_descriptor> >&  assemblyPaths);
    uint64_t findAndConnectAssemblyPaths();



    // Output.
    void write(const string& stage);
    void writeIntermediateStageIfRequested(const string& name);
    void writeFasta(const string& stage) const;
private:
    void writeGfa(const string& fileName) const;
    void writeGfa(ostream&) const;
    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;
    void writeCsv(const string& fileName) const;
    void writeCsv(ostream&) const;



    // Sequence assembly.

    // Assemble sequence for all edges.
public:
    void assembleAll();

    // Clear sequence from all steps of all edges.
    void clearSequence();
private:

    // Assemble sequence for the specified edge.
    void assemble(edge_descriptor);

    // Assemble sequence for step i of the specified edge.
    // This is the lowest level sequence assembly function and is not multithreaded.
    // It runs a LocalAssembly2 on the AnchorPair for that step.
    void assembleStep(edge_descriptor, uint64_t i);

    // Assemble sequence for all edges in the edgesToBeAssembled vector.
    // This fills in the stepsToBeAssembled with all steps of those edges,
    // then assembles each of the steps in parallel.
    void assemble();
    void assembleThreadFunction(uint64_t threadId);
    vector<edge_descriptor> edgesToBeAssembled;
    vector< pair<edge_descriptor, uint64_t> > stepsToBeAssembled;




    // Serialization.
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object<AssemblyGraphBaseClass>(*this);
        ar & nextVertexId;
        ar & nextEdgeId;
    }
    void save(ostream&) const;
    void load(istream&);

    // These do save/load to/from mapped memory.
    // The file name is AssemblyGraph-Stage.
    void save(const string& stage) const;
    void load(const string& stage);
};
