#include "AssemblyGraph.hpp"
#include "performanceLog.hpp"
#include "ReadFollowing.hpp"
#include "StrandSeparation.hpp"
#include "Tangle.hpp"
#include "timestamp.hpp"
using namespace shasta2;



bool AssemblyGraph::readFollowingStrandSymmetric5(
    uint64_t tangleId,
    const Tangle& tangle,
    ostream& html,
    const string& debugOutputBaseName)
{
    SHASTA2_ASSERT(not tangle.entrances.empty());
    SHASTA2_ASSERT(not tangle.exits.empty());

    ReadFollowing::Graph graph(*this, tangleId, tangle, html, debugOutputBaseName);

    if(graph.updateAssemblyGraph1()) {
        return true;
    } else {
        return graph.updateAssemblyGraph2();
    }
}



void AssemblyGraph::separateStrands(const string& debugOutputBaseName)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t maxLengthStart = 10000;
    const uint64_t maxLengthEnd = 100000000;

    for(uint64_t maxLength=maxLengthStart; maxLength<=maxLengthEnd; maxLength*=2) {
        separateStrands(maxLength, debugOutputBaseName);    
    }

    separateStrandsGlobal(debugOutputBaseName);
}



void AssemblyGraph::separateStrandsGlobal(const string& debugOutputBaseName)
{
    const uint64_t maxLength = std::numeric_limits<uint64_t>::max();
    separateStrands(maxLength, debugOutputBaseName);    
}



void AssemblyGraph::separateStrands(
    uint64_t maxLength, 
    const string& debugOutputBaseName)
{
    AssemblyGraph& assemblyGraph = *this;
    const bool debug = false;

    SHASTA2_ASSERT(countZeroLengthSegments() == 0);

    performanceLog << timestamp << "AssemblyGraph::separateStrands1 begins: " <<
        debugOutputBaseName << endl;

    // Create tangles as connected components of the AssemblyGraph,
    // computed considering only short Segments.
    // The vertices of each tangle are returned sorted by id
    // so we can do binary searches in them.
    vector< vector<vertex_descriptor> > tangles;
    vector<uint64_t> tangleRc;
    createTanglesBySegmentLength(maxLength, tangles, tangleRc);

    // The self-complementary tangles are strand contacts.
    vector< vector<vertex_descriptor> > strandContacts;
    for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
        if(tangleRc[tangleId] == tangleId) {
            strandContacts.push_back(tangles[tangleId]);
        }
    }

    if(maxLength == std::numeric_limits<uint64_t>::max()) {
        cout << "Strand separation at infinite maxLength " 
            " found " << strandContacts.size() << " strand contacts." << endl;
    } else {
        cout << "Strand separation at maxLength " << maxLength << 
            " found " << strandContacts.size() << " strand contacts." << endl;
    }

    if(debug) {
        // Write a csv file that can be imported into Bandage to see
        // the tangles.
        vector<uint64_t> strandContactsRc(strandContacts.size());
        std::iota(strandContactsRc.begin(), strandContactsRc.end(), 0);
        writeTangles(strandContacts, strandContactsRc,
            debugOutputBaseName + "-StrandContacts-Bandage.csv");
    }

    // Process each strand contact separately.
    for(uint64_t strandContactId=0; strandContactId<strandContacts.size(); strandContactId++) {
        const StrandSeparation::StrandContact strandContact(
            assemblyGraph,
            strandContacts[strandContactId],
            debugOutputBaseName,
            strandContactId);
    }

    removeZeroLengthSegmentsStrandSymmetric();
    removeIsolatedVertices();

    performanceLog << timestamp << "AssemblyGraph::separateStrands1 ends: " <<
        debugOutputBaseName << endl;
}
