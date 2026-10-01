// Shasta2.
#include "AssemblyGraph.hpp"
#include "color.hpp"
#include "deduplicate.hpp"
#include "DisjointSets.hpp"
#include "html.hpp"
#include "performanceLog.hpp"
#include "StrandSeparation.hpp"
#include "StrandSeparation1.hpp"
#include "Tangle.hpp"
#include "timestamp.hpp"
using namespace shasta2;

// Standard library.
#include "fstream.hpp"



void AssemblyGraph::separateStrands(const string& debugOutputBaseName)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t maxLength = 100000;

    AssemblyGraph& assemblyGraph = *this;

    SHASTA2_ASSERT(countZeroLengthSegments() == 0);

    performanceLog << timestamp << "AssemblyGraph::separateStrands begins: " <<
        debugOutputBaseName << endl;



    // Create tangles as connected components of the AssemblyGraph
    // computed considering only short Segments.
    // The vertices of each tangle are returned sorted by id
    // so we can do binary searches in them.
    vector< vector<vertex_descriptor> > tangles;
    vector<uint64_t> tangleRc;
    createTanglesBySegmentLength(maxLength, tangles, tangleRc);



    // Gather the strand contacts.
    // The Segments of each strand contact are sorted by id
    // as required by the constructor of StrandSeparation::StrandContact.
    vector< vector<Segment> > strandContacts;
    for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
        if(tangleRc[tangleId] == tangleId) {
            vector<vertex_descriptor>& strandContactVertices = tangles[tangleId];
            vector<Segment>& strandContact = strandContacts.emplace_back();
            for(const vertex_descriptor v0: strandContactVertices) {
                BGL_FORALL_OUTEDGES(v0, segment, assemblyGraph, AssemblyGraph) {
                    const vertex_descriptor v1 = target(segment, assemblyGraph);
                    if(binary_search(strandContactVertices.begin(), strandContactVertices.end(), v1, orderById)) {
                        if(assemblyGraph[segment].length() <= maxLength){
                            strandContact.push_back(segment);
                        }
                    }
                }
            }
            sort(strandContact.begin(), strandContact.end(), orderById);
        }
    }



    // Write a csv file that can be loaded in Bandage to show the strand contacts.
    {
        ofstream csv(debugOutputBaseName + "-StrandContacts-Bandage.csv");
        csv << "Segment,StrandContact,Color\n";

        for(uint64_t strandContactId=0; strandContactId<strandContacts.size(); strandContactId++) {
            const string color = randomHslColor(strandContactId, 0.75, 0.5);
            const vector<Segment>& strandContact = strandContacts[strandContactId];

            for(const Segment segment: strandContact) {
                csv << id(segment) << ",";
                csv << strandContactId << ",";
                csv << color << "\n";
            }
        }
    }


    // Process each strand contact separately.
    for(uint64_t strandContactId=0; strandContactId<strandContacts.size(); strandContactId++) {
        const StrandSeparation::StrandContact strandContact(
            assemblyGraph,
            strandContacts[strandContactId],
            debugOutputBaseName,
            strandContactId);
    }


    strandSymmetricCompress();

    performanceLog << timestamp << "AssemblyGraph::separateStrands ends: " <<
        debugOutputBaseName << endl;
}



void AssemblyGraph::separateStrands1(const string& debugOutputBaseName)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t maxLength = 100000;

    AssemblyGraph& assemblyGraph = *this;

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

    // Write a csv file that can be imported into Bandage to see
    // the tangles.
    vector<uint64_t> strandContactsRc(strandContacts.size());
    std::iota(strandContactsRc.begin(), strandContactsRc.end(), 0);
    writeTangles(strandContacts, strandContactsRc,
        debugOutputBaseName + "-StrandContacts-Bandage.csv");

    // Process each strand contact separately.
    for(uint64_t strandContactId=0; strandContactId<strandContacts.size(); strandContactId++) {
        const StrandSeparation1::StrandContact strandContact(
            assemblyGraph,
            strandContacts[strandContactId],
            debugOutputBaseName,
            strandContactId);
    }

    removeIsolatedVertices();

    performanceLog << timestamp << "AssemblyGraph::separateStrands1 ends: " <<
        debugOutputBaseName << endl;
}
