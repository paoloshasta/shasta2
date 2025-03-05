// Shasta
#include "Assembler.hpp"
#include "mode3-LocalAssembly.hpp"
#include "Mode3Assembler.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"
#include <map>



void Assembler::accessMode3Assembler()
{
    shared_ptr<mode3::Anchors> anchorsPointer =
        make_shared<mode3::Anchors>(MappedMemoryOwner(*this), reads(), assemblerInfo->k, markers());
    mode3Assembler = make_shared<Mode3Assembler>(*this,
        assemblerInfo->k, reads(), markers(),
        anchorsPointer, httpServerData.assemblerOptions->assemblyOptions.mode3Options);
}



void Assembler::exploreLocalAnchorGraph(const vector<string>& request, ostream& html)
{
    mode3Assembler->exploreLocalAnchorGraph(request, html);
}
