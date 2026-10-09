#include "ReadFollowing6.hpp"
using namespace shasta2;
using namespace ReadFollowing6;



ReadFollowing6::Tangle::Tangle(
    [[maybe_unused]] AssemblyGraph& assemblyGraph,
    [[maybe_unused]] const vector<AssemblyGraphBaseClass::vertex_descriptor>& tangleVertices,
    [[maybe_unused]] ostream& html
    )
{
    if(html) {
        html << "<h2>ReadFollowing6</h2>";
    }
}
