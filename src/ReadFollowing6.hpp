#pragma once

// ReadFollowing6 code is patterned after StrandSeparation1.
// It uses a bipartite graph.

// Shasta2.
#include "AssemblyGraphBaseClass.hpp"

// Standard library.
#include "iosfwd.hpp"
#include "string.hpp"
#include "vector.hpp"

namespace shasta2 {
    namespace ReadFollowing6 {
        class Tangle;
    }
}



// Note this is unrelated to shasta2::Tangle.
class shasta2::ReadFollowing6::Tangle {
public:

    // Tangle constructor.
    // The tangleVertices must be sorted by id.
    // The Tangle must not be self-complementary - that is,
    // if it includes a vertex, it cannot also include its reverse complement.
    // The code will work on the tangle passed in and also make the
    // necessary changes on its reverse complement to keep the AssemblyGraph
    // strand-symmatric.
    // The debugOutputBaseName and tangleId are only used for debug output.
    Tangle(
        AssemblyGraph&,
        const vector<AssemblyGraphBaseClass::vertex_descriptor>& tangleVertices,
        ostream& html);
};
