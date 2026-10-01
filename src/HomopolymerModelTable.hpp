#pragma once

#include <map>
#include "string.hpp"

// The homopolymerModelTable is a map that contains
// the definitions of shasta2 built-in HomopolymerModels.
// The key of the map is the HomopolymerModel name.
// The name is the same as the name of the csv file,
// with the suffix "HomopolymerModel-" removed
// and without the extension. So for example for a
// HomopolymerModel defined by file shasta2/conf/HomopolymerModel-xyz.csv
// the name is simply "xyz".
// The value of the map is a long string containing
// the entire contents of the csv file that defines
// the HomopolymerModel, including line ends.


namespace shasta2 {
	extern const std::map<string, string> homopolymerModelTable;
}
