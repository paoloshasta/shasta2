#pragma once

#include <map>
#include "string.hpp"

// The configurationTable is a map that contains
// the definitions of shasta2 built-in assembly configurations.
// The key of the map is the configuration name.
// The name is the same as the name of the conf file,
// with the extension ".conf" removed. So for example for an
// assembly configuration defined by file shasta2/conf/xyz.conf
// the configuration name is simply "xyz".
// The value of the map is a long string containing
// the entire contents of the .conf file that defines
// the configuration, including line ends.


namespace shasta2 {
	extern const std::map<string, string> configurationTable;
}
