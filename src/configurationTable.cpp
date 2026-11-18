#include "configurationTable.hpp"
using namespace shasta2;


const std::map<string, string> shasta2::configurationTable  =
{
	{
		"ONT-R10.4.1-sup5.0.0-40x-60kb-Q25-2026.10.13",
		{
			#embed "../conf/ONT-R10.4.1-sup5.0.0-40x-60kb-Q25-2026.10.13.conf"
		}
	},
		
	{
		"ONT-R10.4.1-sup5.0.0-40x-60kb-Q25-pangenome-2026.10.08",
		{
			#embed "../conf/ONT-R10.4.1-sup5.0.0-40x-60kb-Q25-pangenome-2026.10.08.conf"
		}
	}
};

