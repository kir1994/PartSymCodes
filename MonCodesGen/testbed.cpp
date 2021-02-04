#include <iostream>
#include <string>
#include "PartSymMonCodes.h"


int main(int argc, char** argv)
{
	if (argc == 1)
	{
		std::cerr << "Usage: m k t r \n";
		return 0;
	}

	unsigned curArg = 1;

	try {
		unsigned m = std::stoi(argv[curArg++]);
		unsigned k = std::stoi(argv[curArg++]);
		unsigned t = std::stoi(argv[curArg++]);
		unsigned r = std::stoi(argv[curArg++]);
		CPartSymMonCodeGen::Generate(m, k, t, r);
    }
	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
	}

	return 0;
}
