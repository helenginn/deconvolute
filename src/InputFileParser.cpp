/*
 * InputFileParser.cpp
 *
 *  Created on: 13 Jan 2015
 *      Author: helenginn
 */

#include "InputFileParser.h"
#include "FileReader.h"
#include "Deconvoluter.h"
#include <sstream>

InputFileParser::InputFileParser(std::string filename) : FileParser(filename)
{
    // TODO Auto-generated constructor stub
}

InputFileParser::~InputFileParser()
{
	// TODO Auto-generated destructor stub
}

void InputFileParser::parse(bool fromPython)
{
    parameters = ParametersMap();

	std::string fileContents = FileReader::get_file_contents(filename.c_str());
	vector<std::string> fileLines = FileReader::split(fileContents, '\n');

	bool foundCommands = false;
	int continueFrom = 0;

	for (int i = 0; i < fileLines.size(); i++)
	{
		std::string line = fileLines[i];

		if (line.length() == 0)
			continue;

		if (line.substr(0, 1) == "#")
			continue;

		if (line == "COMMANDS")
		{
			log << "Found COMMANDS section" << std::endl;
			continueFrom = i;
			foundCommands = true;
			break;
		}

		if (!checkSpaces(line))
			continue;

		std::string command; std::string rest;
		ParserFunction function = this->splitLine(line, command, rest);

		if (parserMap.count(command) == 0)
		{
            std::cout << "Warning: Line \"" << line << "\" not recognised."
					<< std::endl;
			exit(0);
		}
        
		function(&parameters, command, rest);
	}
    
    Deconvoluter deconvoluter = Deconvoluter();
    deconvoluter.process();
}
