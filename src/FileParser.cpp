/*
 * FileParser.cpp
 *
 *  Created on: 4 Jan 2015
 *      Author: helenginn
 */

#include "FileParser.h"

#include "FileReader.h"

#include <sstream>
#include <vector>
#include <iostream>
#include <locale>
#include <stdio.h>

#define FILE_PARSER_CPP_

ParametersMap FileParser::parameters;
std::ostringstream FileParser::log;
int FileParser::threadsFound = 0;

bool FileParser::hasKey(std::string key)
{
	return (parameters.count(key) > 0);
}

void FileParser::simpleFloat(ParametersMap *map, std::string command,
		std::string rest)
{
	double theFloat = atof(rest.c_str());
    
	log << "Setting double " << command << " to " << theFloat << std::endl;

	(*map)[command] = theFloat;
}

void FileParser::simpleBool(ParametersMap *map, std::string command,
		std::string rest)
{
	bool on = (rest == "ON" || rest == "on" ? 1 : 0);

	log << "Setting bool " << command << " to " << on << std::endl;

	(*map)[command] = on;
}

void FileParser::simpleString(ParametersMap *map, std::string command,
		std::string rest)
{
	(*map)[command] = rest;

	log << "Setting string " << command << " to " << rest << std::endl;

}

void FileParser::simpleInt(ParametersMap *map, std::string command,
		std::string rest)
{
	int theInt = atoi(rest.c_str());

	log << "Setting int " << command << " to " << theInt << std::endl;

	(*map)[command] = theInt;
}

void FileParser::doubleVector(ParametersMap *map, std::string command,
		std::string rest)
{
	vector<std::string> components = FileReader::split(rest, ' ');
	vector<double> doubleVector;

	log << "Setting " << command << " to ";

	for (int i = 0; i < components.size(); i++)
	{
		double theFloat = atof(components[i].c_str());
		log << theFloat << " ";
		doubleVector.push_back(theFloat);
	}

	log << std::endl;

	(*map)[command] = doubleVector;
}

void FileParser::intVector(ParametersMap *map, std::string command,
		std::string rest)
{
	vector<std::string> components = FileReader::split(rest, ' ');
	vector<int> intVector;

	log << "Setting " << command << " to ";

	for (int i = 0; i < components.size(); i++)
	{
		int theInt = atoi(components[i].c_str());
		log << theInt << " ";
		intVector.push_back(theInt);
	}

	log << std::endl;

	(*map)[command] = intVector;
}

void FileParser::generateFunctionList()
{
	parserMap = ParserMap();

    parserMap["VERBOSITY_LEVEL"] = simpleInt;
    
	// Refinement parameters
	
    parserMap["SPACE_GROUP"] = simpleInt; // override, should be untwinned
    parserMap["UNIT_CELL"] = doubleVector; // override
    parserMap["POINT_GROUP"] = simpleInt;
    parserMap["LABIN_AMPLITUDES"] = simpleString;
    parserMap["LABIN_SIGMAS"] = simpleString;
    parserMap["LABIN_PHASES"] = simpleString;
    parserMap["LABOUT_CALC_AMPLITUDES"] = simpleString;
    parserMap["LOW_RESOLUTION"] = simpleFloat;
    parserMap["HIGH_RESOLUTION"] = simpleFloat;
    parserMap["GRIDDING"] = intVector;
    
    parserMap["WANG_ENVELOPE_OUTPUT"] = simpleString;
    parserMap["NCS_DEFINITIONS_FILE"] = simpleString;
    parserMap["TWINNED_MTZ"] = simpleString;
    parserMap["ORIGIN_IN_UNIT_CELL"] = doubleVector;
    parserMap["VIRION_GENEROUS_DIAMETER"] = simpleFloat;
    parserMap["VIRION_INNER_DIAMETER"] = simpleFloat;
    parserMap["VIRION_OUTER_DIAMETER"] = simpleFloat;
    
    parserMap["MAXIMUM_CYCLES"] = simpleInt;
    parserMap["2_FOBS_MINUS_FC"] = simpleBool;
    
    parserMap["TEMP_AVERAGE_MAP"] = simpleString;
    parserMap["TEMP_AVERAGE_MTZ"] = simpleString;
    parserMap["TEMP_AVERAGE_SORTED_MTZ"] = simpleString;
    parserMap["TEMP_FSMELT_MTZ"] = simpleString;
    
    parserMap["SKIP_FIRST_CYCLE"] = simpleBool;
    parserMap["RESUME_FROM_CYCLE"] = simpleInt;
    
	// Indexing parameters;
}

ParserFunction FileParser::splitLine(std::string line, std::string &command,
		std::string &rest)
{
	int space_index = (int)line.find_first_of(" ");

	command = line.substr(0, space_index);

	std::ostringstream stream;

	std::locale theLocale;
	for (std::string::size_type j = 0; j < command.length(); ++j)
		stream << std::toupper(command[j], theLocale);

	std::string upperCommand = stream.str();

	rest = line.substr(space_index + 1, std::string::npos);

    if (parserMap.count(upperCommand) == 0 && upperCommand != "PANEL")
    {
        std::cout << "Error: do not understand command " << upperCommand << std::endl;
        exit(1);
    }
    
    ParserFunction function = parserMap[upperCommand];
	command = upperCommand;

	return function;
}

bool FileParser::checkSpaces(std::string line)
{
	int space_index = (int)line.find_first_of(" ");

	if (space_index == std::string::npos)
	{
		log << "Warning: " << line << " has no assignment" << std::endl;
		return false;
	}

	return true;
}

FileParser::FileParser(void)
{
    
}

FileParser::FileParser(std::string name)
{
	std::cout << "Initialising parser" << std::endl;
    
	this->filename = name;
	generateFunctionList();
}

FileParser::~FileParser()
{
	// TODO Auto-generated destructor stub
}

