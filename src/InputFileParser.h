/*
 * InputFileParser.h
 *
 *  Created on: 13 Jan 2015
 *      Author: helenginn
 */

#ifndef INPUTFILEPARSER_H_
#define INPUTFILEPARSER_H_

#include <boost/shared_ptr.hpp>
#include "FileParser.h"


class InputFileParser : public FileParser
{
private:
public:
    virtual void parse(bool fromPython = false);

	InputFileParser(std::string filename);
	virtual ~InputFileParser();
};

#endif /* INPUTFILEPARSER_H_ */
