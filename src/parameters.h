/*
 * parameters.h
 *
 *  Created on: 17 Nov 2014
 *      Author: helenginn
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#define PARTIAL_CUTOFF 0.2
#define MULTIPLIER 512
#define OFFSET 256

class Miller;


#include "parameters.h"
#include <sstream>
#include <memory>
#include <vector>
#include <map>
#include <boost/variant.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <string>

using std::vector;

typedef enum
{
	RFactorNone, RFactorTypeMerge, RFactorTypePim, RFactorTypeMeas,
} RFactorType;

class Matrix;
class MtzManager;

typedef boost::shared_ptr<Miller> MillerPtr;
typedef boost::shared_ptr<MtzManager>MtzPtr;
typedef std::shared_ptr<Matrix>MatrixPtr;

typedef boost::variant<double, double, std::string, bool, int,
		vector<double>, vector<int> > ParameterVariant;
typedef std::map<std::string, ParameterVariant> ParametersMap;
typedef void (*ParserFunction)(ParametersMap *, std::string, std::string);
typedef std::map<std::string, ParserFunction> ParserMap;

typedef double (StatisticsFunction)(MtzManager *, MtzManager *, int, int *,
		double *, double, double, bool);
typedef double (RFactorFunction)(RFactorType, MtzManager *, int *, double *,
		double, double);

typedef enum
{
    LogLevelNormal = 0,
    LogLevelDetailed = 1,
    LogLevelDebug = 2
} LogLevel;

#endif /* PARAMETERS_H_ */
