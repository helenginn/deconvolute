//
//  Deconvoluter.h
//  deconvolute
//
//  Created by Helen Ginn on 03/08/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __deconvolute__Deconvoluter__
#define __deconvolute__Deconvoluter__

#include <stdio.h>
#include "parameters.h"


class Deconvoluter
{
private:
    std::string scriptCentreOnVirion();
    std::string outputRoot;
    std::string wangEnvName;
    std::string molUnitName;
    std::string tempFftMapName;
    std::string tempAveMapName;
    std::string tempMtzName;
    std::string tempSortedMtzName;
    std::string cadMtzName;
    std::string fsmeltMtzName;
    MtzPtr originalMtz;
    void fastFourierTransform(std::string filename, std::string tempFftMapName, bool calcAmps);
    void gapEnvelope(std::string inputMap);
    void gapAverage(std::string inputMap, std::string inputWang, std::string inputMolUnit);
    void map_to_sf(std::string averageMap, std::string tempMtz);
    void sort_mtz(std::string tempMtz, std::string tempSortedMtz);
    void cad(std::string tempSortedMtz, std::string cadMtz);
    void fsmelt(std::string originalMtz, std::string tempMtz, std::string fsmeltMtz);
    
    int spaceGroup;
    std::vector<int> gridding;
    void runProgram(std::string invocation, std::string script, std::string programName);
    
public:
    Deconvoluter();
    void process();
};

#endif /* defined(__deconvolute__Deconvoluter__) */
