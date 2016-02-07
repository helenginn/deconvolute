//
//  Deconvoluter.cpp
//  deconvolute
//
//  Created by Helen Ginn on 03/08/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "Deconvoluter.h"
#include "FileParser.h"
#include "MtzManager.h"
#include <iostream>
#include "csymlib.h"
#include "misc.h"
#include <cstdlib>
#include <iomanip>

void Deconvoluter::runProgram(std::string invocation, std::string script, std::string programName)
{
    std::string fullInvocation = invocation + " << EOF\n" + script + "EOF\n";
    
    std::cout << "***** Running " << programName << " *****" << std::endl;
    std::cout << fullInvocation << std::endl;
    std::cout << "***** " << programName << " Output *****" << std::endl;
    
    system(fullInvocation.c_str());
    
    std::cout << "***** End of " << programName << " *****" << std::endl;
    std::cout << std::endl;
}

Deconvoluter::Deconvoluter()
{
    outputRoot = FileParser::getKey("OUTPUT_ROOT", "cycle_");
    gridding = FileParser::getKey("GRIDDING", std::vector<int>(3, 360));
    
}

void Deconvoluter::fsmelt(std::string originalMtz, std::string tempMtz, std::string fsmeltMtz)
{
    std::string fObsLab = FileParser::getKey("LABIN_AMPLITUDES", std::string("F"));
    std::string fCalcLab = FileParser::getKey("LABOUT_CALC_AMPLITUDES", std::string("FC"));
    std::string sigFLab = FileParser::getKey("LABIN_SIGMAS", std::string("SIGF"));
    std::string phiLab = FileParser::getKey("LABIN_PHASES", std::string("PHIC"));
    int pointGroup = FileParser::getKey("POINT_GROUP", 23);
    
    std::ostringstream fsmelt_script;
    fsmelt_script << "inco 1 2 " << fObsLab << " " << sigFLab << std::endl;
    fsmelt_script << "inco 2 2 " << fCalcLab << " " << phiLab << std::endl;
    fsmelt_script << "scal 1" << std::endl;
    fsmelt_script << "poin " << pointGroup << std::endl;
    fsmelt_script << "ADD" << std::endl;
    fsmelt_script << "go" << std::endl;
    
    std::string invocation = "fsmelt HKLIN1 " + originalMtz + " HKLIN2 " + tempMtz + " HKLOUT " + fsmeltMtz;
    
    runProgram(invocation, fsmelt_script.str(), "Mtz merge (fsmelt)");
}

void Deconvoluter::cad(std::string tempSortedMtz, std::string tempMtz)
{
    double highRes = FileParser::getKey("HIGH_RESOLUTION", 0.0);
    double lowRes = FileParser::getKey("LOW_RESOLUTION", 0.0);
    
    std::ostringstream cad_script;
    cad_script << "TITL back transform unit cell" << std::endl;
    cad_script << "LABIN FILE_NUM 1 ALL" << std::endl;
    cad_script << "RESOLUTION OVERALL " << lowRes << " " << highRes << std::endl;
    cad_script << "SYMMETRY " << spaceGroup << std::endl;
    cad_script << "end" << std::endl;
    
    std::string invocation = "cad HKLIN1 " + tempSortedMtz + " HKLOUT " + tempMtz;
    runProgram(invocation, cad_script.str(), "Cad");
}

void Deconvoluter::sort_mtz(std::string tempMtz, std::string tempSortedMtz)
{
    std::ostringstream sort_mtz_script;
    sort_mtz_script << "H K L" << std::endl << "END" << std::endl;
    
    std::string invocation = "sortmtz HKLIN " + tempMtz + " HKLOUT " + tempSortedMtz;
    runProgram(invocation, sort_mtz_script.str(), "Sort MTZ");
}

void Deconvoluter::map_to_sf(std::string averageMap, std::string tempMtz)
{
    std::ostringstream map_to_sfScript;
    double highRes = FileParser::getKey("HIGH_RESOLUTION", 0.0);
    double lowRes = FileParser::getKey("LOW_RESOLUTION", 0.0);
    
    map_to_sfScript << "TITL back transform unit cell" << std::endl;
    map_to_sfScript << "RESO " << lowRes << " " << highRes << std::endl;
    map_to_sfScript << "end" << std::endl;
    
    std::string invocation = "map_to_sf MAPIN " + averageMap + " HKLOUT " + tempMtz;
    
    runProgram(invocation, map_to_sfScript.str(), "Map to Structure Factor");
}

std::string Deconvoluter::scriptCentreOnVirion()
{
    std::vector<double> origin = FileParser::getKey("ORIGIN_IN_UNIT_CELL", std::vector<double>(3, 0.0));
    int gridOrigins[3];
    int minOrigins[3]; int maxOrigins[3];
    
    std::cout << "Using origin of (" << origin[0] << ", " << origin[1] << ", " << origin[2] << ") as fractions of unit cell length." << std::endl;
    
    for (int i = 0; i < 3; i++)
    {
        gridOrigins[i] = origin[i] * gridding[i];
        minOrigins[i] = gridOrigins[i] - gridding[i] / 2;
        maxOrigins[i] = gridOrigins[i] + gridding[i] / 2;
    }

    std::ostringstream centreScript;
    
    centreScript << "upda" << std::endl;
    for (int i = 0; i < 3; i++)
    {
        centreScript << minOrigins[i] << " " << maxOrigins[i] << " " << gridding[i] << std::endl;
    }
    centreScript << "GO" << std::endl;
    
    return centreScript.str();
}

void Deconvoluter::gapAverage(std::string inputMap, std::string inputWang, std::string inputMolUnit)
{
    std::ostringstream gapScript;
    std::string ncsDefs = FileParser::getKey("NCS_DEFINITIONS_FILE", std::string(""));
    
    gapScript << "TITLE \"Cyclically average map\"" << std::endl;
    gapScript << "STAT UNKNOWN\n! First of all read in map...\n\nASSI MAP1 MAPI\nRMAP " << inputMap << "\n{ read in cryst a.u. of map }" << std::endl;
    gapScript << "{ read in wang env }\nass env1 mapi\nrzip " << inputWang << std::endl;
    gapScript << "{ read in mol unit env }\nass env2 mapi\nrzip " << inputMolUnit << std::endl;
    
    gapScript << "{ read in operators }\n\n@" << ncsDefs << std::endl;
    gapScript << "{ now av map1 into map2 expanded to cover virion }\n\nycsy\n\nass map1 mapi\nass off  envi\nass env2 envo\nass map3 mapo\naver" << std::endl;
    
    gapScript << scriptCentreOnVirion() << "\ngo\ngo\n" << std::endl;
    gapScript << "{ now flatten asym unit using wang env }\nass  map1 mapi\nass  env1 envi\nass off mapo\nass off envo \nflat \ngo\n\n{ then fold averaged density down into asym unit }\n{ using expaned env as a filter }\n\nass map3 mapi\nass map1 mapo\nass env2 envi\nass off envo\nfold\n\n{ delete map2, env2...for space }\ndmap map2\ndmap env2\ndmap map3\n\n{ then need to expand to fill cell for SFALL }\n\nass map1 mapi\nass map4 mapo\ncmap\n" << std::endl;
    
    gapScript << "upda" << std::endl;
    for (int i = 0; i < 3; i++)
    {
        gapScript << "0  " << gridding[i] << " " << gridding[i] << std::endl;
    }
    gapScript << "rese 2 1 3 go" << std::endl;
    gapScript << "ass map4 mapi\nsgrp map4 1\nlmap inte 6 6 6 \n\ngo\n\nass map4 mapi\nwbin " << tempAveMapName << "\n\n\nstop\n" << std::endl;

    std::string invocation = "gap";
    
    runProgram(invocation, gapScript.str(), "Gap Density Averaging");
}

void Deconvoluter::gapEnvelope(std::string inputMap)
{
     double generousDiameter = FileParser::getKey("VIRION_GENEROUS_DIAMETER", 340.0);
    double outerDiameter = FileParser::getKey("VIRION_OUTER_DIAMETER", 300.0);
    double innerDiameter = FileParser::getKey("VIRION_INNER_DIAMETER", 220.0);
    double sqrtDiameter = generousDiameter / 1.41421;
    
    std::cout << "Using generous diameter of " << generousDiameter << " Å, inner diameter of " << innerDiameter << ", outer diameter of " << outerDiameter << " Å." << std::endl;
    
    
    
    std::ostringstream envScript;
    
    envScript << "#noecho\nSTAT UNKNOWN\nTITLE \"Envelope determination\"{ read in asymmetric unit }\nass map1 mapi\nrmap " << inputMap << std::endl;
    
    envScript << "assi map1 mapi\nassi map2 mapo\ncmap " << std::endl;
    
    envScript << scriptCentreOnVirion() << std::endl;
    
    envScript << "\n{now use planes to cut away the other virions from our chosen one}" << std::endl;
    envScript << "\nncsy\nass map2 mapi\nass env1 envm\nncsy\nMENV GO\n0" << std::endl;
    
    envScript << "plane " << generousDiameter << " 0.0 0.0 0.0 0.0 0.0 1" << std::endl;
    envScript << "plane 0.0 " << generousDiameter << " 0.0 0.0 0.0 0.0 1" << std::endl;
    envScript << "plane 0.0 0.0 " << generousDiameter << " 0.0 0.0 0.0 1" << std::endl;
    envScript << "plane " << -generousDiameter << " 0.0 0.0 0.0 0.0 0.0 1" << std::endl;
    envScript << "plane 0.0 " << -generousDiameter << " 0.0 0.0 0.0 0.0 1" << std::endl;
    envScript << "plane 0.0 0.0 " << -generousDiameter << " 0.0 0.0 0.0 1" << std::endl;
    
    envScript << "plane " << sqrtDiameter << " " << sqrtDiameter << " " << sqrtDiameter << " 0.0 0.0 0.0 1" << std::endl;
    envScript << "plane " << -sqrtDiameter << " " << -sqrtDiameter << " " << -sqrtDiameter << " 0.0 0.0 0.0 1" << std::endl;
    envScript << "plane " << sqrtDiameter << " " << sqrtDiameter << " " << -sqrtDiameter << " 0.0 0.0 0.0 1" << std::endl;
    envScript << "plane " << -sqrtDiameter << " " << -sqrtDiameter << " " << sqrtDiameter << " 0.0 0.0 0.0 1" << std::endl;
    envScript << "plane " << -sqrtDiameter << " " << sqrtDiameter << " " << sqrtDiameter << " 0.0 0.0 0.0 1" << std::endl;
    envScript << "plane " << sqrtDiameter << " " << -sqrtDiameter << " " << -sqrtDiameter << " 0.0 0.0 0.0 1" << std::endl;
    envScript << "plane " << sqrtDiameter << " " << -sqrtDiameter << " " << sqrtDiameter << " 0.0 0.0 0.0 1" << std::endl;
    envScript << "plane " << -sqrtDiameter << " " << sqrtDiameter << " " << -sqrtDiameter << " 0.0 0.0 0.0 1" << std::endl;
    envScript << "go\n\n{now make mask that defines our virion as an inner sphere and outer sphere}" << std::endl;
    
    envScript << "assi env1 mapi\nassi env1 envm\nncsy\nmenv usei go" << std::endl;

    envScript << "sphere 0.0 0.0 0.0 " << outerDiameter / 2 << " 0 1" << std::endl;
    envScript << "sphere 0.0 0.0 0.0 " << innerDiameter / 2 << " 0 0\ngo" << std::endl;
    
    envScript << "ass env1 mapi\nlmas inte 10 10 10 go\n\n{now set density outside env1 in map2 to 0}\n\nncsy\n\nass map2 mapi\nass env1 envi\nflat\nflot 0\ngo\nass map2 mapi\nlmas inte 10 10 10 go\n{now make map2 the W.E.}\nncsy\n{Sig cut it}\nass off envi\nass map2 mapi\nscut .2 5 go\n\nass map2 mapi\nlmas inte 10 10 10 go\n\nncsy\n\n{convolute it}\nass map2 mapi\nconv bfac 200 go\nncsy\n\nass map2 mapi\nlmas inte 10 10 10 go\n\n{histogram it}\nass map2 mapi\nass env2 envo\nhist .75\n\nass env2 mapi\nlmas inte 10 10 10 go\n\nncsy\nass env2 mapi\ntidy go\nass env2 mapi\ntidy go\nass env2 mapi\nlmas inte 10 10 10 go\n\n{now write out mol unit for 1/4 virion}\n\nass env2 mapi\nwzip " << molUnitName << "\n\n{now create wang env in asymmetric unit}\n\nycsy\nass map5 mapo\n\nass env2 mapi\nass off envo\nass env2 envi\n\nfold\nupda\n0 180 360\n0 180 360\n0 180 360\ngo\n\n{lms the final env}\nass map5 mapi\nlmas inte 10 10 10 go\n\nass map5 mapi\nwzip " << wangEnvName << "\n\nstop" << std::endl;
    
    std::string invocation = "gap";
    
    runProgram(invocation, envScript.str(), "Gap Envelope Generation");
}


void Deconvoluter::fastFourierTransform(std::string inputName, std::string tempFftMapName, bool calcAmps)
{
    std::ostringstream fftScript;
    double highRes = FileParser::getKey("HIGH_RESOLUTION", 0.0);
    double lowRes = FileParser::getKey("LOW_RESOLUTION", 0.0);
    
    std::string fObsLab = FileParser::getKey("LABIN_AMPLITUDES", std::string("F"));
    std::string fCalcLab = FileParser::getKey("LABOUT_CALC_AMPLITUDES", std::string("FC"));
    std::string sig1Lab = FileParser::getKey("LABIN_SIGMAS", std::string("SIGF"));
    std::string phiLab = FileParser::getKey("LABIN_PHASES", std::string("PHIC"));
    std::string f1Lab = calcAmps ? fCalcLab : fObsLab;
    
    fftScript << "TITLE new 2fo-fc map" << std::endl;
    
    if (lowRes > 0 && highRes > 0)
    {
        fftScript << "RESOLUTION " << lowRes << " " << highRes << std::endl;
    }
    
    fftScript << "SCAL F1 2.0 0.0" << std::endl << "SCAL F2 1.0 0.0" << std::endl;
    fftScript << "GRID " << gridding[0] << " " << gridding[1] << " " << gridding[2] << std::endl;
    fftScript << "SYMM " << spaceGroup << std::endl;
    fftScript << "FFTSYMMETRY 1" << std::endl;
    fftScript << "XYZLIM 0 " << gridding[0] << " 0 " << gridding[1] << " 0 " << gridding[2] << std::endl;
    fftScript << "LABIN F1=" << f1Lab << " SIG1=" << sig1Lab << " PHI=" << phiLab << std::endl;
    
    std::string invocation = "fft HKLIN " + inputName + " " + " MAPOUT " + tempFftMapName;
    
    runProgram(invocation, fftScript.str(), "Fast Fourier Transform");
}

void Deconvoluter::process()
{
    std::string file = FileParser::getKey("TWINNED_MTZ", std::string(""));
    wangEnvName = FileParser::getKey("WANG_ENVELOPE_OUTPUT", std::string("wang.env"));
    molUnitName = FileParser::getKey("MOL_UNIT_ENVELOPE_OUTPUT", std::string("mol_unit.env"));
    tempFftMapName = FileParser::getKey("TEMP_FFT_MAP", std::string("temp_fft.map"));
    tempAveMapName = FileParser::getKey("TEMP_AVERAGE_MAP", std::string("temp_ave.map"));
    tempMtzName = FileParser::getKey("TEMP_AVERAGE_MTZ", std::string("temp_ave.mtz"));
    tempSortedMtzName = FileParser::getKey("TEMP_AVERAGE_SORTED_MTZ", std::string("temp_sorted_ave.mtz"));
    cadMtzName = FileParser::getKey("TEMP_CAD_MTZ", std::string("temp_cad.mtz"));
    fsmeltMtzName = FileParser::getKey("TEMP_FSMELT_MTZ", std::string("temp_fsmelt.mtz"));
    std::string fObsLab = FileParser::getKey("LABIN_AMPLITUDES", std::string("F"));
    std::string fCalcLab = FileParser::getKey("LABOUT_CALC_AMPLITUDES", std::string("FC"));
    std::string ncsDefs = FileParser::getKey("NCS_DEFINITIONS_FILE", std::string(""));

    std::vector<double> ccAll, ccSinglets, rAll, rSinglets;
    
    if (ncsDefs.length() == 0)
    {
        std::cout << "Warning! NCS definitions file missing. Please specify file in GAP format using keyword NCS_DEFINITIONS_FILE" << std::endl;
        exit(1);
    }
    
    int bins = 20;
    int maxCycles = FileParser::getKey("MAXIMUM_CYCLES", 10);
    
    if (file == "")
    {
        std::cout << "Twinned MTZ has not been provided, please provide file path under keyword TWINNED_MTZ" << std::endl;
        exit(1);
    }
    
    // loading original, twinned MTZ into memory.
    originalMtz = MtzPtr(new MtzManager(file));
    originalMtz->loadReflections(false);
    MtzManager::setReference(&*originalMtz);
    spaceGroup = FileParser::getKey("SPACE_GROUP", originalMtz->getLowGroup()->spg_num);
    
    
    
    std::cout << "Loaded original MTZ file " << file << std::endl;
    bool skipFirst = FileParser::getKey("SKIP_FIRST_CYCLE", false);
    
    int beginning = FileParser::getKey("RESUME_FROM_CYCLE", 0);
    if (beginning > 0) skipFirst = true;
    
    for (int i = beginning; i < maxCycles; i++)
    {
        if (i > beginning || (i == beginning && !skipFirst))
        {
            fastFourierTransform(file, tempFftMapName, (i > 0));
            gapEnvelope(tempFftMapName);
            gapAverage(tempFftMapName, wangEnvName, molUnitName);
            map_to_sf(tempAveMapName, tempMtzName);
            sort_mtz(tempMtzName, tempSortedMtzName);
            cad(tempSortedMtzName, cadMtzName);
            fsmelt(file, cadMtzName, fsmeltMtzName);
        }
        
        MtzPtr nextMtz = MtzPtr(new MtzManager(fsmeltMtzName));
        nextMtz->loadReflections(true); // fc in "intensity" and f in "fc"...
        
        nextMtz->applyScaleFactorsForBins(bins);
        
        nextMtz->individualDetwinningScales((i == maxCycles - 1));
        nextMtz->copyOtherAmplitudesFromReference();
        file = "detwinned_cycle_" + i_to_str(i) + ".mtz";
        nextMtz->writeToFile(file);
        
        std::cout << "******************************" << std::endl;
        std::cout << "**    CORRELATION (ALL)     **" << std::endl;
        std::cout << "******************************" << std::endl;
        std::cout << std::endl << "Correlation between all scaled data and original twinned data" << std::endl;
        std::cout << std::endl << std::setw(15) << "Low res " << std::setw(15) << "High res " << std::setw(15) << "Correl" << std::setw(15) << "Num refl" << std::endl;
        
        ccAll.push_back(nextMtz->correlationWithManager(&*originalMtz, false, false, 0, 0, bins, NULL, false));
        
        std::cout << "*******************************" << std::endl;
        std::cout << "**  CORRELATION (SINGLETS)   **" << std::endl;
        std::cout << "*******************************" << std::endl;
        std::cout << std::endl << "Correlation between singlet scaled data and original twinned data" << std::endl;
        std::cout << std::endl << std::setw(15) << "Low res " << std::setw(15) << "High res " << std::setw(15) << "Correl" << std::setw(15) << "Num refl" << std::endl;
        
        ccSinglets.push_back(nextMtz->correlationWithManager(&*originalMtz, false, false, 0, 0, bins, NULL, true));

        std::cout << "******************************" << std::endl;
        std::cout << "**      R FACTOR (ALL)      **" << std::endl;
        std::cout << "******************************" << std::endl;
        std::cout << std::endl << "R factor between all scaled data and original twinned data" << std::endl;
        std::cout << std::endl << std::setw(15) << "Low res " << std::setw(15) << "High res " << std::setw(15) << "Correl" << std::setw(15) << "Num refl" << std::endl;
        
        rAll.push_back(nextMtz->rSplitWithManager(&*originalMtz, false, false, 0, 0, bins, NULL, false));

        std::cout << "*******************************" << std::endl;
        std::cout << "**   R FACTOR (SINGLETS)     **" << std::endl;
        std::cout << "*******************************" << std::endl;
        std::cout << std::endl << "R factor between singlet scaled data and original twinned data" << std::endl;
        std::cout << std::endl << std::setw(15) << "Low res " << std::setw(15) << "High res " << std::setw(15) << "Correl" << std::setw(15) << "Num refl" << std::endl;
        
        rSinglets.push_back(nextMtz->rSplitWithManager(&*originalMtz, false, false, 0, 0, bins, NULL, true));
    }
    
    std::cout << "*******************************" << std::endl;
    std::cout << "**   END OF DECONVOLUTION    **" << std::endl;
    std::cout << "*******************************" << std::endl;
    std::cout << std::endl << "Summary of deconvolution:" << std::endl << std::endl;
    std::cout << "Cycle\tCCall\tCCsinglets\tRall\tRsinglets" << std::endl;
    
    for (int i = 0; i < ccAll.size(); i++)
    {
        std::cout << i << "\t" << ccAll[i] << "\t" << ccSinglets[i] << "\t" << rAll[i] << "\t" << rSinglets[i] << std::endl;
    }
}