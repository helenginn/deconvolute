//
//  main.cpp
//  deconvolute
//
//  Created by Helen Ginn on 18/07/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include <iostream>
#include "InputFileParser.h"
#include "MtzManager.h"

int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "********************" << std::endl;
    std::cout << "* DECONVOLUTE V1.0 *" << std::endl;
    std::cout << "********************" << std::endl << std::endl;
    
    if (argc < 3)
    {
        std::cout << "arguments: -i <input_script>" << std::endl;
        exit(1);
    }
    
    if (strcmp(argv[1], "-i") == 0)
    {
        if (argc < 3)
        {
            std::cout << "arguments: -i <input_script>" << std::endl;
            exit(1);
        }
        
        std::cout << "Loading input script " << std::string(argv[2]) << std::endl;
        
        InputFileParser *parser = new InputFileParser(std::string(argv[2]));
        parser->parse(false);
        
        delete parser;
    }
    
    if (strcmp(argv[1], "-cc") == 0 || strcmp(argv[1], "-r") == 0 || strcmp(argv[1], "-ri") == 0)
    {
        if (argc < 3)
        {
            std::cout
            << "arguments: -cc / -r / -ri <file1> {[singlets (0/1)] {[lowRes] [highRes] {[bins]}}} [bfac]."
            << std::endl;
            exit(1);
        }
        
        bool rsplit = (strcmp(argv[1], "-r") == 0);
        bool rsplitIntensity = (strcmp(argv[1], "-ri") == 0);
        
        double highRes = 0;
        double lowRes = 0;
        int bins = 20;
        bool singlets = false;
        double bFac = 0;
        
        if (argc >= 4)
        {
            singlets = atoi(argv[3]);
        }
        if (argc >= 6)
        {
            lowRes = atof(argv[4]);
            highRes = atof(argv[5]);
        }
        if (argc >= 7)
        {
            bins = atoi(argv[6]);
        }
        if (argc >= 8)
        {
            bFac = atof(argv[7]);
        }
        
        std::cout << "Singlets: " << singlets << std::endl;
        
        MtzManager *mtz1 = new MtzManager(std::string(argv[2]));
        mtz1->loadReflections(false);
        mtz1->applyBFactor(bFac);
        
        MtzManager *mtz2 = new MtzManager(std::string(argv[2]));
        mtz2->loadReflections(true);
        
        if (bFac != 0)
        {
            double scale = mtz1->gradientAgainstManager(*mtz2);
            mtz1->applyScaleFactor(scale);
        }
        
        if (rsplit)
        {
            mtz1->rSplitWithManager(mtz2, 1, 0, lowRes, highRes, bins, NULL, singlets);
        }
        else if (rsplitIntensity)
        {
            mtz1->rSplitIntensityBinsWithManager(mtz2, 1, 0, lowRes, highRes, bins, NULL, singlets);
        }
        else
        {
            mtz1->correlationWithManager(mtz2, 1, 0, lowRes, highRes, bins, NULL, singlets);
            
        }
        
        if (bFac != 0)
        {
            mtz1->writeToFile("bFac-" + mtz1->getFilename());
        }
        
        delete mtz1;
        delete mtz2;
        
        exit(1);
    }

    
    return 0;
}
