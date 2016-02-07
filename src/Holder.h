/*
 * Holder.h
 *
 *  Created on: 25 Oct 2014
 *      Author: helenginn
 */

#ifndef HOLDER_H_
#define HOLDER_H_

#include "csymlib.h"

class Reflection;

#include "parameters.h"
#include "Miller.h"
#include "csymlib.h"

class Miller;

class Reflection
{
private:
	vector<MillerPtr> millers;
	int refl_id;
	int inv_refl_id;
    double unitCell[6];
	double refIntensity;
	double refSigma;
	double resolution;
    bool flagScaled;
    CCP4SPG *spaceGroup;
    int spgNum;
    static bool hasSetup;
    static bool setupUnitCell;
    int activeAmbiguity;
    vector<int> reflectionIds;
public:
    Reflection(float *unitCell = NULL, CSym::CCP4SPG *group = NULL);
    void setUnitCell(float *unitCell);
	virtual ~Reflection();

    MillerPtr acceptedMiller(int num);
	MillerPtr miller(int i);
    void printDescription();
	void addMiller(MillerPtr miller);
	int millerCount();
	Reflection *copy(bool copyMillers);
	Reflection *copy();
    
	static int indexForReflection(int h, int k, int l, CSym::CCP4SPG *lowspgroup, bool inverted);

	void reflectionDescription();
	void calculateResolution(MtzManager *mtz);
	void clearMillers();
    void removeMiller(int index);
    bool isTwinned();

	double meanIntensityWithExclusion(std::string *filename, int start = 0, int end = 0);
	double meanIntensity(int start = 0, int end = 0);
	double meanSigma();
	double meanSigma(bool friedel);

	bool betweenResolutions(double lowAngstroms, double highAngstroms);
	bool anyAccepted();
	int acceptedCount();
	int rejectCount();
    double meanWeight();
	
    int reflectionIdForMiller(vec miller);
    void generateReflectionIds();
    void setAdditionalWeight(double weight);
    
    void incrementAmbiguity();
    int ambiguityCount();
    MatrixPtr matrixForAmbiguity(int i);

    void anglesAsRadians(double *alpha, double *beta, double *gamma);
    double unitCellVolume();
    void reciprocalAxes(vec *aStar, vec *bStar, vec *cStar);
    void realSpaceAxes(vec *a, vec *b, vec *c);
    void setSpaceGroup(int spaceGroupNum);
    void setSpaceGroup(CSym::CCP4SPG *ccp4spg);
    
    void resetFlip();
    void setFlip(int i);
    void setFlipAsActiveAmbiguity();
    
    void setFlagScaled(bool scaled)
    {
        flagScaled = scaled;
    }
    
    bool isFlagScaled()
    {
        return flagScaled;
    }
    
    void scale(double scale);
    
    int getActiveAmbiguity()
    {
        return activeAmbiguity;
    }
    
    void setActiveAmbiguity(int newAmbiguity)
    {
        activeAmbiguity = newAmbiguity;
    }

	const vector<MillerPtr>& getMillers() const
	{
		return millers;
	}

	void setMillers(const vector<MillerPtr>& millers)
	{
		this->millers = millers;
	}

	double getRefIntensity() const
	{
		return refIntensity;
	}

	void setRefIntensity(double refIntensity)
	{
		this->refIntensity = refIntensity;
	}

    int getReflId()
    {
        if (activeAmbiguity > ambiguityCount())
        {
            std::cout << "Active ambiguity error" << std::endl;
        }
        
        return reflectionIds[activeAmbiguity];
    }
    
    int getReflId(int ambiguity)
    {
        return reflectionIds[ambiguity];
    }
    
	double getRefSigma() const
	{
		return refSigma;
	}

	void setRefSigma(double refSigma)
	{
		this->refSigma = refSigma;
	}

	double getResolution() const
	{
		return resolution;
	}

	void setResolution(double resolution)
	{
		this->resolution = resolution;
	}
};

#endif /* HOLDER_H_ */
