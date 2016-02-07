/*
 * Miller.h
 *
 *  Created on: 27 Aug 2014
 *      Author: helenginn
 */

#ifndef MILLER_H_
#define MILLER_H_

#include "Matrix.h"
#include <string>

#include <map>
#include "MtzManager.h"
#include "definitions.h"
#include "parameters.h"
#include "Holder.h"
#include <boost/weak_ptr.hpp>

class Image;

typedef enum
{
	CalculationTypeOriginal, CalculationTypeIntegrate
} CalculationType;

class Miller
{
private:
    double h;
    double k;
    double l;
    short int fakeFriedel;
	bool normalised;
    bool correctingPolarisation;
	double polarisationCorrection;
    double polarisationFactor;
	std::map<std::string, bool> rejectedReasons;
	double partialCutoff;
	double bFactor;
	double scale;
    double gainScale;
	double lastX;
	double lastY;
    double lastBandwidth;
    double lastWavelength;
    double lastRlpSize;
    double lastMosaicity;
    int slices;
    double trickyRes;
    double phase;
    double fc;
    double maxSlices;
    
    double lastRadius;
    double lastVolume;
    double lastSurfaceArea;
    bool sizeChanged;
    
    double latestHRot;
    double latestKRot;
    double bFactorScale;
    bool excluded;
    bool rejected;
    bool calculatedRejected;
    double additionalWeight;
    double denormaliseFactor;
    
	std::pair<double, double> shift;
    double resol;
 
    boost::weak_ptr<Miller> selfPtr;
    MatrixPtr flipMatrix;
public:
    int getH();
    int getK();
    int getL();
    
    vec hklVector(bool shouldFlip = true);
    void setFlipMatrix(MatrixPtr flipMat);
	MatrixPtr matrix;
	Reflection *parentReflection;
    MtzManager *mtzParent;

	Miller(MtzManager *parent, int _h = 0, int _k = 0, int _l = 0);
	MillerPtr copy(void);
	void printHkl(void);
	static double scaleForScaleAndBFactor(double scaleFactor, double bFactor, double resol, double exponent_exponent = 1);

    void setData(double _intensity, double _sigma, double _partiality,
			double _wavelength);
	void setParent(Reflection *reflection);
	void setFree(bool newFree);
	bool activeWithAngle(double degrees);
	void setRejected(std::string reason, bool rejection);
	bool isRejected(std::string reason);
	void makeScalesPermanent();
    void integrateIntensity(double hRot, double kRot);
    vec getTransformedHKL(double hRot, double kRot);
    double getWeight();
    
	bool accepted(void);
	bool free(void);
	void flip(void);

    bool isRejected();
    double getBFactorScale();
	double intensity(void);
	double getSigma(void);
	double resolution();

	void applyScaleFactor(double scaleFactor);
	
    void makeComplexShoebox(double wavelength, double bandwidth, double mosaicity, double rlpSize);
    
    static double averageRawIntensity(vector<MillerPtr> millers);

	virtual ~Miller();
    
    double getPhase()
    {
        return phase;
    }
    
    double getFc()
    {
        return fc;
    }
    
    void setPhase(double thePhase)
    {
        phase = thePhase;
    }
    
    void setFc(double fcalc)
    {
        this->fc = fcalc;
    }
    
    void setExcluded(bool exc = true)
    {
        excluded = exc;
    }
    
    bool isExcluded()
    {
        return excluded;
    }
    
    void setAdditionalWeight(double weight)
    {
        additionalWeight = weight;
    }
    
    double getAdditionalWeight()
    {
        return additionalWeight;
    }
    
    MtzManager *&getMtzParent()
    {
        return mtzParent;
    }

    void setMtzParent(MtzManager *mtz)
    {
        mtzParent = mtz;
    }
    
	void setPartiality(double partiality)
	{
		this->partiality = partiality;
	}
    
    double getRawestIntensity() const
    {
        return rawIntensity;
    }

	double getRawIntensity() const
	{
		return rawIntensity * scale * gainScale;
	}

    void setScaledIntensity(double scaledIntensity)
    {
        this->rawIntensity = scaledIntensity;
        scale = 1;
    }
    
	void setRawIntensity(double rawIntensity)
	{
		this->rawIntensity = rawIntensity;
	}

	void setSigma(double sigma)
	{
		this->sigma = sigma;
	}
    
    void setGainScale(double newScale)
    {
        if (newScale == newScale)
            gainScale = newScale;
    }

	const std::string& getFilename() const
	{
		return filename;
	}

	void setFilename(const std::string& filename)
	{
		this->filename = filename;
	}

	double getCountingSigma() const
	{
        return countingSigma;
	}

	void setCountingSigma(double countingSigma)
	{
		this->countingSigma = countingSigma;
	}

	bool isNormalised() const
	{
		return normalised;
	}

	void setNormalised(bool normalised)
	{
		this->normalised = normalised;
	}

	MatrixPtr getMatrix()
	{
		return matrix;
	}

	void setMatrix(MatrixPtr matrix)
	{
		this->matrix = matrix;
	}

	void setPolarisationCorrection(double polarisationCorrection)
	{
		this->polarisationCorrection = polarisationCorrection;
	}

	double getLastWavelength() const
	{
		return lastWavelength;
	}

	double getLastBandwidth()
    {
        return lastBandwidth;
    }
    
    double getLastMosaicity()
    {
        return lastMosaicity;
    }
    
    double getLastRlpSize()
    {
        return lastRlpSize;
    }

	void setRejected(bool rejected)
	{
        setRejected("merge", rejected);        
	}

	double getLastX() const
	{
		return lastX;
	}

	double getLastY() const
	{
		return lastY;
	}
    
    std::pair<double, double> position()
    {
        return std::make_pair(lastX, lastY);
    }

	void setLastX(double lastX)
	{
		this->lastX = lastX;
	}

	void setLastY(double lastY)
	{
		this->lastY = lastY;
	}

	double getPartialCutoff() const
	{
		return partialCutoff;
	}

	void setPartialCutoff(double partialCutoff)
	{
		this->partialCutoff = partialCutoff;
	}

	double getBFactor() const
	{
		return bFactor;
	}

	void setBFactor(double factor)
	{
		if (factor == factor)
			bFactor = factor;
        
        bFactorScale = 0;
	}

	double getScale() const
	{
		return scale;
	}

	void setScale(double scale)
	{
		if (scale == scale)
			this->scale = scale;
        
    /*    if ((h == 3 && k == -4 && l == 12) || (h == -4 && k == 3 && l == -12))
        {
            std::cout << "Changing (" << h << ", " << k << ", " << l << ") to " << scale << std::endl;
        }*/
	}

	double getResolution()
	{
        return resolution();
	}

	void setResolution(double resol)
	{
		this->resol = resol;
	}

	std::pair<double, double>& getShift()
	{
		return shift;
	}

	void setShift(const std::pair<int, int>& shift)
	{
		this->shift = shift;
	}
        
    void setSelf(MillerPtr ptr)
    {
        selfPtr = ptr;
    }

protected:
	double rawIntensity;
	double sigma;
	double countingSigma;
	double partiality;
	double wavelength;
	std::string filename;

	void rotateMatrix(double hRot, double kRot, MatrixPtr *newMatrix);
};

#endif /* MILLER_H_ */
