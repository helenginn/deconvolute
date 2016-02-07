/*
 * Miller.cpp
 *
 *  Created on: 27 Aug 2014
 *      Author: helenginn
 */

#include "parameters.h"
#include "Miller.h"

#include "definitions.h"
#include "Vector.h"
#include <cmath>
#include <cmath>
#include <memory>
#include "FileParser.h"

Miller::Miller(MtzManager *parent, int _h, int _k, int _l)
{
    h = _h;
    k = _k;
    l = _l;
    lastX = 0;
    lastY = 0;
    rawIntensity = 0;
    gainScale = 1;
    sigma = 1;
    wavelength = 0;
    partiality = -1;
    filename = "";
    countingSigma = 0;
    latestHRot = 0;
    latestKRot = 0;
    normalised = true;
    correctingPolarisation = FileParser::getKey("POLARISATION_CORRECTION", false);
    polarisationFactor = FileParser::getKey("POLARISATION_FACTOR", 0.0);
    polarisationCorrection = 0;
    lastWavelength = 0;
    rejectedReasons["merge"] = false;
    scale = 1;
    bFactor = 0;
    bFactorScale = 0;
    resol = 0;
    shift = std::make_pair(0, 0);
    selfPtr = MillerPtr();
    fakeFriedel = -1;
    rejected = false;
    calculatedRejected = true;
    denormaliseFactor = 1;
    excluded = false;
    lastRlpSize = 0;
    lastMosaicity = 0;
    lastVolume = 0;
    lastSurfaceArea = 0;
    slices = FileParser::getKey("PARTIALITY_SLICES", 8);
    trickyRes = FileParser::getKey("CAREFUL_RESOLUTION", 8);
    maxSlices = FileParser::getKey("MAX_SLICES", 100);
    
    partialCutoff = FileParser::getKey("PARTIALITY_CUTOFF",
                                       PARTIAL_CUTOFF);
    
    mtzParent = parent;
    parentReflection = NULL;
    matrix = MatrixPtr();
    flipMatrix = MatrixPtr(new Matrix());
}

vec Miller::hklVector(bool shouldFlip)
{
    vec newVec = new_vector(h, k, l);
    
    if (shouldFlip)
    {
        flipMatrix->multiplyVector(&newVec);
    }
    
    return newVec;
}

int Miller::getH()
{
    return hklVector().h;
}

int Miller::getK()
{
    return hklVector().k;

}

int Miller::getL()
{
    return hklVector().l;
}

void Miller::setFlipMatrix(MatrixPtr flipMat)
{
    flipMatrix = flipMat;
}

void Miller::setParent(Reflection *reflection)
{
    parentReflection = reflection;
}


void Miller::setData(double _intensity, double _sigma, double _partiality,
                     double _wavelength)
{
    rawIntensity = _intensity;
    sigma = _sigma;
    partiality = _partiality;
    wavelength = _wavelength;
}

void Miller::printHkl(void)
{
    std::cout << "h k l " << h << " " << k << " " << l << std::endl;
}

bool Miller::accepted(void)
{
        return true;
}

bool Miller::isRejected()
{
    if (calculatedRejected)
        return rejected;
    
    rejected = false;
    
    for (std::map<std::string, bool>::iterator it = rejectedReasons.begin();
         it != rejectedReasons.end(); ++it)
    {
        if (rejectedReasons[it->first])
            rejected = true;
    }
    
    calculatedRejected = true;
    
    return rejected;
}

double Miller::resolution()
{
    if (resol == 0)
    {
        vec newVec = getTransformedHKL(0, 0);
        
        resol = length_of_vector(newVec);
    }
    
    return resol;
}

double Miller::getBFactorScale()
{
    if (bFactorScale != 0)
        return bFactorScale;
    
    double resn = getResolution();
    
    double sinThetaOverLambda = 1 / (2 / resn);
    
    double factor = exp(- bFactor * pow(sinThetaOverLambda, 2));
    
    bFactorScale = factor;
    
    return factor;
}

double Miller::scaleForScaleAndBFactor(double scaleFactor, double bFactor,
                                      double resn, double exponent_exponent)
{
    double four_d_squared = 4 * pow(resn, 2);
    
    double right_exp = pow(1 / four_d_squared, exponent_exponent);
    
    double four_d_to_exp = pow(2, exponent_exponent) * right_exp;
    
    double exponent = pow(bFactor, exponent_exponent) * four_d_to_exp;
    
    double scale = scaleFactor * exp(pow(exponent, exponent_exponent));
    
    return scale;
}

double Miller::intensity()
{
    
    double modifier = scale * gainScale;
    modifier /= getBFactorScale();
    
    return modifier * rawIntensity;
    
}

void Miller::makeScalesPermanent()
{
    double scale = scaleForScaleAndBFactor(this->scale, this->bFactor,
                                          1 / this->resol);
    
    this->rawIntensity *= scale;
    this->sigma *= scale;
    
    this->scale = 1;
    bFactor = 0;
}

double Miller::getSigma(void)
{
    // bigger sigma - larger error
    
    double correction = gainScale * scale;
    /*
    double bFactorScale = getBFactorScale();
    correction /= bFactorScale;
    
    if (bFactorScale == 0) correction = 0;*/
    
    return sigma * correction;
}

vec Miller::getTransformedHKL(double hRot, double kRot)
{
    MatrixPtr newMatrix = MatrixPtr();
    rotateMatrix(hRot, kRot, &newMatrix);
    
    vec hkl = new_vector(h, k, l);
    
    newMatrix->multiplyVector(&hkl);

    return hkl;
}

double Miller::getWeight()
{
    if (!this->accepted())
    {
        return 0;
    }
    
    double weight = 1;
    double sigma = this->getSigma();
    
    weight = 1 / (sigma * sigma);
    return weight;
}

void Miller::rotateMatrix(double hRot, double kRot, MatrixPtr *newMatrix)
{
    (*newMatrix) = matrix->copy();
    
    double hRad = hRot * M_PI / 180;
    double kRad = kRot * M_PI / 180;
    
    (*newMatrix)->rotate(hRad, kRad, 0);
}

MillerPtr Miller::copy(void)
{
    MillerPtr newMiller = MillerPtr(new Miller(mtzParent));
    
    newMiller->h = h;
    newMiller->k = k;
    newMiller->l = l;
    newMiller->rawIntensity = rawIntensity;
    newMiller->sigma = sigma;
    newMiller->matrix = matrix;
    newMiller->filename = std::string(filename);
    newMiller->countingSigma = countingSigma;
    newMiller->lastX = lastX;
    newMiller->lastY = lastY;
    newMiller->shift = shift;
    
    return newMiller;
}

void Miller::flip(void)
{
    int tmp = l;
    l = k;
    k = tmp;
}

void Miller::applyScaleFactor(double scaleFactor)
{
    setScale(scale * scaleFactor);
}

void flattenAngle(double *radians)
{
    while (*radians < 0)
    {
        *radians += M_PI;
    }
    
    while (*radians > 2 * M_PI)
    {
        *radians -= M_PI;
    }
}

void Miller::setRejected(std::string reason, bool rejection)
{
    rejectedReasons[reason] = rejection;
    
    if (!rejection)
    {
        rejectedReasons.erase(reason);
    }
    
    calculatedRejected = false;
}

bool Miller::isRejected(std::string reason)
{
    if (rejectedReasons.count(reason) == 0)
        return false;
    
    return rejectedReasons[reason];
}

double Miller::averageRawIntensity(vector<MillerPtr> millers)
{
    double allIntensities = 0;
    int num = 0;
    
    for (int i = 0; i < millers.size(); i++)
    {
        MillerPtr miller = millers[i];
        
        allIntensities += miller->getRawIntensity();
        num++;
    }
    
    allIntensities /= num;
    
    return allIntensities;
}

Miller::~Miller()
{
    
}

// Vector stuff

