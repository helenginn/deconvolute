/*
 * Holder.cpp
 *
 *  Created on: 25 Oct 2014
 *      Author: helenginn
 */

#include "Holder.h"

#include <vector>
#include <cmath>
#include "csymlib.h"

#include "FileParser.h"
#include "StatisticsManager.h"

bool Reflection::hasSetup = false;
bool Reflection::setupUnitCell = false;

/**
 Hard-coding the matrix ambiguities for various twinnable space groups.
*/
MatrixPtr Reflection::matrixForAmbiguity(int i)
{
    if (i >= ambiguityCount())
        std::cout << "Ambiguity issue!" << std::endl;
    
    if (i == 0)
    {
        MatrixPtr identity = MatrixPtr(new Matrix());
        
        return identity;
    }
    
    if (i == 1)
    {
        if (ambiguityCount() == 2 || ambiguityCount() == 4)
        {
            MatrixPtr khMinusL = MatrixPtr(new Matrix());
            (*khMinusL)[0] = 0;
            (*khMinusL)[4] = 1;
            (*khMinusL)[1] = 1;
            (*khMinusL)[5] = 0;
            (*khMinusL)[10] = -1;
            
            return khMinusL;
        }
        
        if (ambiguityCount() == 3)
        {
            // return -h -k -l
            MatrixPtr minusHminusKL = MatrixPtr(new Matrix());
            (*minusHminusKL)[0] = -1;
            (*minusHminusKL)[5] = -1;
            (*minusHminusKL)[10] = 1;
            
            return minusHminusKL;
        }
    }
    
    if (i == 2)
    {
        if (ambiguityCount() == 3)
        {
            if (spgNum == 149 || spgNum == 151 || spgNum == 153)
            {
                // return k h -l
                MatrixPtr khMinusL = MatrixPtr(new Matrix());
                (*khMinusL)[0] = 0;
                (*khMinusL)[4] = 1;
                (*khMinusL)[1] = 1;
                (*khMinusL)[5] = 0;
                (*khMinusL)[10] = -1;
                
                return khMinusL;
            }
            
            if (spgNum == 152 || spgNum == 152 || spgNum == 154)
            {
                // return -k -h -l
                MatrixPtr minusAllHKL = MatrixPtr(new Matrix());
                (*minusAllHKL)[0] = -1;
                (*minusAllHKL)[5] = -1;
                (*minusAllHKL)[10] = -1;
                
                return minusAllHKL;
            }
        }

        if (ambiguityCount() == 4)
        {
            // return k h -l
            MatrixPtr khMinusL = MatrixPtr(new Matrix());
            (*khMinusL)[0] = 0;
            (*khMinusL)[4] = 1;
            (*khMinusL)[1] = 1;
            (*khMinusL)[5] = 0;
            (*khMinusL)[10] = -1;
            
            return khMinusL;
        }
    }
    
    if (i == 3)
    {
        if (ambiguityCount() == 4)
        {
            // return -k -h -l
            MatrixPtr minusAllHKL = MatrixPtr(new Matrix());
            (*minusAllHKL)[0] = -1;
            (*minusAllHKL)[5] = -1;
            (*minusAllHKL)[10] = -1;
            
            return minusAllHKL;
        }

    }
    
    return MatrixPtr(new Matrix());
}

/**
 Hard-coding the number of matrix ambiguities for the twinnable space groups.
 */
int Reflection::ambiguityCount()
{
    if (spgNum >= 195 && spgNum <= 199)
        return 2;
    
    if (spgNum >= 168 && spgNum <= 173)
        return 2;
    
    if (spgNum == 146)
        return 2;
    
    if (spgNum >= 149 && spgNum <= 154)
        return 3;
    
    if (spgNum >= 143 && spgNum <= 145)
        return 4;
    
    return 1;
}

void Reflection::setSpaceGroup(CSym::CCP4SPG *ccp4spg)
{
    spgNum = ccp4spg->spg_num;
    spaceGroup = ccp4spg;
}

/**
 Takes the h,k,l values, flips the Miller and folds back onto the asymmetric unit. If the Miller indices are the same, it is twinned. (Not suitable for non-hemihedral twinning at the moment)
 */
bool Reflection::isTwinned()
{
    if (millerCount() > 0)
    {
        int h, k, l = 0;
        int trueH = miller(0)->getH();
        int trueK = miller(0)->getK();
        int trueL = miller(0)->getL();
        
        miller(0)->flip();
        MillerPtr firstMiller = miller(0);
        ccp4spg_put_in_asu(spaceGroup, firstMiller->getH(),
                           firstMiller->getK(), firstMiller->getL(), &h, &k, &l);
        miller(0)->flip();
        
        return (h == trueH && k == trueK && l == trueL);

    }
    else
    {
        return false;
    }
}

void Reflection::setSpaceGroup(int spaceGroupNum)
{
    spgNum = spaceGroupNum;
    spaceGroup = ccp4spg_load_by_ccp4_num(spgNum);
}

/** Easy to search integer for finding a particular Miller index. Should be replaced with the normal CCP4 sort order if this does not reduce speed too much */
int Reflection::reflectionIdForMiller(vec miller)
{
    int h = miller.h;
    int k = miller.k;
    int l = miller.l;
    
    int index = (h + OFFSET) * pow((double) MULTIPLIER, (int) 2)
    + (k + OFFSET) * MULTIPLIER + (l + OFFSET);
    
    return index;
}

/**
 Generates twinned/untwinned reflection IDs for easy searching.
 */
void Reflection::generateReflectionIds()
{
    if (millerCount() == 0)
    {
        std::cout << "Warning! Miller count is 0" << std::endl;
    }
    
    int h = miller(0)->getH();
    int k = miller(0)->getK();
    int l = miller(0)->getL();
    
    vec hkl = new_vector(h, k, l);
    
    for (int i = 0; i < ambiguityCount(); i++)
    {
        MatrixPtr ambiguityMat = matrixForAmbiguity(i);
        
        int asuH, asuK, asuL;
        
        ccp4spg_put_in_asu(spaceGroup, hkl.h, hkl.k, hkl.l, &asuH, &asuK, &asuL);
        vec hklAsu = new_vector(asuH, asuK, asuL);
        
        ambiguityMat->multiplyVector(&hklAsu);
       
        int newId = reflectionIdForMiller(hklAsu);
        
        reflectionIds.push_back(newId);
    }
}

void Reflection::setUnitCell(float *theUnitCell)
{
    this->unitCell[0] = theUnitCell[0];
    this->unitCell[1] = theUnitCell[1];
    this->unitCell[2] = theUnitCell[2];
    this->unitCell[3] = theUnitCell[3];
    this->unitCell[4] = theUnitCell[4];
    this->unitCell[5] = theUnitCell[5];
  
    setupUnitCell = true;
}

Reflection::Reflection(float *unitCell, CSym::CCP4SPG *newspg)
{
    spgNum = 0;
    
    if (newspg != NULL)
        setSpaceGroup(newspg->spg_num);

    if (unitCell != NULL)
    {
        setUnitCell(unitCell);
    }
    
    // TODO Auto-generated constructor stub
    
    refIntensity = 0;
    refSigma = 0;
    refl_id = 0;
    inv_refl_id = 0;
    resolution = 0;
    activeAmbiguity = 0;
    flagScaled = false;
}

MillerPtr Reflection::acceptedMiller(int num)
{
    int accepted = 0;
    
    for (int i = 0; i < millers.size(); i++)
    {
        if (miller(i)->accepted() && accepted == num)
            return miller(i);
        
        if (miller(i)->accepted())
            accepted++;
    }
    
    return MillerPtr();
}


MillerPtr Reflection::miller(int i)
{
    return millers[i];
}

void Reflection::addMiller(MillerPtr miller)
{
    miller->setResolution(resolution);
    millers.push_back(miller);
    
    if (reflectionIds.size() == 0)
    {
        generateReflectionIds();
    }
}

bool Reflection::betweenResolutions(double lowAngstroms, double highAngstroms)
{
    double minD, maxD = 0;
    StatisticsManager::convertResolutions(lowAngstroms,
                                          highAngstroms, &minD, &maxD);
    
    if (resolution > maxD || resolution < minD)
        return false;
    
    return true;
}

int Reflection::millerCount()
{
    return (int)millers.size();
}

void Reflection::removeMiller(int index)
{
    millers.erase(millers.begin() + index);
}

Reflection *Reflection::copy()
{
    return copy(false);
}

Reflection *Reflection::copy(bool copyMillers)
{
    Reflection *newReflection = new Reflection();
    
    newReflection->spgNum = spgNum;
    newReflection->activeAmbiguity = activeAmbiguity;
    newReflection->reflectionIds = reflectionIds;
    newReflection->refIntensity = refIntensity;
    newReflection->refSigma = refSigma;
    newReflection->refl_id = refl_id;
    newReflection->inv_refl_id = inv_refl_id;
    newReflection->resolution = resolution;
    
    for (int i = 0; i < millerCount(); i++)
    {
        MillerPtr newMiller;
        
        if (copyMillers)
        {
            newMiller = miller(i)->copy();
        }
        else
        {
            newMiller = miller(i);
        }
        
        newReflection->addMiller(newMiller);
    }

    return newReflection;
}

double Reflection::meanIntensity(int start, int end)
{
    if (end == 0)
        end = acceptedCount();
    
    double total_intensity = 0;
    double weight = 0;
    
    for (int i = start; i < end; i++)
    {
        MillerPtr miller = this->acceptedMiller(i);
        
        total_intensity += miller->intensity();
        weight += 1;
    }
    
    total_intensity /= weight;
    
    return total_intensity;
}

double Reflection::meanIntensityWithExclusion(std::string *filename, int start, int end)
{
    if (filename == NULL)
        return meanIntensity(start, end);
    
    if (end == 0)
        end = acceptedCount();
    double total_intensity = 0;
    double weight = 0;
    int accepted = acceptedCount();
    
    for (int i = start; i < end; i++)
    {
        MillerPtr miller = this->acceptedMiller(i);
        
        if (accepted > 2 && miller->getFilename() == *filename)
            continue;
            
        total_intensity += miller->intensity();
        weight += 1;
    }
    
    total_intensity /= weight;
    
    return total_intensity;
}
double Reflection::meanWeight()
{
    int num = (int)millers.size();
    int count = 0;
    
    double total_weight = 0;
    
    for (int i = 0; i < num; i++)
    {
        MillerPtr miller = millers[i];
        
        if (miller->accepted())
        {
            total_weight += 1 / pow(miller->getSigma(), 2);
            count++;
        }
    }
    
    total_weight /= count;
    
    return total_weight;
}

double Reflection::meanSigma()
{
    int num = (int)millers.size();
    int count = 0;
    
    double total_sigi = 0;
    
    for (int i = 0; i < num; i++)
    {
        MillerPtr miller = millers[i];
        
        if (miller->accepted())
        {
            total_sigi += miller->getSigma();
            count++;
        }
    }
    
    total_sigi /= count;
    
    return total_sigi;
}

void Reflection::anglesAsRadians(double *alpha, double *beta, double *gamma)
{
    *alpha = unitCell[3] * M_PI / 180;
    *beta = unitCell[4] * M_PI / 180;
    *gamma = unitCell[5] * M_PI / 180;
}

double Reflection::unitCellVolume()
{
    double alpha, beta, gamma;
    anglesAsRadians(&alpha, &beta, &gamma);
    
    double a = unitCell[0];
    double b = unitCell[1];
    double c = unitCell[2];
    
    double cosAlpha = cos(alpha);
    double cosBeta = cos(beta);
    double cosGamma = cos(gamma);
    
    double abc = a * b * c;
    double add = 1 - pow(cosAlpha, 2) - pow(cosBeta, 2) - pow(cosGamma, 2) +
                     2 * cosAlpha * cosBeta * cosGamma;
    
    double volume = abc * sqrt(add);
    
    return volume;
}

void Reflection::realSpaceAxes(vec *avec, vec *bvec, vec *cvec)
{
    double alpha, beta, gamma;
    anglesAsRadians(&alpha, &beta, &gamma);
    double cosAlpha = cos(alpha);
    double cosBeta = cos(beta);
    double sinGamma = sin(gamma);
    double cosGamma = cos(gamma);
    
    double a = unitCell[0];
    double b = unitCell[1];
    double c = unitCell[2];
    
    double a0 = a;
    double b0 = b * cosBeta;
    double b1 = b * sinGamma;
    double c0 = c * cosBeta;
    double c1 = (b * c * cosAlpha - b * cosGamma * cosBeta) / (b * sinGamma);
    double c2 = sqrt(c * c - c1 * c1 - c0 * c0);
    
    *avec = new_vector(a0, 0, 0);
    *bvec = new_vector(b0, b1, 0);
    *cvec = new_vector(c0, c1, c2);
}

void Reflection::reciprocalAxes(vec *aStar, vec *bStar, vec *cStar)
{
    double volume = unitCellVolume();
    vec avec, bvec, cvec;
    realSpaceAxes(&avec, &bvec, &cvec);
    double alpha, beta, gamma;
    anglesAsRadians(&alpha, &beta, &gamma);
    double sinAlpha = sin(alpha);
    double sinBeta = sin(beta);
    double sinGamma = sin(gamma);
    
    *aStar = cross_product_for_vectors(bvec, cvec);
    double scalar = sinAlpha / volume;
    multiply_vector(aStar, scalar);

    *bStar = cross_product_for_vectors(avec, cvec);
    scalar = sinBeta / volume;
    multiply_vector(bStar, scalar);

    *cStar = cross_product_for_vectors(avec, bvec);
    scalar = sinGamma / volume;
    multiply_vector(cStar, scalar);
}

void Reflection::calculateResolution(MtzManager *mtz)
{
    int h = millers[0]->getH();
    int k = millers[0]->getK();
    int l = millers[0]->getL();
    
    // calculate resolution
    
    vec aStar, bStar, cStar;
    reciprocalAxes(&aStar, &bStar, &cStar);
    
    double dotAA = dot_product_for_vectors(aStar, aStar);
    double dotBB = dot_product_for_vectors(bStar, bStar);
    double dotCC = dot_product_for_vectors(cStar, cStar);

    double dotAB = dot_product_for_vectors(aStar, bStar);
    double dotAC = dot_product_for_vectors(aStar, cStar);
    double dotBC = dot_product_for_vectors(bStar, cStar);

    double term1 = h * h * dotAA;
    double term2 = k * k * dotBB;
    double term3 = l * l * dotCC;
    
    double term4 = 2 * h * k * dotAB;
    double term5 = 2 * h * l * dotAC;
    double term6 = 2 * k * l * dotBC;
    
    double dStarSqr = term1 + term2 + term3 + term4 + term5 + term6;
    this->resolution = sqrt(dStarSqr);
    
    for (int i = 0; i < millerCount(); i++)
    {
        miller(i)->setResolution(resolution);
    }
    
//    std::cout << h << " " << k << " " << l << " " << resolution << std::endl;
}


void Reflection::incrementAmbiguity()
{
    int count = ambiguityCount();
    int newActive = activeAmbiguity + 1;
    activeAmbiguity = (newActive % count);
}

void Reflection::setFlipAsActiveAmbiguity()
{
    setFlip(activeAmbiguity);
}

void Reflection::resetFlip()
{
    for (int j = 0; j < millerCount(); j++)
    {
        miller(j)->setFlipMatrix(MatrixPtr(new Matrix()));
    }
}

void Reflection::setFlip(int i)
{
    for (int j = 0; j < millerCount(); j++)
    {
        MatrixPtr ambiguityMat = matrixForAmbiguity(i);
        
        miller(j)->setFlipMatrix(ambiguityMat);
    }
}


void Reflection::reflectionDescription()
{
    int acceptedCount = 0;
    
    std::ostringstream logged;
    
    for (int i = 0; i < millerCount(); i++)
    {
        MillerPtr miller = this->miller(i);
        logged << miller->getH() << "\t" << miller->getK() << "\t" << miller->getL() << "\t"
        << miller->getRawIntensity() << "\t"
        << "\t" << miller->getSigma() << "\t" << miller->getFilename()
        << std::endl;
        if (miller->accepted())
            acceptedCount++;
    }
    logged << std::endl;
    
    std::cout << logged.str();
}

void Reflection::clearMillers()
{
    millers.clear();
    vector<MillerPtr>().swap(millers);
    
}

Reflection::~Reflection()
{
    clearMillers();
}

int Reflection::indexForReflection(int h, int k, int l, CSym::CCP4SPG *spgroup,
                               bool inverted)
{
    int _h, _k, _l;
    
    ccp4spg_put_in_asu(spgroup, h, k, l, &_h, &_k, &_l);
    
    int multiplier = MULTIPLIER;
    int offset = OFFSET;
    
    int index = (_h + OFFSET) * pow((double) MULTIPLIER, (int) 2)
    + (_k + OFFSET) * MULTIPLIER + (_l + OFFSET);
    if (inverted)
        index = (_h + OFFSET) * pow((double) MULTIPLIER, (int) 2)
        + (_l + OFFSET) * MULTIPLIER + (_k + OFFSET);
    
    if (spgroup->spg_num == 197)
    {
        if (inverted == 0)
            index = (_h + offset) * pow((double) multiplier, (int) 2)
            + (_k + offset) * multiplier + (_l + offset);
        else
            index = (_h + offset) * pow((double) multiplier, (int) 2)
            + (_l + offset) * multiplier + (_k + offset);
    }
    else if (spgroup->spg_num == 146)
    {
        if (inverted == 0)
            index = (_h + offset) * pow((double) multiplier, (int) 2)
            + (_k + offset) * multiplier + (_l + offset);
        else
            index = (_k + offset) * pow((double) multiplier, (int) 2)
            + (_h + offset) * multiplier + ((-_l) + offset);
    }
    
    return index;
}

int Reflection::rejectCount()
{
    int count = 0;
    
    for (int i = 0; i < millerCount(); i++)
    {
        if (miller(i)->isRejected())
        {
            count++;
        }
    }
    
    return count;
}

int Reflection::acceptedCount()
{
    int count = 0;
    
    for (int i = 0; i < millerCount(); i++)
    {
        if (miller(i)->accepted())
        {
            count++;
        }
    }
    
    return count;
}

bool Reflection::anyAccepted()
{
    for (int i = 0; i < millers.size(); i++)
    {
        if (miller(i)->accepted())
            return true;
    }
    
    return false;
}

void Reflection::printDescription()
{
    std::cout << "Mean intensity " << this->meanIntensity() << std::endl;
    std::cout << "Miller corrected intensities: ";
    
    for (int i = 0; i < millerCount(); i++)
    {
        std::cout << miller(i)->intensity() << ", ";
    }
    
    std::cout << std::endl;
}


void Reflection::setAdditionalWeight(double weight)
{
    for (int i = 0; i < millerCount(); i++)
    {
        miller(i)->setAdditionalWeight(weight);
    }
}

void Reflection::scale(double scale)
{
    for (int i = 0; i < millerCount(); i++)
    {
        miller(i)->applyScaleFactor(scale);
    }
}

