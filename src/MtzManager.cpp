#include "MtzManager.h"

#include <string>
#include <iostream>
#include <csymlib.h>
#include "Vector.h"
#include "FileReader.h"
#include <cerrno>
#include <fstream>
#include <iostream>
#include <boost/variant.hpp>
#include "StatisticsManager.h"
#include "Matrix.h"
#include "Holder.h"
#include "Miller.h"
#include <ccp4_parser.h>

#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "definitions.h"
#include "FileParser.h"


using namespace CMtz;

using namespace CSym;

MtzManager *MtzManager::referenceManager;
double MtzManager::highRes = 3.5;
double MtzManager::lowRes = 0;

std::string MtzManager::filenameRoot()
{
    vector<std::string> components = FileReader::split(filename, '.');
    
    std::string root = "";
    
    for (int i = 0; i < components.size() - 1; i++)
    {
        root += components[i] + ".";
    }
    
    return root.substr(0, root.size() - 1);
}

double MtzManager::extreme_index(MTZ *mtz, int max)
{
    double extreme = 0;
    MTZCOL *col_h = MtzColLookup(mtz, "H");
    MTZCOL *col_k = MtzColLookup(mtz, "K");
    MTZCOL *col_l = MtzColLookup(mtz, "L");
    
    if (max == 0)
    {
        extreme = col_h->min;
        extreme = (extreme > col_k->min) ? col_k->min : extreme;
        extreme = (extreme > col_l->min) ? col_l->min : extreme;
    }
    else
    {
        extreme = col_h->max;
        extreme = (extreme < col_k->max) ? col_k->max : extreme;
        extreme = (extreme < col_l->max) ? col_l->max : extreme;
    }
    
    return extreme;
}

void MtzManager::findMultiplier(MTZ *mtz, int *multiplier, int *offset)
{
    int hkl_min = extreme_index(mtz, 0);
    int hkl_max = extreme_index(mtz, 1);
    
    *offset = (hkl_min < 0) ? -hkl_min : 0;
    
    *multiplier = hkl_max - ((hkl_min < 0) ? hkl_min : 0);
    
    *multiplier = MULTIPLIER;
    *offset = OFFSET;
}

void MtzManager::hkls_for_reflection(MTZ *mtz, float *adata, int *h, int *k,
                                     int *l, int *multiplier, int *offset)
{
    if (multiplier == 0)
        findMultiplier(mtz, multiplier, offset);
    
    MTZCOL *col_h = MtzColLookup(mtz, "H");
    MTZCOL *col_k = MtzColLookup(mtz, "K");
    MTZCOL *col_l = MtzColLookup(mtz, "L");
    
    *h = adata[col_h->source - 1];
    *k = adata[col_k->source - 1];
    *l = adata[col_l->source - 1];
}

int MtzManager::index_for_reflection(int h, int k, int l, bool inverted)
{
    int _h = h;
    int _k = k;
    int _l = l;
    
    if (inverted)
    {
        if (low_group->spg_num == 197)
        {
            _h = k;
            _k = h;
            _l = -l;
        }
        
        if (low_group->spg_num == 146)
        {
            _h = k;
            _k = h;
            _l = -l;
        }
        
        if (low_group->spg_num == 4 || low_group->spg_num == 21)
        {
            _h = l;
            _k = k;
            _l = h;
        }
    }
    
    ccp4spg_put_in_asu(low_group, _h, _k, _l, &h, &k, &l);
    
    int index = (h + OFFSET) * pow((double) MULTIPLIER, (int) 2)
    + (k + OFFSET) * MULTIPLIER + (l + OFFSET);
    
    return index;
}

bool MtzManager::reflection_comparison(Reflection *i, Reflection *j)
{
    return (i->getReflId() < j->getReflId());
}

void MtzManager::insertionSortReflections(void)
{
    std::sort(reflections.begin(), reflections.end(), reflection_comparison);
}

void MtzManager::sortLastReflection(void)
{
    int i, j;
    Reflection *tmp;
    
    i = (int)reflections.size() - 1;
    j = i;
    tmp = reflections[i];
    while (j > 0 && tmp->getReflId() < reflections[j - 1]->getReflId())
    {
        reflections[j] = reflections[j - 1];
        j--;
    }
    reflections[j] = tmp;
    
    return;
}

MtzManager::MtzManager(std::string nameOfFile)
{
    lastReference = NULL;
    failedCount = 0;
    filename = nameOfFile;
    reflections.resize(0);
    low_group = NULL;
    inverse = false;
    flipped = false;
    finalised = false;
    scoreType = ScoreTypeCorrelation;
    trust = TrustLevelBad;
    freePass = false;
    rejected = false;
    previousReference = NULL;
    previousAmbiguity = -1;
    activeAmbiguity = 0;
    bFactor = 0;
    requiredPriority = FileParser::getKey("VERBOSITY_LEVEL", 0);
    scale = 0;
    
    matrix = MatrixPtr();
}

MtzPtr MtzManager::copy()
{
    MtzPtr newManager = MtzPtr(new MtzManager());
    newManager->filename = filename;
    double lowNum = low_group->spg_num;
    newManager->setSpaceGroup(lowNum);
    newManager->flipped = flipped;
    newManager->matrix = matrix;
    
    for (int i = 0; i < 3; i++)
    {
        newManager->cellDim[i] = cellDim[i];
        newManager->cellAngles[i] = cellAngles[i];
        
    }
    
    for (int i = 0; i < reflections.size(); i++)
    {
        Reflection *newReflection = reflections[i]->copy();
        newManager->reflections.push_back(newReflection);
    }
    
    return newManager;
}


void MtzManager::setUnitCell(double a, double b, double c, double alpha, double beta,
                             double gamma)
{
    cellDim[0] = a;
    cellDim[1] = b;
    cellDim[2] = c;
    cellAngles[0] = alpha;
    cellAngles[1] = beta;
    cellAngles[2] = gamma;
}

void MtzManager::setUnitCell(vector<double> unitCell)
{
    cellDim[0] = unitCell[0];
    cellDim[1] = unitCell[1];
    cellDim[2] = unitCell[2];
    cellAngles[0] = unitCell[3];
    cellAngles[1] = unitCell[4];
    cellAngles[2] = unitCell[5];
}

void MtzManager::clearReflections()
{
    for (int i = 0; i < reflections.size(); i++)
    {
        delete reflections[i];
    }
    
    reflections.clear();
    vector<Reflection *>().swap(reflections);
}

void MtzManager::removeReflection(int i)
{
    delete reflections[i];
    
    reflections.erase(reflections.begin() + i);
}

void MtzManager::addReflection(Reflection *reflection)
{
    Reflection *newReflection = NULL;
    int insertionPoint = findReflectionWithId(reflection->getReflId(), &newReflection, true);
    
    reflections.insert(reflections.begin() + insertionPoint, reflection);
    
 //   reflections.push_back(reflection);
 //   this->sortLastReflection();
}

void MtzManager::addReflections(vector<Reflection *>reflections)
{
    for (int i = 0; i < reflections.size(); i++)
    {
        addReflection(reflections[i]);
    }
}

void MtzManager::setMatrix(double *components)
{
    matrix = MatrixPtr(new Matrix(components));
}

void MtzManager::setMatrix(MatrixPtr newMat)
{
    matrix = newMat;
}

void MtzManager::getUnitCell(double *a, double *b, double *c, double *alpha,
                             double *beta, double *gamma)
{
    *a = cellDim[0];
    *b = cellDim[1];
    *c = cellDim[2];
    *alpha = cellAngles[0];
    *beta = cellAngles[1];
    *gamma = cellAngles[2];
}

void MtzManager::copySymmetryInformationFromManager(MtzPtr toCopy)
{
    double a, b, c, alpha, beta, gamma;
    toCopy->getUnitCell(&a, &b, &c, &alpha, &beta, &gamma);
    this->setUnitCell(a, b, c, alpha, beta, gamma);
    
    std::cout << "Setting unit cell to " << a << ", " << b << ", " << c << ", " << alpha << ", " << beta << ", " << gamma <<std::endl;
    
    int spgnum = toCopy->getLowGroup()->spg_ccp4_num;
    CCP4SPG *spg = ccp4spg_load_by_ccp4_num(spgnum);
    this->setLowGroup(spg);
}

void MtzManager::loadReflections(bool calculated)
{
    if (filename.length() == 0)
    {
        std::cerr
        << "Cannot load reflections as no filename has been specified."
        << std::endl;
        return;
    }
    
    std::cout << "Loading " << filename << std::endl;
    
    MTZ *mtz = MtzGet(filename.c_str(), 0);
    
    int spgnum = MtzSpacegroupNumber(mtz);
    
    low_group = ccp4spg_load_by_ccp4_num(spgnum);
    
    if (low_group == NULL)
    {
        return;
    }
    
    if (mtz == NULL)
        return;
    
    float *refldata = (float *) CCP4::ccp4_utils_malloc(
                                                        (mtz->ncol_read + 1) * mtz->nref_filein * sizeof(float));
    memset(refldata, '\0',
           (mtz->ncol_read + 1) * mtz->nref_filein * sizeof(float));
    
    float *adata = (float *) CCP4::ccp4_utils_malloc(
                                                     (mtz->ncol_read) * sizeof(float));
    
    MtzRrefl(mtz->filein, mtz->ncol_read * mtz->nref_filein, refldata);
    
    std::string fObsLab = FileParser::getKey("LABIN_AMPLITUDES", std::string("F"));
    std::string fCalcLab = FileParser::getKey("LABOUT_CALC_AMPLITUDES", std::string("FC"));
    std::string sig1Lab = FileParser::getKey("LABIN_SIGMAS", std::string("SIGF"));
    std::string phiLab = FileParser::getKey("LABIN_PHASES", std::string("PHIC"));
    
    std::string chosenF = calculated ? fCalcLab : fObsLab;
    std::string otherF = calculated ? fObsLab : fCalcLab;
    
    MTZCOL *col_f = MtzColLookup(mtz, chosenF.c_str());
    MTZCOL *other_f = MtzColLookup(mtz, otherF.c_str());
    
    if (col_f == NULL)
    {
        std::cout << "Warning: intensity column not labelled " + chosenF + ", please specify amplitude label in keyword LABIN_AMPLITUDES" << std::endl;
        exit(1);
    }
    
    if (other_f == NULL && calculated)
    {
        std::cout << "Warning: intensity column not labelled " + otherF + ", please specify amplitude label in keyword LABIN_AMPLITUDES" << std::endl;
        exit(1);
    }
        
    MTZCOL *col_sigf = MtzColLookup(mtz, sig1Lab.c_str());
    
    if (col_sigf == NULL)
    {
        std::cout << "Warning: sigma column not labelled " + sig1Lab + ", please specify sigf label in keyword LABIN_SIGMAS" << std::endl;
        exit(1);
    }
    
    MTZCOL *col_phic = MtzColLookup(mtz, phiLab.c_str());
    
    if (col_phic == NULL)
    {
        std::cout << "Warning: phase column not labelled " + phiLab + ", please specify phase label in keyword LABIN_PHASES" << std::endl;
        exit(1);
    }
    
    int multiplier = MULTIPLIER;
    
    MTZXTAL **xtals = MtzXtals(mtz);
    float cell[6];
    double coefhkl[6];
    ccp4_lrcell(xtals[0], cell);
    MtzHklcoeffs(cell, coefhkl);
    
    setUnitCell(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]);
    
    int num = 0;
    
    for (int i = 0; i < mtz->nref_filein * mtz->ncol_read; i += mtz->ncol_read)
    {
        memcpy(adata, &refldata[i], mtz->ncol_read * sizeof(float));
        
        int h;
        int k;
        int l;
        int offset;
        
        hkls_for_reflection(mtz, adata, &h, &k, &l, &multiplier, &offset);
        
        int reflection_reflection_index = index_for_reflection(h, k, l, false);
       
        float intensity = adata[col_f->source - 1];
        
        float sigma = 1;
        if (col_sigf != NULL)
            sigma = adata[col_sigf->source - 1];
        
        if (sigma == -100)
        {
            sigma = isnan(' ');
        }
        
        float wavelength = 1;
        float partiality = 1;
        float shiftX = 0;
        float shiftY = 0;
        float otherF = 0;
        float phase = 0;
        
        if (other_f != NULL)
            otherF = adata[other_f->source - 1];
        if (col_phic != NULL)
            phase = adata[col_phic->source - 1];
        
        MillerPtr miller = MillerPtr(new Miller(this, h, k, l));
        miller->setData(intensity, sigma, partiality, wavelength);
        miller->setCountingSigma(sigma);
        miller->setFilename(filename);
        miller->setPhase(phase);
        miller->setFc(otherF);
        miller->setShift(std::make_pair(shiftX, shiftY));
        miller->matrix = this->matrix;
        
        Reflection *prevReflection;
        
        this->findReflectionWithId(reflection_reflection_index, &prevReflection);
        
        if (prevReflection != NULL)
        {
            /** Exclude unobserved reflections by checking for nan */
            if (adata[col_f->source - 1] == adata[col_f->source - 1])
            {
                prevReflection->addMiller(miller); // TODO
                miller->setParent(prevReflection);
            }
            
            // reflection is a repeat so set flag.
        }
        
        if (prevReflection == NULL)
        {
            Reflection *newReflection = new Reflection();
            reflections.push_back(newReflection);
            newReflection->setUnitCell(cell);
            newReflection->setSpaceGroup(low_group);
            
            reflections[num]->addMiller(miller);
            miller->setParent(reflections[num]);
            reflections[num]->calculateResolution(this);
            
            this->sortLastReflection();
            
            num++;
        }
    }
    
    free(refldata);
    free(adata);
    
    //	insertionSortReflections();
    
    std::ostringstream log;
    
    log << "Loaded " << mtz->nref_filein << " reflections." << std::endl;
    
    MtzFree(mtz);
}

void MtzManager::getRefReflections(vector<Reflection *> *refPointer,
                               vector<Reflection *> *matchPointer)
{
    if (lastReference != referenceManager)
    {
        refReflections.clear();
        matchReflections.clear();
        
        for (int i = 0; i < reflectionCount(); i++)
        {
            int reflId = reflection(i)->getReflId();
            
            Reflection *refReflection = NULL;
            referenceManager->findReflectionWithId(reflId, &refReflection);
            
            if (refReflection != NULL)
            {
                matchReflections.push_back(reflection(i));
                refReflections.push_back(refReflection);
            }
        }
        
        lastReference = referenceManager;
    }
    
    refPointer->reserve(refReflections.size());
    refPointer->insert(refPointer->begin(), refReflections.begin(),
                       refReflections.end());
    
    matchPointer->reserve(matchReflections.size());
    matchPointer->insert(matchPointer->begin(), matchReflections.begin(),
                         matchReflections.end());
    
}

void MtzManager::setReference(MtzManager *reference)
{
    MtzManager::referenceManager = reference;
}

void MtzManager::setFilename(std::string name)
{
    filename = name;
}

std::string MtzManager::getFilename(void)
{
    return filename;
}

void MtzManager::setSpaceGroup(int spgnum)
{
    low_group = ccp4spg_load_by_ccp4_num(spgnum);
}

int MtzManager::findReflectionWithId(int refl_id, Reflection **reflection, bool insertionPoint)
{
    if (reflectionCount() == 0)
    {
        *reflection = NULL;
        return 0;
    }
    
    int lower = 0;
    int higher = reflectionCount() - 1;
    int new_bound = (higher + lower) / 2;
    
    if ((refl_id < this->reflection(lower)->getReflId())
        || (refl_id > this->reflection(higher)->getReflId()))
    {
        if (insertionPoint)
        {
            if (refl_id < this->reflection(lower)->getReflId())
                return 0;
            
            if (refl_id > this->reflection(higher)->getReflId())
                return reflectionCount();
        }
        
        *reflection = NULL;
        return -1;
    }
    
    while (this->reflection(new_bound)->getReflId() != refl_id)
    {
        if (new_bound == higher || new_bound == lower)
        {
            if (this->reflection(higher)->getReflId() == refl_id)
            {
                (*reflection) = this->reflection(higher);
                return -1;
            }
            
            if (this->reflection(lower)->getReflId() == refl_id)
            {
                (*reflection) = this->reflection(lower);
                return -1;
            }
            
            if (insertionPoint)
            {
                int lowest = 0;
                
                int start = lower - 2 >= 0 ? lower - 2 : 0;
                
                for (int i = start; i < higher + 2 && i < reflectionCount(); i++)
                {
                    if (this->reflection(i)->getReflId() < refl_id)
                        lowest = i;
                }
                
                return lowest + 1;
            }
            
            *reflection = NULL;
            return -1;
        }
        
        if (this->reflection(new_bound)->getReflId() > refl_id)
        {
            higher = new_bound;
        }
        else if (this->reflection(new_bound)->getReflId() < refl_id)
        {
            lower = new_bound;
        }
        
        new_bound = (higher + lower) / 2;
    }
    
    (*reflection) = this->reflection(new_bound);
    
    return -1;
}

void MtzManager::findCommonReflections(MtzManager *other,
                                       vector<Reflection *> &reflectionVector1, vector<Reflection *> &reflectionVector2,
                                       int *num, bool twinnedOnly)
{
    matchReflections.clear();
    refReflections.clear();
    
    previousReference = other;
    previousAmbiguity = activeAmbiguity;
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        if (twinnedOnly && !reflection(i)->isTwinned())
            continue;
        
        int refl_id = reflection(i)->getReflId();
        
        Reflection *otherReflection = NULL;
        
        other->findReflectionWithId(refl_id, &otherReflection);
        
        if (otherReflection != NULL && otherReflection->millerCount() > 0)
        {
            reflectionVector1.push_back(reflection(i));
            matchReflections.push_back(reflection(i));
            reflectionVector2.push_back(otherReflection);
            refReflections.push_back(otherReflection);
        }
    }
    
    if (num != NULL)
    {
        *num = (int)reflectionVector1.size();
    }
}

void MtzManager::applyScaleFactorsForBins(int binCount)
{
    double highRes = FileParser::getKey("HIGH_RESOLUTION", 0.0);
    double lowRes = FileParser::getKey("LOW_RESOLUTION", 0.0);
    
    vector<double> bins;
    StatisticsManager::generateResolutionBins(lowRes, highRes, binCount, &bins);
    
    for (int shell = 0; shell < bins.size() - 1; shell++)
    {
        double low = bins[shell];
        double high = bins[shell + 1];
        
        vector<Reflection *> refReflections, imgReflections;
        
        this->findCommonReflections(referenceManager, imgReflections, refReflections,
                                    NULL);
        
        double weights = 0;
        double refMean = 0;
        double imgMean = 0;
        int count = 0;
        
        for (int i = 0; i < imgReflections.size(); i++)
        {
            if (!imgReflections[i]->anyAccepted())
                continue;
            
            if (imgReflections[i]->betweenResolutions(low, high))
            {
                if (imgReflections[i]->isTwinned())
                {
                    weights++;
                    double refIntensity = refReflections[i]->meanIntensity();
                    double imgIntensity = imgReflections[i]->meanIntensity();
                    
                    if (refIntensity != refIntensity || imgIntensity != imgIntensity)
                        continue;
                    
                    refMean += refIntensity;
                    imgMean += imgIntensity;
                    count++;
                }
            }
        }
        
        refMean /= weights;
        imgMean /= weights;
        
        double ratio = refMean / imgMean;
        
        applyScaleFactor(ratio, low, high);
    }
}

double MtzManager::gradientAgainstManager(MtzManager &otherManager,
                                          bool leastSquares, double lowRes, double highRes)
{
    vector<Reflection *> reflections1;
    vector<Reflection *> reflections2;
    
    vector<double> ints1, ints2;
    int num = 0;
    
    double minD = 0;
    double maxD = 0;
    StatisticsManager::convertResolutions(lowRes, highRes, &minD, &maxD);
    
    MtzManager::findCommonReflections(&otherManager, reflections1, reflections2, &num);
    
    if (num <= 1)
        return 1;
    
    double x_squared = 0;
    double x_y = 0;
    
    for (int i = 0; i < num; i++)
    {
        for (int j = 0; j < reflections1[i]->millerCount(); j++)
        {
            if (!reflections1[i]->miller(j)->accepted())
                continue;
            
            if (reflections1[i]->getResolution() < minD
                || reflections1[i]->getResolution() > maxD)
                continue;
            
            double mean1 = reflections1[i]->miller(j)->intensity();
            double mean2 = reflections2[i]->meanIntensity();
            
            double part1 = reflections1[i]->meanWeight();
            
            if (mean1 != mean1 || mean2 != mean2)
                continue;
            
            x_squared += mean1 * mean1 * part1;
            x_y += mean1 * mean2 * part1;
            
     //       if (reflections1[i]->miller(j)->getPartiality() < 0.5)
      //          continue;
            
            ints1.push_back(mean1);
            ints2.push_back(mean2);
            
        }
    }
    
    double grad = x_y / x_squared;
    
    if (grad < 0)
        grad = -1;
    
    if (leastSquares)
        grad = minimize_gradient_between_vectors(&ints1, &ints2);
    
    return grad;
}

void MtzManager::clearScaleFactor()
{
    for (int i = 0; i < reflections.size(); i++)
    {
        for (int j = 0; j < reflections[i]->millerCount(); j++)
        {
            reflections[i]->miller(j)->setScale(1);
            reflections[i]->miller(j)->setBFactor(0);
        }
    }
}

void MtzManager::makeScalesPermanent()
{
    for (int i = 0; i < reflections.size(); i++)
    {
        for (int j = 0; j < reflection(i)->millerCount(); j++)
        {
            reflection(i)->miller(j)->makeScalesPermanent();
        }
    }
}

void MtzManager::applyBFactor(double bFactor)
{
 //   if (bFactor < 0)
 //       bFactor = 0 - bFactor;
    
    for (int i = 0; i < reflections.size(); i++)
    {
        for (int j = 0; j < reflections[i]->millerCount(); j++)
        {
            reflection(i)->miller(j)->setBFactor(bFactor);
        }
    }
}

void MtzManager::applyScaleFactor(double scaleFactor,
                                  double lowRes, double highRes, bool absolute)
{
    if (scaleFactor == scaleFactor)
    {
        if (absolute)
            scale = scaleFactor;
        else
            scale *= scaleFactor;
    }
    
    for (int i = 0; i < reflections.size(); i++)
    {
        for (int j = 0; j < reflections[i]->millerCount(); j++)
        {
            if (reflection(i)->betweenResolutions(lowRes, highRes))
            {
                if (absolute)
                    reflection(i)->miller(j)->setScale(scale);
                else
                    reflection(i)->miller(j)->applyScaleFactor(scaleFactor);
                
                if (i == 0 && j == 0)
                    logged << reflection(i)->miller(j)->getScale() << std::endl;
            }
        }
    }
    
}

double MtzManager::averageIntensity(void)
{
    double total_intensity = 0;
    double total = 0;
    
    for (int i = 0; i < reflections.size(); i++)
    {
        double intensity = reflections[i]->meanIntensity();
        if (intensity == intensity)
        {
            double weight = reflections[i]->meanWeight();
            total_intensity += intensity * weight;
            total += weight;
        }
    }
    
    total_intensity /= total;
    return total_intensity;
}


void MtzManager::writeToFile(std::string newFilename, bool announce, bool shifts, bool includeAmbiguity)
{
    int columns = 7;
    if (includeAmbiguity)
        flipToActiveAmbiguity();
    
    if (shifts) columns += 2;
    
    float cell[6], wavelength;
    float *fdata = new float[columns];
    
    /* variables for symmetry */
    CCP4SPG *mtzspg = low_group;
    float rsm[192][4][4];
    char ltypex[2];
    
    /* variables for MTZ data structure */
    MTZ *mtzout;
    MTZXTAL *xtal;
    MTZSET *set;
    MTZCOL *colout[9];
    
    
    /*  Removed: General CCP4 initializations e.g. HKLOUT on command line */
    
    cell[0] = cellDim[0];
    cell[1] = cellDim[1];
    cell[2] = cellDim[2];
    cell[3] = cellAngles[0];
    cell[4] = cellAngles[1];
    cell[5] = cellAngles[2];
    wavelength = 1;
    
    mtzout = MtzMalloc(0, 0);
    ccp4_lwtitl(mtzout, "Written from Helen's XFEL tasks ", 0);
    mtzout->refs_in_memory = 0;
    mtzout->fileout = MtzOpenForWrite(newFilename.c_str());
    
    // then add symm headers...
    for (int i = 0; i < mtzspg->nsymop; ++i)
        CCP4::rotandtrn_to_mat4(rsm[i], mtzspg->symop[i]);
    strncpy(ltypex, mtzspg->symbol_old, 1);
    ccp4_lwsymm(mtzout, mtzspg->nsymop, mtzspg->nsymop_prim, rsm, ltypex,
                mtzspg->spg_ccp4_num, mtzspg->symbol_old, mtzspg->point_group);
    
    // then add xtals, datasets, cols
    xtal = MtzAddXtal(mtzout, "XFEL crystal", "XFEL project", cell);
    set = MtzAddDataset(mtzout, xtal, "Dataset", wavelength);
    colout[0] = MtzAddColumn(mtzout, set, "H", "H");
    colout[1] = MtzAddColumn(mtzout, set, "K", "H");
    colout[2] = MtzAddColumn(mtzout, set, "L", "H");
    colout[3] = MtzAddColumn(mtzout, set, "F", "F");
    colout[4] = MtzAddColumn(mtzout, set, "FC", "F");
    colout[5] = MtzAddColumn(mtzout, set, "SIGF", "Q");
    colout[6] = MtzAddColumn(mtzout, set, "PHIC", "P");
    
    int num = 0;
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        for (int j = 0; j < reflections[i]->millerCount(); j++)
        {
            if (!reflection(i)->miller(j))
                std::cout << "!miller(j) in mtz manager" << std::endl;
            
            if (reflections[i]->miller(j)->isRejected())
                continue;
            
            double intensity = reflections[i]->meanIntensity();
            double sigma = reflections[i]->miller(j)->getCountingSigma();
            double calcF = reflections[i]->miller(j)->getFc();
            double phases = reflections[i]->miller(j)->getPhase();
            
            if (calcF != calcF)
            {
                continue;
            }
            
            num++;
            
            int h = reflections[i]->miller(j)->getH();
            int k = reflections[i]->miller(j)->getK();
            int l = reflections[i]->miller(j)->getL();
            int _h, _k, _l;
            ccp4spg_put_in_asu(low_group, h, k, l, &_h, &_k, &_l);
            
            // set fdata
            fdata[0] = h;
            fdata[1] = k;
            fdata[2] = l;
            fdata[3] = calcF;
            fdata[4] = intensity;
            fdata[5] = sigma;
            fdata[6] = phases;
            
            if (shifts)
            {
                fdata[7] = reflections[i]->miller(j)->getShift().first;
                fdata[8] = reflections[i]->miller(j)->getShift().second;
            }
            
            ccp4_lwrefl(mtzout, fdata, colout, columns, num);
        }
    }
    
    MtzPut(mtzout, " ");
    MtzFree(mtzout);
    
    std::ostringstream logged;
    logged << "Written to file " << newFilename << std::endl;
    std::cout << logged.str();
    
    if (includeAmbiguity)
        resetFlip();
    
    delete [] fdata;
}


MtzManager::~MtzManager(void)
{
    for (int i = 0; i < reflections.size(); i++)
    {
        delete reflections[i];
    }
    
    reflections.clear();
    vector<Reflection *>().swap(reflections);
    
    if (low_group != NULL)
        ccp4spg_free(&low_group);
    
}

void MtzManager::description(void)
{
    logged << "Filename: " << filename << std::endl;
    logged << "Number of reflections: " << reflectionCount() << std::endl;
    logged << "Average intensity: " << this->averageIntensity() << std::endl;
    
    sendLog();
}

void MtzManager::sendLog(LogLevel priority)
{
    if (priority > requiredPriority)
        std::cout << logged.str();
    
    logged.str("");
    logged.clear();
}

int MtzManager::ambiguityCount()
{
    int count = reflection(0)->ambiguityCount();
    
    return count;
}

void MtzManager::incrementActiveAmbiguity()
{
    int count = ambiguityCount();
    
    activeAmbiguity++;
    activeAmbiguity = (activeAmbiguity % count);
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        reflection(i)->setActiveAmbiguity(activeAmbiguity);
    }
}

void MtzManager::flipToActiveAmbiguity()
{
    for (int i = 0; i < reflectionCount(); i++)
    {
        reflection(i)->setFlip(activeAmbiguity);
    }
    
    this->insertionSortReflections();
}


void MtzManager::resetFlip()
{
    for (int i = 0; i < reflectionCount(); i++)
    {
        reflection(i)->resetFlip();
    }
    
    this->insertionSortReflections();
}

void MtzManager::setActiveAmbiguity(int newAmbiguity)
{
    activeAmbiguity = newAmbiguity;
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        reflection(i)->setActiveAmbiguity(newAmbiguity);
    }
}

double MtzManager::maxResolution()
{
    double maxResolution = 0;
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        if (reflection(i)->getResolution() > maxResolution)
        {
            maxResolution = reflection(i)->getResolution();
        }
    }
    
    return maxResolution;
}


double MtzManager::correlation(bool silent, double lowResolution,
                               double highResolution)
{
    if (highResolution == -1)
        highResolution = maxResolution();
    
    double correlation = StatisticsManager::cc_pearson(this,
                                                       MtzManager::referenceManager, true, NULL,
                                                       NULL, lowResolution, highResolution, false);
    
    return correlation;
}

double MtzManager::rSplitWithManager(MtzManager *otherManager, bool printHits,
                                     bool silent, double lowRes, double highRes, int bins,
                                     vector<boost::tuple<double, double, double, int> > *correlations,
                                     bool shouldLog)
{
    StatisticsFunction *function = StatisticsManager::r_split;
    
    return statisticsWithManager(otherManager, function, printHits, silent,
                                 lowRes, highRes, bins, correlations, shouldLog);
}

double MtzManager::rSplitIntensityBinsWithManager(MtzManager *otherManager, bool printHits,
                                     bool silent, double lowInt, double highInt, int bins,
                                     vector<boost::tuple<double, double, double, int> > *correlations,
                                     bool shouldLog)
{
    StatisticsFunction *function = StatisticsManager::r_split_by_intensity;
    
    return statisticsWithManager(otherManager, function, printHits, silent,
                                 lowInt, highInt, bins, correlations, shouldLog);
}

double MtzManager::correlationWithManager(MtzManager *otherManager,
                                          bool printHits, bool silent, double lowRes, double highRes, int bins,
                                          vector<boost::tuple<double, double, double, int> > *correlations,
                                          bool shouldLog)
{
    StatisticsFunction *function = StatisticsManager::cc_pearson;
    
    return statisticsWithManager(otherManager, function, printHits, silent,
                                 lowRes, highRes, bins, correlations, shouldLog);
}

double MtzManager::statisticsWithManager(MtzManager *otherManager,
                                         StatisticsFunction *function, bool printHits, bool silent, double lowRes,
                                         double highRes, int bins,
                                         vector<boost::tuple<double, double, double, int> > *correlations,
                                         bool shouldLog)
{
    return statisticsWithManager(otherManager, function, NULL, RFactorNone,
                                 printHits, silent, lowRes, highRes, bins, correlations, shouldLog);
}


double MtzManager::statisticsWithManager(MtzManager *otherManager,
                                         StatisticsFunction *function, RFactorFunction *rFactorFunction,
                                         RFactorType rFactor, bool printHits, bool silent, double lowRes,
                                         double highRes, int bins,
                                         vector<boost::tuple<double, double, double, int> > *correlations,
                                         bool shouldLog)
{
    vector<double> shells;
    
    if (bins == 0)
        bins = 10;
    
    double maxResolution = 0;
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        if (reflection(i)->getResolution() > maxResolution
            && reflection(i)->getResolution() < 1 / 1.2)
        {
            Reflection *bestReflection = reflection(i);
            
            maxResolution = bestReflection->getResolution();
        }
    }
    
    if (highRes == 0 || highRes < 1 / maxResolution)
        highRes = 1 / maxResolution;
    
    int hits = 0;
    double multiplicity = 0;
    
    double statistic = 0;
    
    
    if (bins > 1 || !silent)
    {
        if (function == StatisticsManager::r_split_by_intensity)
        {
            double interval = (lowRes - highRes) / (double)bins;
            
            for (int i = highRes; i < lowRes; i += interval)
            {
                shells.push_back(i);
            }
            
            shells.push_back(lowRes);
        }
        else
        {
            StatisticsManager::generateResolutionBins(lowRes, highRes, bins,
                                                      &shells);
        }
    }
    
    
    if (rFactor == RFactorNone)
        statistic = function(this, otherManager, !printHits, &hits,
                             &multiplicity, shells[0], shells[shells.size() - 1], shouldLog);
    else
        statistic = rFactorFunction(rFactor, this, &hits, &multiplicity, lowRes,
                                    highRes);
   
    if (bins > 1 || !silent)
    {
        for (int i = 0; i < shells.size() - 1; i++)
        {
            double statistic = 0;
            int hits = 0;
            
            if (rFactor == RFactorNone)
                statistic = function(this, otherManager, 1, &hits,
                                     &multiplicity, shells[i], shells[i + 1], shouldLog);
            else
                statistic = rFactorFunction(rFactor, this, &hits, &multiplicity,
                                            shells[i], shells[i + 1]);
            
            if (!silent)
            {
                std::cout << "	" << shells[i] <<  "	" << shells[i + 1] << "	"
                << statistic <<  "	" << hits
                << std::endl;
            }
            
            if (correlations != NULL)
            {
                boost::tuple<double, double, double, int> result = boost::make_tuple(
                                                                                     shells[i], shells[i + 1], statistic, hits);
                correlations->push_back(result);
            }
        }
    }
    
    if (!silent)
    {
        double statistic = 0;
        
        if (rFactor == RFactorNone)
            statistic = function(this, otherManager, 1,  &hits,
                                 &multiplicity, shells[0], shells[shells.size() - 1], shouldLog);
        else
            statistic = rFactorFunction(rFactor, this, &hits, &multiplicity,
                                        lowRes, highRes);
        
        std::cout << "N: *** Overall ***" << std::endl;
        std::cout << "N: " << lowRes << "\t" << highRes << "\t" << statistic
        << "\t" << hits << "\t" << multiplicity << std::endl;
    }
    
    return statistic;
}

double MtzManager::twinnedIntensity(Reflection *refl)
{
    int invReflId = refl->getReflId(1);
    
    Reflection *holder = NULL;
    this->findReflectionWithId(invReflId, &holder);
    
    if (holder == NULL)
        return nan(" ");
    
    double one = refl->meanIntensity();
    double two = holder->meanIntensity();
    
    return sqrt((one * one + two * two) / 2);
}

void MtzManager::individualDetwinningScales(bool final)
{
    bool twoFobsMinusFc = FileParser::getKey("2_FOBS_MINUS_FC", true);
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        int reflid = reflection(i)->getReflId(0);
        int invReflId = reflection(i)->getReflId(1);
        
        Reflection *refHolder = NULL;
        referenceManager->findReflectionWithId(reflid, &refHolder);
        
        if (reflection(i)->isTwinned())
        {
            if (twoFobsMinusFc == false)
            {
                reflection(i)->setFlagScaled(true);
                continue;
            }

            double fobs = refHolder->meanIntensity();
            double fc = reflection(i)->meanIntensity();
            double _2fofc = 2 * fobs - fc;
            

            if (_2fofc == _2fofc && !final)
                reflection(i)->miller(0)->setScaledIntensity(_2fofc);
            else if (final)
                reflection(i)->miller(0)->setScaledIntensity(fobs);
            
            reflection(i)->setFlagScaled(true);
            
            continue;
        }
        
        Reflection *holder = NULL;
        this->findReflectionWithId(invReflId, &holder);
        
        if (holder == NULL)
            continue;
        
        if (reflection(i)->isFlagScaled() || holder->isFlagScaled())
            continue;
    
        reflection(i)->setFlagScaled(true);
        holder->setFlagScaled(true);
        
        double refIntensity = refHolder->meanIntensity();
        double oneInt = reflection(i)->meanIntensity();
        double twoInt = holder->meanIntensity();
        
        double numerator = sqrt(2) * refIntensity;
        double denominator = sqrt(pow(oneInt, 2) + pow(twoInt, 2));
        
        double scale = numerator / denominator;

        reflection(i)->scale(scale);
        holder->scale(scale);
    }
    
    clearReflectionScalingFlags();
}

void MtzManager::clearReflectionScalingFlags()
{
    for (int i = 0; i < reflectionCount(); i++)
    {
        reflection(i)->setFlagScaled(false);
    }
    
    std::cout << std::endl;
}

void MtzManager::copyOtherAmplitudesFromReference()
{
    for (int i = 0; i < reflectionCount(); i++)
    {
        int reflid = reflection(i)->getReflId();
        Reflection *other = NULL;
        referenceManager->findReflectionWithId(reflid, &other);
        
        if (other)
        {
            for (int j = 0; j < reflection(i)->millerCount(); j++)
            {
                reflection(i)->miller(j)->setFc(other->meanIntensity());
            }
        }
    }
}
