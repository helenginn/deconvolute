#ifndef mtz_manager
#define mtz_manager

#include <string>
#include <iostream>
#include "cmtzlib.h"
#include "csymlib.h"
#include <vector>

//#include <tr1/memory>
#include "definitions.h"
#include "parameters.h"
#include <tuple>
#include <iomanip>

using namespace CMtz;
using namespace CSym;

class Reflection;
class Miller;
class Matrix;

typedef enum
{
	ScoreTypeMinimizeRSplit = 0,
	ScoreTypeCorrelation = 1,
	ScoreTypePartialityCorrelation = 2,
	ScoreTypeMinimizeRSplitLog = 3,
	ScoreTypeCorrelationLog = 4,
	ScoreTypeAngle = 5,
	ScoreTypePartialityLeastSquares = 6,
	ScoreTypePartialityGradient = 7,
    ScoreTypeSymmetry = 8,
    ScoreTypeStandardDeviation = 9,
    ScoreTypeMinimizeRMeas = 10,
} ScoreType;

typedef enum
{
	TrustLevelGood, TrustLevelAverage, TrustLevelBad
} TrustLevel;

class Miller;

class MtzManager
{

protected:
    std::string filename;
	CCP4SPG *low_group;
	static bool reflection_comparison(Reflection *i, Reflection *j);

	double extreme_index(MTZ *mtz, int max);
	void hkls_for_reflection(MTZ *mtz, float *adata, int *h, int *k, int *l,
			int *multiplier, int *offset);
	int index_for_reflection(int h, int k, int l, bool inverted);
	void findMultiplier(MTZ *mtz, int *multiplier, int *offset);

	// minimisation stuff

    int requiredPriority;
	vector<Reflection *> reflections;
	vector<Reflection *> refReflections;
	vector<Reflection *> matchReflections;
    MtzManager *previousReference;
    int previousAmbiguity;
    int activeAmbiguity;
    
    double scale;
	bool finalised;
	bool inverse;
	bool flipped;
    int failedCount;
	static double lowRes;
	static double highRes;
	bool freePass;
	bool rejected;

	TrustLevel trust;
	ScoreType scoreType;

    MatrixPtr baseMatrix;
	static MtzManager *referenceManager;
	MtzManager *lastReference;
    
    double cellDim[3];
    double cellAngles[3];

    static double unitCellScore(void *object);
    double wavelengthStandardDeviation();
    std::ostringstream logged;
public:
    vector<double> superGaussianTable;
    double bFactor;
    
    MtzManager(std::string nameOfFile = "");
	virtual ~MtzManager(void);
	MtzPtr copy();
	void loadParametersMap();

    void addReflections(vector<Reflection *>reflections);
	void clearReflections();
	void addReflection(Reflection *reflection);
	void removeReflection(int i);
	
    MatrixPtr matrix;
    bool checkUnitCell(double trueA, double trueB, double trueC, double tolerance);
    
	void setFilename(std::string name);
    std::string getFilename(void);
	void description(void);
    
    void setDefaultMatrix();
	void setMatrix(double *components);
	void setMatrix(MatrixPtr newMat);
	void insertionSortReflections(void);
	void sortLastReflection(void);
    void incrementActiveAmbiguity();
    double maxResolution();

    void copyOtherAmplitudesFromReference();
    std::string filenameRoot();
	void setSpaceGroup(int spgnum);
    virtual void loadReflections(bool calculated);
	static void setReference(MtzManager *reference);
	int findReflectionWithId(int refl_id, Reflection **reflection, bool insertionPoint = false);
	void findCommonReflections(MtzManager *other,
			vector<Reflection *> &reflectionVector1, vector<Reflection *> &reflectionVector2,
			int *num = NULL, bool singlets = false);
	double meanCorrectedIntensity(Reflection *reflection);
	double gradientAgainstManager(MtzManager &otherManager, bool leastSquares = false, double lowRes = 0, double highRes = 0);
	void applyBFactor(double bFactor);
	void applyScaleFactor(double scaleFactor, double lowRes = 0, double highRes = 0, bool absolute = false);
	void applyScaleFactorsForBins(int bins);
	void clearScaleFactor();
	void makeScalesPermanent();
	double averageIntensity(void);
	void getRefReflections(vector<Reflection *> *refPointer,
			vector<Reflection *> *matchPointer);
	
	void setUnitCell(double a, double b, double c, double alpha, double beta,
			double gamma);
	void setUnitCell(vector<double> unitCell);
	void getUnitCell(double *a, double *b, double *c, double *alpha, double *beta,
			double *gamma);
	void copySymmetryInformationFromManager(MtzPtr toCopy);
	
	void writeToFile(std::string newFilename, bool announce = false, bool shifts = false, bool includeAmbiguity = false);
	void writeToDat();
    void sendLog(LogLevel priority = LogLevelNormal);

    double rSplitIntensityBinsWithManager(MtzManager *otherManager, bool printHits,
                                          bool silent, double lowInt, double highInt, int bins,
                                          vector<boost::tuple<double, double, double, int> > *correlations,
                                          bool shouldLog);
    
	double correlationWithManager(MtzManager *otherManager, bool printHits = false,
			bool silent = true, double lowRes = 0, double highRes = 0, int bins = 20,
			vector<boost::tuple<double, double, double, int> > *correlations = NULL, bool shouldLog = false);

	double rSplitWithManager(MtzManager *otherManager, bool printHits = false,
			bool silent = true, double lowRes = 0, double highRes = 0, int bins = 20,
			vector<boost::tuple<double, double, double, int> > *correlations = NULL, bool shouldLog = false);

	double statisticsWithManager(MtzManager *otherManager, StatisticsFunction *function,
			bool printHits, bool silent, double lowRes, double highRes, int bins,
			vector<boost::tuple<double, double, double, int> > *correlations,
			bool shouldLog);

	double statisticsWithManager(MtzManager *otherManager,
			StatisticsFunction *function, RFactorFunction *rFactorFunction,
			RFactorType rFactor, bool printHits, bool silent, double lowRes,
			double highRes, int bins,
			vector<boost::tuple<double, double, double, int> > *correlations,
			bool shouldLog);

	double correlation(bool silent = true, double lowResolution = 0, double highResolution = -1);
	double rSplit(double low, double high, bool square = false);
    
// detwinning stuff
    
    double twinnedIntensity(Reflection *refl);
    void individualDetwinningScales(bool final = false);
    void clearReflectionScalingFlags();
    
    int ambiguityCount();
    
    void flipToActiveAmbiguity();
    void resetFlip();
    void setAdditionalWeight(double weight);
    
    int getActiveAmbiguity()
    {
        return activeAmbiguity;
    }
    
    void setActiveAmbiguity(int newAmbiguity);

    virtual int reflectionCount()
    {
        return (int)reflections.size();
    }
    
    virtual Reflection *reflection(int i)
    {
        return reflections[i];
    }
    
    ScoreType getScoreType()
	{
		return scoreType;
	}

	void setScoreType(ScoreType newScoreType)
	{
		scoreType = newScoreType;
	}

	static double getHighRes()
	{
		return highRes;
	}

	static void setHighRes(double _highRes)
	{
		highRes = _highRes;
	}

	static double getLowRes()
	{
		return lowRes;
	}

	static void setLowRes(double _lowRes)
	{
		lowRes = _lowRes;
	}

	CCP4SPG*& getLowGroup()
	{
		return low_group;
	}

	void setLowGroup(CCP4SPG*& lowGroup)
	{
		low_group = lowGroup;
	}

	TrustLevel getTrust() const
	{
		return trust;
	}

	void setTrust(TrustLevel trust)
	{
		this->trust = trust;
	}

	bool isFinalised() const
	{
		return finalised;
	}

	void setFinalised(bool finalised)
	{
		this->finalised = finalised;
	}

	bool isFreePass() const
	{
		return freePass;
	}

	void setFreePass(bool freePass)
	{
		this->freePass = freePass;
	}

	bool isFlipped() const
	{
		return flipped;
	}

	MatrixPtr getMatrix()
	{
		return matrix;
	}

	static MtzManager*& getReferenceManager()
	{
        return referenceManager;
	}
};

#endif

