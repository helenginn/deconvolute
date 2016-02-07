#ifndef statistics
#define statistics


#include <string>
#include <iostream>
#include "MtzManager.h"

struct Partial
{
	double partiality;
	double percentage;
	double resolution;
	double wavelength;
};

class StatisticsManager
{
private:

public:
	StatisticsManager(void);
	~StatisticsManager(void);

    vector<MtzPtr> mtzs;
	
	void printGradientsAgainstRef(MtzManager *reference);

	void setMtzs(vector<MtzPtr> mtzs);

	void write_refls(int num);

	static double midPointBetweenResolutions(double minD, double maxD);
	static void generateResolutionBins(double minD, double maxD, int binCount,
			vector<double> *bins);
	static void convertResolutions(double lowAngstroms, double highAngstroms,
			double *lowReciprocal, double *highReciprocal);

	double cc_through_origin(int num1, int num2, int silent, int inverted,
			int *hits);
	double cc_through_origin(int num1, int num2, int silent, int inverted,
			int *hits, double lowResolution, double highResolution, bool log);

	static double cc_pearson(MtzManager *shot1, MtzManager *shot2, int silent,
            int *hits, double *multiplicity, double lowResolution,
			double highResolution, bool log);
	double cc_pearson(int num1, int num2, int silent, int *hits,
			double *multiplicity, double lowResolution = 0, double highResolution =
					0, bool log = false);

	static double r_split(MtzManager *shot1, MtzManager *shot2, int silent,
			int *hits, double *multiplicity, double lowResolution,
			double highResolution, bool log);

    static double r_split_by_intensity(MtzManager *shot1, MtzManager *shot2,
                                             int silent, int *hits, double *multiplicity,
                                             double lowIntensity, double highIntensity, bool twinnedOnly);
    
	vector<vector<double> > gradient_array;

	int mtz_num;
};

#endif // statistics
