#include <string>
#include <iostream>
#include <cmath>
#include "misc.h"
#include "StatisticsManager.h"
#include <new>
#include <fstream>
#include "FileParser.h"
#include "Holder.h"


void StatisticsManager::setMtzs(vector<MtzPtr> newMtzs)
{
	mtzs = newMtzs;
    mtz_num = (int)newMtzs.size();
}

void StatisticsManager::write_refls(int num)
{
	MtzManager *shot = &*(mtzs[num]);

	for (int i = 0; i < (int) (*shot).reflectionCount(); i++)
	{
		if (shot->reflection(i)->millerCount() < 10)
			continue;

		double mean_intensity = shot->reflection(i)->meanIntensity();
		for (int j = 0; j < shot->reflection(i)->millerCount(); j++)
		{
			double percentage = shot->reflection(i)->miller(j)->intensity()
					/ mean_intensity * 100;
			std::cout << percentage << std::endl;
		}

		std::cout << std::endl;
	}
}

double StatisticsManager::midPointBetweenResolutions(double minD, double maxD)
{
	vector<double> bins;
	generateResolutionBins(minD, maxD, 2, &bins);

	return bins[1];
}

void StatisticsManager::generateResolutionBins(double minD, double maxD,
		int binCount, vector<double> *bins)
{
	double minRadius = (minD == 0) ? 0 : 1 / minD;
	double maxRadius = 1 / maxD;

	if (maxD <= 0)
	{
		std::cout << "Warning: maximum resolution set to 0. Ignoring.";
        return;
	}

	double maxVolume = pow(maxRadius, 3);
	double minVolume = pow(minRadius, 3);
	double totalVolume = maxVolume - minVolume;

	double eachVolume = totalVolume / binCount;

	double r1 = minRadius;
	double r2 = 0;

	bins->push_back(1 / r1);

	for (int i = 0; i < binCount; i++)
	{
		double r2_cubed = pow(r1, 3) + eachVolume;
		r2 = cbrt(r2_cubed);

		bins->push_back(1 / r2);

		r1 = r2;
	}
}

double StatisticsManager::cc_pearson(int num1, int num2, int silent,
		int *hits, double *multiplicity, double lowResolution,
		double highResolution, bool shouldLog)
{
	MtzManager *shot1 = &*(mtzs[num1]);
	MtzManager *shot2 = &*(mtzs[num2]);

	return cc_pearson(shot1, shot2, silent, hits, multiplicity,
			lowResolution, highResolution, shouldLog);
}

void StatisticsManager::convertResolutions(double lowAngstroms,
		double highAngstroms, double *lowReciprocal, double *highReciprocal)
{
	*lowReciprocal = (lowAngstroms == 0) ? 0 : 1 / lowAngstroms;
	*highReciprocal = (highAngstroms == 0) ? FLT_MAX : 1 / highAngstroms;
}

double StatisticsManager::cc_pearson(MtzManager *shot1, MtzManager *shot2,
		int silent, int *hits, double *multiplicity,
		double lowResolution, double highResolution, bool twinnedOnly)
{
	vector<Reflection *> reflections1;
	vector<Reflection *> reflections2;
	int num = 0;

	shot1->findCommonReflections(shot2, reflections1, reflections2, &num, twinnedOnly);

    if (reflections1.size() <= 2)
	{
		return -1;
	}

	double invHigh = 1 / highResolution;
	double invLow = 1 / lowResolution;

	if (highResolution == 0)
		invHigh = FLT_MAX;

	if (lowResolution == 0)
		invLow = 0;

	if (!silent)
	{
		for (int i = 0; i < num; i++)
		{
            int h = reflections1[i]->miller(0)->getH();
            int k = reflections1[i]->miller(0)->getK();
            int l = reflections1[i]->miller(0)->getL();
            
			double int1 = reflections1[i]->meanIntensity();
			std::string filename = shot1->getFilename();
			double int2 = reflections2[i]->meanIntensityWithExclusion(&filename);

			if (int1 != int1 || int2 != int2)
				continue;

			if (!(reflections1[i]->betweenResolutions(lowResolution, highResolution)))
				continue;
            
            double resolution = 1 / reflections1[i]->getResolution();

			std::cout << h << "\t" << k << "\t" << l << "\t" << int1 << "\t" << int2 << "\t" << resolution << std::endl;
		}
	}

	double sum_x = 0;
	double sum_y = 0;
	double weight_counted = 0;
	int num_counted = 0;
    
    for (int i = 0; i < num; i++)
	{
		if (!(reflections1[i]->betweenResolutions(lowResolution, highResolution)))
            continue;
        
        double weight = 1;
        
		double mean1 = reflections1[i]->meanIntensity();
		double mean2 = reflections2[i]->meanIntensity();

		if (mean1 != mean1 || mean2 != mean2)
			continue;

		num_counted++;
		sum_x += mean1 * weight;
		sum_y += mean2 * weight;

		weight_counted += weight;
	}

	double mean_x = sum_x / weight_counted;
	double mean_y = sum_y / weight_counted;

	double sum_x_y_minus_mean_x_y = 0;
	double sum_x_minus_mean_x_sq = 0;
	double sum_y_minus_mean_y_sq = 0;

	for (int i = 0; i < num; i++)
	{
		if (!(reflections1[i]->getResolution() > invLow
				&& reflections1[i]->getResolution() < invHigh))
			continue;

		double amp_x = reflections1[i]->meanIntensity();
		std::string filename = shot1->getFilename();
		double amp_y = reflections2[i]->meanIntensityWithExclusion(&filename);

        double weight = 1;

        if (weight < 0)
            continue;

		if (amp_x != amp_x || amp_y != amp_y)
			continue;

		if (weight != weight)
			continue;

		sum_x_y_minus_mean_x_y += weight * (amp_x - mean_x) * (amp_y - mean_y);
		sum_x_minus_mean_x_sq += weight * pow(amp_x - mean_x, 2);
		sum_y_minus_mean_y_sq += weight * pow(amp_y - mean_y, 2);
	}

	sum_x_y_minus_mean_x_y /= weight_counted;
	sum_x_minus_mean_x_sq /= weight_counted;
	sum_y_minus_mean_y_sq /= weight_counted;

	if (hits != NULL)
		*hits = num_counted;

	double r = sum_x_y_minus_mean_x_y
			/ (sqrt(sum_x_minus_mean_x_sq) * sqrt(sum_y_minus_mean_y_sq));

	if (r < 0)
		r = 0;
	if (r != r)
		r = -1;

	return r;
}

double StatisticsManager::r_split_by_intensity(MtzManager *shot1, MtzManager *shot2,
                                  int silent, int *hits, double *multiplicity,
                                  double lowIntensity, double highIntensity, bool twinnedOnly)
{
    double sum_numerator = 0;
    double sum_denominator = 0;
    double mult = 0;
    int count = 0;
    
    double dMin = 0;
    double dMax = 0;
    
    std::string filename = shot1->getFilename();
    
    for (int i = 0; i < shot1->reflectionCount(); i++)
    {
        Reflection *reflection = shot1->reflection(i);
        
        if (twinnedOnly && !reflection->isTwinned())
        {
            continue;
        }
        
        int reflid = reflection->getReflId();
        
        Reflection *reflection2 = NULL;
        shot2->findReflectionWithId(reflid, &reflection2);
        
        if (reflection2 == NULL)
            continue;
        
        double int1 = reflection->meanIntensity();
        double int2 = reflection2->meanIntensity();
        
        if (int1 < lowIntensity || int1 > highIntensity)
            continue;
        
        if (int1 != int1 || int2 != int2)
            continue;
        
        double res = reflection->getResolution();

        if (!silent)
        {
            std::cout << int1 << "\t" << int2 << "\t" << res << std::endl;
        }

        
        count++;
        sum_numerator += fabs(int1 - int2);
        sum_denominator += fabs((int1 + int2) / 2);
        
        
        
//        mult += reflection->acceptedCount() + reflection2->acceptedCount();
    }
    
    mult /= count;
    
    if (hits != NULL)
        *hits = count;
    
    if (multiplicity != NULL)
        *multiplicity = mult;
    
    double r_split = sum_numerator / (sum_denominator);
    
    return r_split;
}

double StatisticsManager::r_split(MtzManager *shot1, MtzManager *shot2,
		int silent, int *hits, double *multiplicity,
		double lowResolution, double highResolution, bool twinnedOnly)
{
	double sum_numerator = 0;
	double sum_denominator = 0;
	double mult = 0;
	int count = 0;

	double dMin = 0;
	double dMax = 0;

	convertResolutions(lowResolution, highResolution, &dMin, &dMax);
    std::string filename = shot1->getFilename();

	for (int i = 0; i < shot1->reflectionCount(); i++)
	{
        Reflection *reflection = shot1->reflection(i);
        
        if (twinnedOnly && !reflection->isTwinned())
        {
            continue;
        }
        
        int reflid = reflection->getReflId();

		Reflection *reflection2 = NULL;
		shot2->findReflectionWithId(reflid, &reflection2);

		if (reflection2 == NULL)
			continue;
        
		double int1 = reflection->meanIntensity();
		double int2 = reflection2->meanIntensity();

		if (int1 != int1 || int2 != int2)
			continue;

		double res = reflection->getResolution();

		if (res < dMin || res > dMax)
			continue;
        
		count++;
		sum_numerator += fabs(int1 - int2);
		sum_denominator += fabs((int1 + int2) / 2);

		mult += reflection->acceptedCount() + reflection2->acceptedCount();
	}

	mult /= count;

	if (hits != NULL)
		*hits = count;

	if (multiplicity != NULL)
		*multiplicity = mult;

	double r_split = sum_numerator / (sum_denominator);

	return r_split;
}

StatisticsManager::StatisticsManager(void)
{
	mtz_num = 0;
	mtzs.clear();
}

StatisticsManager::~StatisticsManager(void)
{

}

