/*
 * Vector.h
 *
 *  Created on: 27 Aug 2014
 *      Author: helenginn
 */

#ifndef VECTOR_H_
#define VECTOR_H_
#include <vector>

#include <boost/tuple/tuple.hpp>
#include "parameters.h"

typedef struct
{
	double h;
	double k;
	double l;
} vec;

vec cross_product_for_vectors(vec vec1, vec vec2);
vec new_vector(double h, double k, double l);
double dot_product_for_vectors(vec vec1, vec vec2);
double length_of_vector(vec vect);
void multiply_vector(vec *vec1, double scale);
double angleBetweenVectors(vec vec1, vec vec2);
vec copy_vector(vec old_vec);
void add_vector_to_vector(vec *vec1, vec vec2);
vec vector_between_vectors(vec vec1, vec vec2);
void take_vector_away_from_vector(vec vec1, vec *vec2);
void scale_vector_to_distance(vec *vec, double new_distance);

double cdf(double x, double mean, double sigma);
double _cdf(double x);
double normal_distribution(double x, double mean, double sigma);
double super_gaussian(double x, double mean, double sigma_0, double exponent);

void regression_line(vector<boost::tuple<double, double, double> > values, double &intercept, double &gradient);
double correlation_between_vectors(vector<double> *vec1,
		vector<double> *vec2, vector<double> *weights, int exclude);
double correlation_between_vectors(vector<double> *vec1,
		vector<double> *vec2, vector<double> *weights);
double correlation_between_vectors(vector<double> *vec1,
		vector<double> *vec2);
double correlation_through_origin(vector<double> *vec1,
		vector<double> *vec2, vector<double> *weights = NULL);
double least_squares_between_vectors(vector<double> *vec1,
		vector<double> *vec2, double slope);
double gradient_between_vectors(vector<double> *vec1,
		vector<double> *vec2);
double minimize_gradient_between_vectors(vector<double> *vec1,
		vector<double> *vec2);
double weighted_mean(vector<double> *means, vector<double> *weights = NULL);
double median(vector<double> *means);
void histogram_gaussian(vector<double> *means, vector<int> *freq, double &mean, double &stdev);
double standard_deviation(vector<double> *values, vector<double> *weights = NULL);
double r_factor_between_vectors(vector<double> *vec1,
		vector<double> *vec2, vector<double> *weights, double scale);

double cartesian_to_distance(double x, double y);
double cartesian_to_angle(double x, double y);

#endif /* VECTOR_H_ */
