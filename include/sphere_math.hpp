/*
 * =====================================================================================
 *
 *       Filename:  sphere_math.hpp
 *
 *    Description:  Sphere mathematics functions and objects
 *
 *        Version:  1.0
 *        Created:  06/02/2015 10:50:07 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */


#define PI 3.14159265
#define SMIDGE 0.001

#ifndef SPHERE_MATH_HPP
#define SPHERE_MATH_HPP


struct vec3 {
	vec3(double x, double y, double z);
	vec3();
	double x, y, z;
	double magnitude();
	vec3 normalize();
	void equals(vec3 vec);
	vec3 scalar_multiply(double b);
	vec3 operator+(vec3 b);
	vec3 operator-(vec3 b);
};


double rand_range(double min, double max);
double dot_product(double v[], double u[], int n);
double dot_product(vec3 v, vec3 u);
double distance(vec3 p1, vec3 p2);
double intersect_correct(int i, int j);
vec3 axis_rotation(vec3 pos, vec3 axis, double theta);
vec3 cross_product(vec3 i, vec3 j);

#endif
