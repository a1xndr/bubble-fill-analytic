/*
 * =====================================================================================
 *
 *       Filename:  sphere_math.cpp
 *
 *    Description:  Functions for working with spheres in 3d including collision checks,
 *    distance between two point calculations,etc
 *
 *        Version:  1.0
 *        Created:  06/02/2015 10:29:27 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexander Oleinik 
 *   Organization:  
 *
 * =====================================================================================
 */

#include "sphere_math.hpp"
#include <stdlib.h>
#include <math.h>
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  rand_range
 *  Description:  Choose a random double bewtween min and max. Make sure srand is called 
 *  with a seed before using this function.
 * =====================================================================================
 */
double rand_range(double min, double max)
{
    return min + static_cast<double>(rand()/static_cast<double>(RAND_MAX) * (max - min));
}

double poisson(double min, double max)
{
    return min + static_cast<double>(rand()/static_cast<double>(RAND_MAX) * (max - min));
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  dot_product
 *  Description:  Dot product of the first n indices of two vectors
 * =====================================================================================
 */
double dot_product(double v[], double u[], int n)
{
    double result = 0.0;
    for (int i = 0; i < n; i++)
        result += v[i]*u[i];
    return result;
}

double dot_product(vec3 v, vec3 u)
{
    return v.x*u.x + v.y*u.y + v.z*u.z;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  distance
 *  Description:  Calculate the distance between the centers of two spheres
 *
 * =====================================================================================
 */
double distance(vec3 p1, vec3 p2)
{
		double distance = pow(	pow(p2.x - p1.x , 2)
				+pow(p2.y - p1.y , 2)
				+pow(p2.z - p1.z , 2),0.5);
		return distance;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  axis_rotation
 *  Description:  Rotates a position about an arbitrary axis for and angle theta
 *  Thoroughly described in http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
 *  Accepts non-normalized axes
 *
 * =====================================================================================
 */
vec3 axis_rotation(vec3 pos, vec3 axis, double theta)
{
	vec3 tform_position;

	double x = pos.x;
	double y = pos.y;
	double z = pos.z;

	double u = axis.x/axis.magnitude();
	double v = axis.y/axis.magnitude();
	double w = axis.z/axis.magnitude();

    tform_position.x = 	u*(u*x+v*y+w*z)*(1-cos(theta)) +
    					x*cos(theta) + (-w*y +v*z)*sin(theta);

    tform_position.y = 	v*(u*x+v*y+w*z)*(1-cos(theta)) +
    					y*cos(theta) + (w*x -u*z)*sin(theta);

    tform_position.z = 	w*(u*x+v*y+w*z)*(1-cos(theta)) +
    					z*cos(theta) + (-v*x +u*y)*sin(theta);

    return tform_position;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  axis_rotation
 *  Description:  Rotates a position about an arbitrary axis for and angle theta
 *  Thoroughly described in http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
 *  Accepts non-normalized axes
 *
 * =====================================================================================
 */
vec3 cross_product(vec3 i, vec3 j)
{
	vec3 k(	i.y*j.z-i.z*j.y,
			i.z*j.x-i.x*j.z,
			i.x*j.y-i.y*j.x);
    return k;
}

vec3 scalar_product(vec3 i, double j)
{
	vec3 k(	i.x*j,
	        i.y*j,
		i.z*j);
    return k;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  equals
 *  Description:  Constructor with component intialization
 * =====================================================================================
 */
vec3::vec3(double X, double Y, double Z)
{
	x=X;
	y=Y;
	z=Z;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  equals
 *  Description:  Default Constructor
 * =====================================================================================
 */
vec3::vec3()
{
	x=0;
	y=0;
	z=0;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  magnitude
 *  Description:  Returns magnitude of the vector
 * =====================================================================================
 */
double vec3::magnitude()
{
	return distance(vec3(0,0,0),vec3(x,y,z));
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  normalize
 *  Description:  Normalizes a vector and returns the result as a new vector
 * =====================================================================================
 */
vec3 vec3::normalize()
{
	vec3 norm(  this->x/this->magnitude(),
		    this->y/this->magnitude(),
		    this->z/this->magnitude());
	return norm;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  equals
 *  Description:  Set the vector equal to another
 * =====================================================================================
 */
void vec3::equals(vec3 vec){
	x = vec.x;
	y = vec.y;
	z = vec.z;
	return ;
}

vec3 vec3::operator+(vec3& b){
    vec3 sum(	this->x+b.x,
		this->y+b.y,
		this->z+b.z);
    return sum;
}

vec3 vec3::operator-(vec3& b){
    vec3 diff(	this->x-b.x,
		this->y-b.y,
		this->z-b.z);
    return diff;
}
vec3 vec3::operator*(double b){
    vec3 mult(	this->x*b,
		this->y*b,
		this->z*b);
    return mult;
}
vec3 vec3::scalar_multiply(double b){
    vec3 scalarmult(this->x*b,
		this->y*b,
		this->z*b);
    return scalarmult;
}
//double intersect_correct(int i, int j)
//{
//    //The target distnance between the two spheres: i.e. the sum of the radii
//    double target_distance = spheres[i].radius + spheres[j].radius;
//    //The length we need to correct for
//    double error = target_distance-distance(i,j);
//    //The vector
//    //double;
//    /* :TODO:06/01/2015 02:08:36 PM::  Use Linear Algebra Lib to implement vector
//     * subtraction*/
//    return 0;
//}
