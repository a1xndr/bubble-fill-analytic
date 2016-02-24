/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  There are 3 States in this process.
 *    ------------------ 1 ---------------------------
 *    In the first state, the bubble is free to grow and move. IE 
 *
 *        Version:  1.0
 *        Created:  01/27/2016 10:52:28 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexander Oleinik (ax), alxndr@bu.edu
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream> 
#include <string>
#include <cmath>
#include <Eigen/LU>
#include "bubble.hpp"
#include "sphere.hpp"
#include "sphere_math.hpp"

const int spheres_max = 1000000;
const int bubbles_max = 1000000;
const int x_max=10;
const int y_max=10;
const int z_max=10;
double sphere_volume=0;

sphere spheres[spheres_max];
bubble bubbles[bubbles_max];


int check_intersect(double radius, vec3 pos, double i, double k, bool fastx, bool fasty, bool fastz){
    for(int j=i-1; j >=0; j--)
    {
    	if(fastx && fabs(pos.x-spheres[j].pos.x) > radius + spheres[j].radius)continue;
    	if(fasty && fabs(pos.y-spheres[j].pos.y) > radius + spheres[j].radius)continue;
    	if(fastz && fabs(pos.z-spheres[j].pos.z) > radius + spheres[j].radius)continue;
    	//Square of distance between two centers comparison against square of radius.
    	if( distance(pos,spheres[j].pos)
	    < spheres[j].radius + radius - SMIDGE)return j;
    }
    for(int j=k-1; j >=0; j--)
    {
    	if(fastx && fabs(pos.x-bubbles[j].pos.x) > radius + bubbles[j].radius)continue;
    	if(fasty && fabs(pos.y-bubbles[j].pos.y) > radius + bubbles[j].radius)continue;
    	if(fastz && fabs(pos.z-bubbles[j].pos.z) > radius + bubbles[j].radius)continue;
    	//Square of distance between two centers comparison against square of radius.
    	if( distance(pos,bubbles[j].pos)
	    < bubbles[j].radius + radius - SMIDGE)return j;
    }
    return -1;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  stage1(bubble *s, int i)
 *  Description:  Grows a sphere not in contact with any other bubble, until
 *  a collision
 * =====================================================================================
 */
sphere stage1(bubble *b, int num_spheres, int num_bubbles)
{   
    /*-----------------------------------------------------------------------------
     *  First we need to find the sphere that we will be "growing" up to
     *-----------------------------------------------------------------------------*/
    double min_distance = x_max;
    int index = -1;
    bool bubble = false;
    for(int i=0; i<num_spheres; i++)
    {
        double d = distance(b->pos, spheres[i].pos) 
            -b->radius - spheres[i].radius; 
        if( d < min_distance )
        {
            index = i;
            min_distance = d;
        }
    }
    for(int i=0; i<num_bubbles; i++)
    {
        double d = distance(b->pos, bubbles[i].pos) 
            -b->radius - bubbles[i].radius; 
        if( d < min_distance )
        {
            bubble=true;
            index = i;
            min_distance = d;
        }
    }

    //The vector from the collision sphere to the active sphere
    vec3 v1; 
    if(bubble)
    {
        v1 = bubbles[index].pos - b->pos; 
        b->radius = v1.magnitude()- bubbles[index].pos.magnitude();
    } else {
        v1 = spheres[index].pos - b->pos; 
        b->radius = v1.magnitude()- spheres[index].pos.magnitude(); 
    }
    if(index==-1)return bubbles[0];
    if(bubble)  return bubbles[index];
    else        return spheres[index];
}

sphere stage2(bubble * s0, sphere * s1, int num_spheres, int num_bubbles)
{
    double min_radius;
    int index = -1;
    bool bubble = false;
    vec3 t = (s0->pos - s1->pos).normalize();
    vec3 c1 = s1->pos+scalar_product(t,s1->radius);
    vec3 vc1 = c1  - s1->pos; 
    for(int i; i <num_spheres+num_bubbles; i++)
    {
        sphere sc;  //candidate for s2
        if(i>num_spheres) sc = bubbles[i-num_spheres];
        else sc = spheres[i];

        vec3 v1s = s1->pos - sc.pos;
        double r0c = ( pow(s1->radius, 2) - pow(sc.radius, 2) 
            + pow(v1s.magnitude(), 2) + 2*dot_product(vc1, v1s)      ) / 
            (2*(s1->radius + sc.radius - dot_product(vc1, v1s)));
        if(r0c < min_radius)
        {
            min_radius = r0c;
            if(i>num_spheres)
            {
                bubble = true;
                index = i-num_spheres;
            }
            else index =i;
        }
    }
    if(index==-1)return bubbles[0];
    if(bubble)return bubbles[index];
    else return spheres[index];
}


sphere stage3(bubble * s0, sphere * s1, sphere * s2, int num_spheres, int num_bubbles)
{
    double min_radius;
    int index = -1;
    bool bubble = false;

    vec3 vc1, b1, b2;
    double r1 = s1->radius;
    double r2 = s2->radius;
//    double * rs = s0->radius;

    vec3 x2 = s2->pos;
    vec3 x1 = s1->pos;   
    
    vec3 c1 = s1->pos+((s0->pos - s1->pos).normalize()).scalar_multiply(s1->radius);
    vc1 = c1  - s1->pos; 
    
    b1 = s0->pos    - x1; 
    b2 = x2         - x1;

    vec3 vc1xb2 = cross_product(vc1, b2);
    Eigen::Matrix4f A;
    Eigen::Vector4f b;
    for(int i=0; i <num_spheres+num_bubbles; i++)
    {
        sphere sc;  //candidate for s2
        if(i>num_spheres) sc = bubbles[i-num_spheres];
        else sc = spheres[i];
        double rs = sc.radius;
        vec3 bs = sc.pos - x1;
        A   << (b2-bs).x, (b2-bs).y, (b2-bs).z, (r2-rs)
            , b2.x         , b2.y         , b2.z         , (r2-r1)
            , bs.x         , bs.y         , bs.z         , (rs-r1)
            , vc1xb2.x     , vc1xb2.y     , vc1xb2.z     , 0      ;

        b   << (pow(b2.magnitude(), 2) - pow(bs.magnitude(), 2) 
                    + pow(rs, 2) - pow(r2,2))/2.0
            , (pow(b2.magnitude(), 2) + pow(r1, 2) - pow(r2,2))/2.0
            , (pow(b1.magnitude(), 2) + pow(r1, 2) - pow(rs,2))/2.0
            , 0 ;
        Eigen::Vector4f soln = A.inverse() * b;

        if(soln[3] < min_radius)
        {
            min_radius = soln[3];
            if(i>num_spheres)
            {
                bubble = true;
                index = i-num_spheres;
            }
            else index =i;
        }
    }
    if(index==-1)return bubbles[0];
    if(bubble)return bubbles[index];
    else return spheres[index];
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  read_sphere_file
 *  Description:  Reads sphere coordinates from file
 * =====================================================================================
 */
int read_sphere_coords(std::string path){
    std::ifstream file(path);
    int count = 0;
//    double pore_volume = x_max*y_max+z_max - sphere_volume;
    if(file.is_open())
    {
        std::string line;
        while(getline(file, line))
        {
            std::istringstream iss(line);
            double d;
            double params[5];
            int paramcount=0;;
            while (iss>>d){
                params[paramcount]=d;
                paramcount++;
            }
            spheres[count]=(sphere){params[1], (vec3){params[2],params[3],params[4]}};
            sphere_volume+=(4.0/3.0)*PI*pow(spheres[count].radius,3);
            count++;
        }
    }
    return count;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  main loop
 * =====================================================================================
 */
int main(int argc, char * argv[])
{
    if(argc!=3)
    {
        std::cout << 
            "Usage: ./bubble-fill [input-sphere-file] [output-bubble-file]" <<
            std::endl;
        return -1;
    }
    std::string sphere_file_path = argv[1];
    std::string bubble_file_path = argv[2];
    int num_spheres = read_sphere_coords(sphere_file_path);
    for(int i=0; i<1000000; i++)
    {
        bubble b;
        double r, x, y, z;
        while(true)
        {
            double r= 0.001;
            double x = rand_range(0.001,9.999);
            double y = rand_range(0.001,9.999);
            double z = rand_range(0.001,9.999);
            if(check_intersect(r, (vec3){x, y, z}, num_spheres, i, true, true, true)
                    ==-1) break;
        }
        b.radius = r;
        b.pos.x = x;
        b.pos.y = y;
        b.pos.z = z;
        sphere c1 = stage1(&b, num_spheres, i);
        sphere c2 = stage2(&b, &c1, num_spheres, i);
        sphere c3 = stage3(&b, &c1, &c2, num_spheres, i);
        bubbles[i]=b;
        std::cout<< b.radius << " " << 
            b.pos.x << " " << b.pos.y << " " << b.pos.z << std::endl;
    }
}
