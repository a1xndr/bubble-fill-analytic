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
#include <string>
#include <cmath>
#include "bubble.hpp"
#include "sphere.hpp"
#include "sphere_math.hpp"

const int spheres_max = 1000000;
const int bubbles_max = 1000000;
const int x_max=10;
const int y_max=10;
const int z_max=10;

sphere spheres[spheres_max];
bubble bubbles[bubbles_max];

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
    int index;
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
    if(bubble)  return bubbles[index];
    else        return spheres[index];
}

sphere stage2(bubble * s0, sphere * s1, int num_spheres, int num_bubbles)
{
    double min_radius;
    int index;
    bool bubble;

    vec3 c1  = s1->pos +(s0->pos - s1->pos).normalize() * s1->radius;
    vec3 vc1 = c1  - s1->pos; 
    for(int i; i <num_spheres+num_bubbles; i++)
    {
        sphere sc;  //candidate for s2
        if(i>num_spheres) sc = bubbles[i-num_spheres];
        else sc = spheres[i];

        vec3 v1s = s1->pos - sc->pos;
        double r0c = ( pow(s1->radius, 2) - pow(sc.radius, 2) 
            + pow(v1s.magnitude(), 2) + 2*dot_product(vc1, v1s)      ) / 
            (2*(s1->radius + sc->radius - dot_product(vc1, v1s)));
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
    if(bubble)return bubbles[i];
    else return spheres[i];
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
    double pore_volume = X_MAX+Y_MAX+Z_MAX - sphere_volume;
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
}
