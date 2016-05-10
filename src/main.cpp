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

#define DEBUG(y, x) do { if(y<=debug_level){for(int dinc=0; dinc<y; dinc++){std::cout<< "     ";} std::cout << x;} } while (0)

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream> 
#include <vector> 
#include <string>
#include <cmath>
#include <Eigen/LU>
#include <Eigen/QR>
#include "bubble.hpp"
#include "sphere.hpp"
#include "sphere_math.hpp"

const int spheres_max = 1000000;
const int bubbles_max = 1000000;
const int x_max=10;
const int y_max=10;
const int z_max=10;
double sphere_volume=0;

int debug_level = 3;

sphere spheres[spheres_max];
bubble bubbles[bubbles_max];
std::vector<sphere> p_spheres; //pseudospheres



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
 *  
 *  a collision
 * =====================================================================================
 */
sphere stage1(bubble *s0, int num_spheres, int num_bubbles)
{   
    /*-----------------------------------------------------------------------------
     *  First we need to find the sphere that we will be "growing" up to
     *-----------------------------------------------------------------------------*/
    double min_distance = x_max;
    double min_radius = x_max;
    int index = -1;
    bool bubble = false;
    for(int i=0; i<num_spheres+num_bubbles+2; i++)
    {
        sphere sc;  //candidate for s2
        if(i>num_spheres) sc = bubbles[i-num_spheres];
        else sc = spheres[i];
        double d = distance(s0->pos, sc.pos) 
            -s0->radius - sc.radius; 
        if( d < min_distance )
        {
            index = i;
            if(i>num_spheres)
            {
                bubble = true;
                index = i-num_spheres;
            }
            min_distance = d;
        }
    }

    vec3 min_vec = (vec3){-s0->pos.x, 0, 0};
    min_vec = s0->pos.y > min_vec.magnitude() ? (vec3){0, -s0->pos.y, 0}: min_vec; 
    min_vec = s0->pos.z > min_vec.magnitude() ? (vec3){0, 0, -s0->pos.z}: min_vec; 
    min_vec = x_max - s0->pos.x > min_vec.magnitude() ? (vec3){x_max - s0->pos.x, 0, 0}: min_vec ; 
    min_vec = y_max - s0->pos.y > min_vec.magnitude() ? (vec3){0, y_max - s0->pos.y, 0}: min_vec ; 
    min_vec = z_max - s0->pos.z > min_vec.magnitude() ? (vec3){0, 0, z_max - s0->pos.z}: min_vec ; 
    if(min_vec.magnitude()<min_distance)
    {
        min_radius=min_vec.magnitude();
        s0->radius = min_distance;
        p_spheres.push_back(sphere{0, s0->pos-min_vec});
        return p_spheres.back();
    }
    //The vector from the collision sphere to the active sphere
    vec3 v1; 
    if(bubble)
    {
        v1 = bubbles[index].pos - s0->pos; 
        s0->radius = v1.magnitude()-bubbles[index].radius;
    } else {
        v1 = spheres[index].pos - s0->pos; 
        s0->radius = v1.magnitude()-spheres[index].radius; 
    }
    if(num_bubbles==5)
    {
      
    }
    if(index==-1)return bubbles[0];
    if(bubble)
    {
        s0->neighboors.push_back(num_spheres+index);
        return bubbles[index];
    } 
    else 
    {
        s0->neighboors.push_back(index);
        return spheres[index];
    }
}

sphere stage2(bubble * s0, sphere * s1, int num_spheres, int num_bubbles)
{
    double min_radius=10;
    int index = -1;
    bool bubble = false;
    vec3 t = (s0->pos - s1->pos).normalize();
    vec3 c1 = s1->pos+scalar_product(t,s1->radius);
    vec3 v1c = c1  - s1->pos; 

    for(int i=0; i <=num_spheres+num_bubbles; i++)
    {

        bool skip=false;
        for(int j: s0->neighboors){
            if(j==i)
            {
                skip=true;
                break;
            }
        }
        if(skip) continue;

        sphere sc;  //candidate for s2
        //!!!!!!!!!
        if(i>=num_spheres) sc = bubbles[i-num_spheres];
        else sc = spheres[i];
        vec3 vs1 = s1->pos - sc.pos;
        double r0c = ( pow(s1->radius, 2) - pow(sc.radius, 2) 
            + pow(vs1.magnitude(), 2) + 2*dot_product(v1c, vs1)) / 
            (2*(sc.radius - s1->radius - dot_product(v1c, vs1)/s1->radius));
        if(r0c < min_radius && r0c>s0->radius)

        {
	    
            min_radius = r0c;
            if(i>=num_spheres)
            {
                bubble = true;
                index = i-num_spheres;
            }
            else index =i;
        }
    }

    if(index==-1)return bubbles[0];
    vec3 min_vec = (vec3){-s0->pos.x, 0, 0};
    min_vec = s0->pos.y > min_vec.magnitude() ? (vec3){0, -s0->pos.y, 0}: min_vec; 
    min_vec = s0->pos.z > min_vec.magnitude() ? (vec3){0, 0, -s0->pos.z}: min_vec; 
    min_vec = x_max - s0->pos.x > min_vec.magnitude() ? (vec3){x_max - s0->pos.x, 0, 0}: min_vec ; 
    min_vec = y_max - s0->pos.y > min_vec.magnitude() ? (vec3){0, y_max - s0->pos.y, 0}: min_vec ; 
    min_vec = z_max - s0->pos.z > min_vec.magnitude() ? (vec3){0, 0, z_max - s0->pos.z}: min_vec ; 
    if(min_vec.magnitude()<min_radius)
    {
        min_radius=min_vec.magnitude();
        s0->radius = min_radius;
        s0->pos.equals(s1->pos + v1c.scalar_multiply((min_radius+s1->radius)/s1->radius));
        p_spheres.push_back(sphere{0, s0->pos-min_vec});
        return p_spheres.back();
    }

    sphere sc;
    if(bubble)sc = bubbles[index];
    else sc=spheres[index];
    s0->radius = min_radius;
    s0->pos.equals(s1->pos + v1c.scalar_multiply((min_radius+s1->radius)/s1->radius));
    if(bubble)
    {
        s0->neighboors.push_back(num_spheres+index);
        return bubbles[index];
    } 
    else 

    {
        s0->neighboors.push_back(index);
        return spheres[index];
    }
}


sphere stage3(bubble * s0, sphere * s1, sphere * s2, int num_spheres, int num_bubbles)
{
    double min_radius=10;
    double min_dist=10;
    int index = -1;
    bool bubble = false;

    vec3 vc1, b1, b2;
    vec3 beta;
    double r1 = s1->radius;
    double r2 = s2->radius;
//    double * rs = s0->radius;

    vec3 x2 = s2->pos;
    vec3 x1 = s1->pos;   
    
    vc1 = ((s0->pos - s1->pos).normalize()).scalar_multiply(s1->radius);
    
    b1 = s0->pos    - x1; 
    b2 = x2         - x1;

    vec3 vc1xb2 = cross_product(vc1, b2);
    Eigen::Matrix4f A;
    Eigen::Vector4f b;
    Eigen::Vector4f ps1;
    ps1<< s1->pos.x, s1->pos.y, s1->pos.z , 0;
    for(int i=0; i <num_spheres+num_bubbles; i++)
    {
        bool skip=false;
        for(int j: s0->neighboors){
            if(j==i)
            {
                skip=true;
                break;
            }
        }
        if(skip) continue;
        sphere sc;  //candidate for s2
        if(i>num_spheres) sc = bubbles[i-num_spheres];
        else sc = spheres[i];
        double rs = sc.radius;
        vec3 bs = sc.pos - x1;
       
        DEBUG(2,"FOR BUBBLE "<<i <<"\n");
        DEBUG(3, "sc is at " << sc.pos.x << " "<< sc.pos.y << " "<< sc.pos.z << std::endl);
        DEBUG(3, "s1 is at " << x1.x << " "<< x1.y << " "<< x1.z << std::endl);
        DEBUG(3, "s2 is at " << x2.x << " "<< x2.y << " "<< x2.z << std::endl);
        DEBUG(3, "vc1 is  " << vc1.x << " "<< vc1.y << " "<< vc1.z << std::endl);
        A   <<   (bs-b2).x     , (bs-b2).y    ,(bs-b2).z     , (rs-r2)
                , b2.x         , b2.y         , b2.z         , (r2-r1)
                , bs.x         , bs.y         , bs.z         , (rs-r1)
                , vc1xb2.x     , vc1xb2.y     , vc1xb2.z     , 0      ;
        DEBUG(3,  "A: \n" << A <<std::endl);
        b   << (pow(r2, 2) - pow(rs,2) -pow(b2.magnitude(), 2) +pow(bs.magnitude(), 2))/2.0
            , (pow(b2.magnitude(), 2) + pow(r1, 2) - pow(r2,2))/2.0
            , (pow(bs.magnitude(), 2) + pow(r1, 2) - pow(rs,2))/2.0
            , 0 ;
        DEBUG(3, "b: \n " <<b <<std::endl);
        Eigen::Vector4f soln = A.colPivHouseholderQr().solve(b);
        
        DEBUG(3,"soln before radius correction \n" << soln <<std::endl << "corresponding to " << ps1+soln << std::endl);
	if(soln[3]==0)
        {
            soln[3] 
                = pow(pow(bs.x-soln[0],2)+pow(bs.y-soln[1],2)+pow(bs.z-soln[2],2),0.5)-sc.radius;
        }
        DEBUG(3,"soln \n" << soln <<std::endl << "corresponding to " << ps1+soln << std::endl);
        if(soln[3] < min_radius && soln[3]>s0->radius)
        //if(soln.segment(0,2).norm()<min_dist && soln[3]>s0->radius)
        {
            min_dist =soln.segment(0,2).norm(); 
            min_radius = soln[3];
	    beta={soln[0], soln[1], soln[2]};
            if(i>num_spheres)
            {
                bubble = true;
                index = i-num_spheres;
            }
            else index =i;
        }
    }
    if(index==-1)return bubbles[0];
    sphere sc;
    if(bubble)sc = bubbles[index];
    else sc=spheres[index];
    s0->radius = min_radius;
    s0->pos.equals(s1->pos + beta);
    if(bubble)
    {
        s0->neighboors.push_back(num_spheres+index);
        return bubbles[index];
    } 
    else 
    {
        s0->neighboors.push_back(index);
        return spheres[index];
    }
}

sphere stage4(bubble * s0, sphere * s1, sphere * s2, sphere * s3, int num_spheres, int num_bubbles)
{
    double min_radius;
    int index = -1;
    bool bubble = false;

    vec3 vc1, b1, b2, b3;
    vec3 beta;
    double r1 = s1->radius;
    double r2 = s2->radius;
    double r3 = s3->radius;
//    double * rs = s0->radius;

    vec3 x1 = s1->pos;   
    vec3 x2 = s2->pos;
    vec3 x3 = s3->pos;   
    
    vec3 c1 = s1->pos+((s0->pos - s1->pos).normalize()).scalar_multiply(s1->radius);
    vc1 = c1  - s1->pos; 
    
    b1 = s0->pos    - x1; 
    b2 = x2         - x1;
    b3 = x3         - x1;

    vec3 vc1xb2 = cross_product(vc1, b2);
    Eigen::Matrix4f A;
    Eigen::Vector4f b;
    for(int i=0; i <num_spheres+num_bubbles; i++)
    {
        bool skip=false;
        for(int j: s0->neighboors){
            if(j==i)
            {
                skip=true;
                break;
            }
        }
        if(skip) continue;
        sphere sc;  //candidate for s2
        if(i>num_spheres) sc = bubbles[i-num_spheres];
        else sc = spheres[i];
        double rs = sc.radius;
        vec3 bs = sc.pos - x1;
        A   << (b2-b3).x   , (b2-b3).y    , (b2-b3).z    , (r2-r3)
            , b2.x         , b2.y         , b2.z         , (r2-r1)
            , b3.x         , b3.y         , b3.z         , (r3-r1)
            , bs.x         , bs.y         , bs.z         , (rs-r1);

        b   << (pow(b2.magnitude(), 2) - pow(b3.magnitude(), 2) 
                    + pow(rs, 2) - pow(r2,2))/2.0
            , (pow(b2.magnitude(), 2) + pow(r1, 2) - pow(r2,2))/2.0
            , (pow(b3.magnitude(), 2) + pow(r1, 2) - pow(rs,2))/2.0
            , (pow(bs.magnitude(), 2) + pow(r1, 2) - pow(r3,2))/2.0;
        Eigen::Vector4f soln = A.colPivHouseholderQr().solve(b);
        if(soln[3] < min_radius && soln[3]>s0->radius)
        {
            min_radius = soln[3];
	    beta={soln[0], soln[1], soln[2]};
            if(i>num_spheres)
            {
                bubble = true;
                index = i-num_spheres;
            }
            else index =i;
        }
    }
    if(index==-1)return bubbles[0];
    sphere sc;
    if(bubble)sc = bubbles[index];
    else sc=spheres[index];
    s0->radius = min_radius;
    s0->pos = s1->pos + beta;
    if(bubble)
    {
        s0->neighboors.push_back(num_spheres+index);
        return bubbles[index];
    } 
    else 
    {
        s0->neighboors.push_back(index);
        return spheres[index];
    }
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  read_sphere_file
 *  Description:  Reads sphere coordinates from file
 * =====================================================================================
 * 
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


int check_errors(bubble * b, int num_spheres, int num_bubbles)
{
    
    int err= check_intersect(b->radius, b->pos, num_spheres, num_bubbles, true, true, true);
    if(err!=-1)
    {
        if(err>num_spheres)
        {
            err-=num_spheres;
            std::cout << "Error: This bubble intersected with bubble " << err << std::endl;
        }
        else
        {
            std::cout << "Error: This bubble intersected with sphere " << err << std::endl;
        }
        return err;
    }
    return -1;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  cage
 *  Description:  Generate a sphere cage around the container
 * =====================================================================================
 */
int sphere_cage_gen(int num_spheres, double r)
{
    for(int i =0; i<=x_max/(2*r) + 1; i ++)
    {
    for(int j =0; j<=y_max/(2*r) + 1; j ++)
    {
        num_spheres++;
        spheres[num_spheres]=(sphere){r, (vec3){i*2*r, j*2*r, -r}};
        num_spheres++;
        spheres[num_spheres]=(sphere){r, (vec3){i*2*r, i*2*r, z_max+r}};
        
        num_spheres++;
        spheres[num_spheres]=(sphere){r, (vec3){i*2*r, -r, j*2*r}};
        
        num_spheres++;
        spheres[num_spheres]=(sphere){r, (vec3){i*2*r, y_max+r, j*2*r}};
        
        num_spheres++;
        spheres[num_spheres]=(sphere){r, (vec3){-r, i*2*r, j*2*r}};
        
        num_spheres++;
        spheres[num_spheres]=(sphere){r, (vec3){x_max+r, i*2*r, j*2*r}};
    }
    }
    return num_spheres;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  main loop
 * =====================================================================================
 */
int main(int argc, char * argv[])
{
    bool debug = true;
    std::cout << std::fixed;
    if(argc!=4)
    {
        //std::cout << 
        //    "Usage: ./bubble-fill [input-sphere-file] [output-bubble-file]" <<
        //    std::endl;
        argv[1]="spheres";
	argv[2]="bubbles";
	argv[3]="0";
    }
    std::string sphere_file_path = argv[1];
    std::string bubble_file_path = argv[2];
    debug_level = atoi(argv[3]);
    int num_spheres = read_sphere_coords(sphere_file_path);
    //num_spheres =sphere_cage_gen(num_spheres, 0.1);

    
    //Initialize the null sphere
    for(int i=1; i<20; i++)
    { 
        bubble b;
        double r, x, y, z;
        while(true)
        {
            r = 0.00;
            x = rand_range(0.001,9.999);
            y = rand_range(0.001,9.999);
            z = rand_range(0.001,9.999);
            if(check_intersect(r, (vec3){x, y, z}, num_spheres, i, true, true, true)
                    ==-1) break;
        }
        b.id=i;
        b.radius = r;
        b.pos.x = x;
        b.pos.y = y;
        b.pos.z = z;
         //if(debug) std::cout<< "Pre : "<< b.radius << " " << 
         //         b.pos.x << " " << b.pos.y << " " << b.pos.z << std::endl;
        sphere c1 = stage1(&b, num_spheres, i);
        DEBUG(1,"S1: " << b.radius << " " << 
                  b.pos.x << " " << b.pos.y << " " << b.pos.z << std::endl);
        sphere c2 = stage2(&b, &c1, num_spheres, i);
          if(b.radius<0 || b.pos.x<0 || b.pos.y<0 || b.pos.z<0 || b.pos.x>10 & b.pos.y>10 || b.pos.z>10){i--; continue;}
        DEBUG(1,"S2: " << b.radius << " " << 
                  b.pos.x << " " << b.pos.y << " " << b.pos.z << std::endl);
        sphere c3 = stage3(&b, &c1, &c2, num_spheres, i);

        DEBUG(1,"S3: " << b.radius << " " << 
                  b.pos.x << " " << b.pos.y << " " << b.pos.z << std::endl);
        sphere c4 = stage4(&b, &c1, &c2, &c3, num_spheres, i);

        DEBUG(1,"S4: " << b.radius << " " << 
                  b.pos.x << " " << b.pos.y << " " << b.pos.z << std::endl);
        /*if(check_intersect(b.radius, b.pos, num_spheres, i, true, true, true)!=-1)
        {
            i--;
            continue;
        }*/
          if(b.radius>0 && b.pos.x>0 && b.pos.y>0 && b.pos.z>0 && b.pos.x<10 & b.pos.y<10 && b.pos.z<10)
          {
            std::cout << b.radius << " " <<  b.pos.x << " " << b.pos.y << " " << b.pos.z << std::endl;
            for(int j: b.neighboors){
                if(j-num_spheres>0)
                {
                    std::cout<< j-num_spheres << " " ;
                }
            }
            std::cout<< std::endl ;
            bubbles[i]=b;
          }
          else {
         DEBUG(1, "ERROR::OOB: " << b.radius << " " << 
                  b.pos.x << " " << b.pos.y << " " << b.pos.z << std::endl);
            i--;
          }
//        std::cout<< std::endl ;
    }
}
