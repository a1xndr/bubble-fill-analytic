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
    if(index==-1){
        DEBUG(2,"No contact"<<std::endl);
        return bubbles[0];
    }
    if(bubble)
    {
        DEBUG(2,"Contact with bubble: " << index <<std::endl);
        s0->neighboors.push_back(num_spheres+index);
        return bubbles[index];
    } 
    else 
    {
        DEBUG(2,"Contact with sphere: " << index <<std::endl);
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

    if(index==-1){
        DEBUG(2,"No contact" <<std::endl);
        return bubbles[0];
    }
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
        DEBUG(2,"Contact with bubble: " << index <<std::endl);
        s0->neighboors.push_back(num_spheres+index);
        return bubbles[index];
    } 
    else 
    {
        DEBUG(2,"Contact with sphere: " << index <<std::endl);
        s0->neighboors.push_back(index);
        return spheres[index];
    }
}


sphere stage3(bubble * s0, sphere * s1, sphere * s2, int num_spheres, int num_bubbles)
{
    int index = -1;
    
    double radius_min=x_max;
    double r1 = s1->radius;
    double r2 = s2->radius;
    
    bool bubble = false;
    
    vec3 vc1, b1, b2;

    Eigen::Vector3f beta_min, beta;
    Eigen::Vector3f b0e, b2e;
    
    Eigen::Matrix3f B;
    Eigen::Vector3f r;
    Eigen::Vector3f s;

    vec3 x2 = s2->pos;
    vec3 x1 = s1->pos;   
    
    vc1 = ((s0->pos - s1->pos).normalize()).scalar_multiply(s1->radius);
    
    b1 = s0->pos    - x1; 
    b2 = x2         - x1;

    b0e <<   b1.x      , b1.y     , b1.z; 
    b2e <<   b2.x      , b2.y     , b2.z;

    
    DEBUG(3, "b0e " << b0e << std::endl);
    DEBUG(3, "b2e is " << b2e << std::endl);
    Eigen::Vector3f b0_proj = b2e * b0e.dot(b2e) / (pow(b2e.norm(),2));
    DEBUG(3, "b0_proj is " << b0_proj << std::endl);
    
    vec3 vc1xb2 = cross_product(vc1, b2);
    for(int i=0; i <num_spheres+num_bubbles; i++)
    {
        bool skip=false;
        bool check_both=true;
        for(int j: s0->neighboors){
            if(j==i)
            {
                skip=true;
                break;
            }
        }
        if(skip) continue;
        sphere sc;  //candidate for s2
        if(i>num_spheres){
            DEBUG(2,"FOR BUBBLE "<<i-num_spheres <<"\n");
            sc = bubbles[i-num_spheres];
        }
        else{
            DEBUG(2,"FOR SPHERE "<<i <<"\n");
            sc = spheres[i];
        }
        double rs = sc.radius;
        vec3 bs = sc.pos - x1;
       
        DEBUG(3, "sc is at " << sc.pos.x << " "<< sc.pos.y << " "<< sc.pos.z << std::endl);
        DEBUG(3, "s1 is at " << x1.x << " "<< x1.y << " "<< x1.z << std::endl);
        DEBUG(3, "s2 is at " << x2.x << " "<< x2.y << " "<< x2.z << std::endl);
        DEBUG(3, "vc1 is  " << vc1.x << " "<< vc1.y << " "<< vc1.z << std::endl);
     
        /*-----------------------------------------------------------------------------
         *  Here we have a mix of quadratic and linear equations
         *  Our linear condition is:
         *
         *  mat3(B) * vec3(beta) = vec3(s) - vec3(r)*rf 
         *  
         *  The quadratic condition is 
         *  
         *  |vec3(beta)|**2 = (rf+r1)^2
         *
         *  We are solving for rf and then calculating x, y, z components of beta
         *
         *  The quadratic equation is of form ar^2 + br + c = 0 where
         *
         *  a = |inv(B)vec3(r)|**2 - 1
         *  b = 2*inv(B)vec3(r).inv(B)vec3(s) - 2*r1
         *  c = |inv(B)vec3(s)|**2 - r1**2
         *
         *  so
         *            -b +/- sqrt(b^2-4ac)
         *      rf = -----------------------
         *                  2a
         *-----------------------------------------------------------------------------*/
        
        B   <<    b2.x         , b2.y         , b2.z         
                , bs.x         , bs.y         , bs.z         
                , vc1xb2.x     , vc1xb2.y     , vc1xb2.z     ;

        r   <<   r2-r1      , rs-r1     , 0;
        s   <<   (pow(b2.magnitude(), 2) + pow(r1, 2) - pow(r2,2))/2.0
            ,    (pow(bs.magnitude(), 2) + pow(r1, 2) - pow(rs,2))/2.0
            ,    0 ;

        DEBUG(4,"B is: " << B <<std::endl);
        DEBUG(4,"r is: " << r <<std::endl);
        DEBUG(4,"s is: " << s <<std::endl);
        double a, b, c;     //Coefficients of the quadratic
        a = (B.inverse()*r).squaredNorm() -1;
        b = -2*(B.inverse()*r).dot(B.inverse()*s)-2*r1;
        c = (B.inverse()*s).squaredNorm() - r1*r1;

        double sol1, sol2, sol;  //Solutions of the quadratic
        sol1 = (-b + sqrt(b*b-4*a*c))/(2*a);
        sol2 = (-b - sqrt(b*b-4*a*c))/(2*a);

        DEBUG(3,"Radii solutions: " << sol1 << " " << sol2 <<std::endl);
        
        if(sol1<0)
        {
            sol=sol2;
        }
        else if(sol2<0)
        {
            sol=sol1;
        }
        else if(sol1<0 && sol2<0)
        {
            continue;
        }
        else {
            DEBUG(3,"Both radii are positive " <<std::endl);
            sol=std::max(sol1, sol2);
            if(sol <= radius_min && sol > s0->radius)
            {
                /*-----------------------------------------------------------------------------
                *  Now that we have the radius, we find the beta x, y and z  coordinates
                *  by plugging r into the linear system B*beta = vec3(s) - vec3(r)*sol
                *-----------------------------------------------------------------------------*/
                beta = B.colPivHouseholderQr().solve(s-r*sol);
                DEBUG(3,"Beta solution: " << beta <<std::endl);
                if((beta-b0_proj).dot(b0e-b0_proj)<=SMIDGE){
                    DEBUG(3,"Rejected because of negative dot product: " << (beta-b0_proj).dot(b0e-b0_proj) <<std::endl);
                }
                else {
                beta_min = beta;
                radius_min = sol;


                if(i>num_spheres)
                {
                    bubble = true;
                    index = i-num_spheres;
                }
                else index =i;
                }
            }
            sol=std::min(sol1, sol2);
        }

        //Check if  radii are possible solutions
        if(sol <= radius_min && sol > s0->radius)
        {
            /*-----------------------------------------------------------------------------
            *  Now that we have the radius, we find the beta x, y and z  coordinates
            *  by plugging r into the linear system B*beta = vec3(s) - vec3(r)*sol
            *-----------------------------------------------------------------------------*/
            beta = B.colPivHouseholderQr().solve(s-r*sol);
            DEBUG(3,"Beta solution: " << beta <<std::endl);
            if((beta-b0_proj).dot(b0e-b0_proj)<=SMIDGE){
                
                DEBUG(3,"Rejected because of negative dot product: " << (beta-b0_proj).dot(b0e-b0_proj) <<std::endl);
                continue;
            }
            beta_min = beta;
            radius_min = sol;


            if(i>num_spheres)
            {
                bubble = true;
                index = i-num_spheres;
            }
            else index =i;
        }
    }
    if(index==-1)
    {
        DEBUG(2,"No contact: " <<std::endl);
        return bubbles[0];
    }
    sphere sc;
    if(bubble)sc = bubbles[index];
    else sc=spheres[index];
    s0->radius = radius_min;
    vec3 bm;
    bm.x = beta_min(0);
    bm.y = beta_min(1);
    bm.z = beta_min(2);
    DEBUG(2,"bm: " << bm.x << " " << bm.y << " " << bm.z <<std::endl);
    DEBUG(2,"s0.pos: " << s0->pos.x << " " << s0->pos.y << " " << s0->pos.z <<std::endl);
    s0->pos.x = s1->pos.x + bm.x;
    s0->pos.y = s1->pos.y + bm.y;
    s0->pos.z = s1->pos.z + bm.z;
    DEBUG(2,"s0.pos: " << s0->pos.x << " " << s0->pos.y << " " << s0->pos.z <<std::endl);
    if(bubble)
    {
        DEBUG(2,"Contact with bubble: " << index <<std::endl);
        s0->neighboors.push_back(num_spheres+index);
        return bubbles[index];
    } 
    else 
    {
        DEBUG(2,"Contact with sphere: " << index <<std::endl);
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

int check_oob(bubble *b)
{
    if(b->radius>0 && b->pos.x>0 && b->pos.y>0 && b->pos.z>0 && b->pos.x<10 & b->pos.y<10 && b->pos.z<10) return false;
    return true;
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
//        std::cout << spheres[num_spheres].radius << " " <<  spheres[num_spheres].pos.x << " " << spheres[num_spheres].pos.y << " " << spheres[num_spheres].pos.z << std::endl;
        num_spheres++;
        spheres[num_spheres]=(sphere){r, (vec3){i*2*r, j*2*r, z_max+r}};
//        std::cout << spheres[num_spheres].radius << " " <<  spheres[num_spheres].pos.x << " " << spheres[num_spheres].pos.y << " " << spheres[num_spheres].pos.z << std::endl;
        
        num_spheres++;
        spheres[num_spheres]=(sphere){r, (vec3){i*2*r, -r, j*2*r}};
//        std::cout << spheres[num_spheres].radius << " " <<  spheres[num_spheres].pos.x << " " << spheres[num_spheres].pos.y << " " << spheres[num_spheres].pos.z << std::endl;
        
        num_spheres++;
        spheres[num_spheres]=(sphere){r, (vec3){i*2*r, y_max+r, j*2*r}};
//        std::cout << spheres[num_spheres].radius << " " <<  spheres[num_spheres].pos.x << " " << spheres[num_spheres].pos.y << " " << spheres[num_spheres].pos.z << std::endl;
        
        num_spheres++;
        spheres[num_spheres]=(sphere){r, (vec3){-r, i*2*r, j*2*r}};
//        std::cout << spheres[num_spheres].radius << " " <<  spheres[num_spheres].pos.x << " " << spheres[num_spheres].pos.y << " " << spheres[num_spheres].pos.z << std::endl;
        
        num_spheres++;
        spheres[num_spheres]=(sphere){r, (vec3){x_max+r, i*2*r, j*2*r}};
//        std::cout << spheres[num_spheres].radius << " " <<  spheres[num_spheres].pos.x << " " << spheres[num_spheres].pos.y << " " << spheres[num_spheres].pos.z << std::endl;
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
    if(argc<3)
    {
        argv[1]="spheres";
	argv[2]="bubbles";
	argv[3]="0";
    }
    else if(argc<4)
    {
	argv[3]="0";
    }
    std::string sphere_file_path = argv[1];
    std::string bubble_file_path = argv[2];
    debug_level = atoi(argv[3]);
    int num_spheres = read_sphere_coords(sphere_file_path);
    //num_spheres =sphere_cage_gen(num_spheres, 0.25);
   srand (123); 
    bubbles[0].radius=-20;
    //Initialize the null sphere
    for(int i=1; i<3; i++)
    { 
        bool skip = false;
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
        if(c1.radius==-20 || check_oob(&b)) skip=true;
        DEBUG(1,"S1: " << b.radius << " " << 
                  b.pos.x << " " << b.pos.y << " " << b.pos.z << std::endl);
        if(!skip){
            sphere c2 = stage2(&b, &c1, num_spheres, i);
            if(c2.radius==-20 || check_oob(&b)) skip=true;
            DEBUG(1,"S2: " << b.radius << " " << 
                    b.pos.x << " " << b.pos.y << " " << b.pos.z << std::endl);
            if(!skip)
            {
                sphere c3 = stage3(&b, &c1, &c2, num_spheres, i);
                if(c3.radius==-20 || check_oob(&b)) skip=true;
                DEBUG(1,"S3: " << b.radius << " " << 
                        b.pos.x << " " << b.pos.y << " " << b.pos.z << std::endl);
                if(!skip) 
                {
                    sphere c4 = stage4(&b, &c1, &c2, &c3, num_spheres, i);
                    DEBUG(1,"S4: " << b.radius << " " << 
                            b.pos.x << " " << b.pos.y << " " << b.pos.z << std::endl);
                }
            }
        }

        //sphere c4 = stage4(&b, &c1, &c2, &c3, num_spheres, i);

        //DEBUG(1,"S4: " << b.radius << " " << 
        //          b.pos.x << " " << b.pos.y << " " << b.pos.z << std::endl);
          if(b.radius>0 && b.pos.x>0 && b.pos.y>0 && b.pos.z>0 && b.pos.x<10 & b.pos.y<10 && b.pos.z<10)
          {
            if(check_intersect(b.radius, b.pos, num_spheres, i, true, true, true)!=-1)
            {
            DEBUG(1, "ERROR::INTERSECT: " << check_intersect(b.radius, b.pos, num_spheres, i, true, true, true) <<std::endl;);
                break;
            i--;
            continue;
            }
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
