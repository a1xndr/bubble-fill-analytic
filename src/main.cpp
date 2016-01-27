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



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  stage1(bubble *s, int i)
 *  Description:  Grows a sphere not in contact with any other bubble, until
 *  a collision
 * =====================================================================================
 */
sphere stage1(bubble *s, int i)
{
    
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  read_sphere_file
 *  Description:  Reads sphere coordinates from file
 * =====================================================================================
 */
int read_sphere_coords(){
    std::ifstream file("spheres.list");
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
