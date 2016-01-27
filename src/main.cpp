/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  
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
