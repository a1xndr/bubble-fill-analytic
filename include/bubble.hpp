/*
 * =====================================================================================
 *
 *       Filename:  sphere.hpp
 *
 *    Description:  Object definition for sphere
 *
 *        Version:  1.0
 *        Created:  06/02/2015 10:42:48 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef BUBBLE_HPP
#define BUBBLE_HPP

#include "sphere.hpp"
#include <vector>
struct bubble: sphere {
    int id=-1;
    std::vector<int> neighboors;
    
};


#endif

