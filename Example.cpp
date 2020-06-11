/*
 * Example.cpp
 *
 *  Created on: Jun 11, 2020
 *      Author: mirsandiharyo
 */


#include <iostream>
#include "CustomArray.h"

int main() {
    std::cout << "testing custom array" << std::endl;

    csar::variables<double, 32> data;
    data.allocate(2, 2, 2);

    for(uint32_t i = 0; i < data.NI_init(); ++i)
        for(uint32_t j = 0; j < data.NJ_init(); ++j)
            for(uint32_t k = 0; k < data.NK_init(); ++k)
                data(i, j, k) = 0.5 * (i + j + k);

    for(uint32_t n = 0; n < data.size(); ++n)
        std::cout << data(n) << " ";
    std::cout << std::endl;

    csar::variables<double, 32> data2;
    data2.allocate(2, 2, 2);
    data2 = data;
    for(uint32_t n = 0; n < data2.size(); ++n)
        std::cout << data2(n) << " ";
    std::cout << std::endl;

    csar::variables<int> *vec;
    vec = new csar::variables<int>;
    vec->allocate(10);
    for(uint i = 0; i < vec->NI_init(); ++i)
    	vec->at(i) = i;
    for(uint n = 0; n < vec->size(); ++n)
        std::cout << vec->at(n) << " ";
    std::cout << std::endl;
    delete vec;

    std::cout << "finished" << std::endl;
    return 0;
}

