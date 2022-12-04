//
// Created by Арсений Плахотнюк on 04.12.2022.
//

#ifndef HOMETASK2_LINEAR_HPP
#define HOMETASK2_LINEAR_HPP
#include <vector>
#include "comp_math/matrix/ThreeDiagonalMatrix.hpp"
#include "comp_math/Exceptions/SlaeBaseException.hpp"
#include <sstream>

double linear_interpolation(std::vector<double> x, std::vector<double> y, double point){
    if(x[0] <= point && point < x[x.size() - 1]){
        for(int i = 0; i < x.size() - 1; ++i){
            if(x[i] <= point && point < x[i+1])
                return y[i] + (y[i+1] - y[i])/(x[i+1] - x[i])*(point - x[i]);
        }
    }
    else if (x[x.size() - 1] == point){
        return y[y.size() - 1];
    }
    else{
        std::stringstream buff;
        buff << "Value is out of interpolation range!!! " << ". File: " << __FILE__ << ". Line: " << __LINE__;
        throw Slae::SlaeBaseExceptionCpp(buff.str());
    }
}

#endif //HOMETASK2_LINEAR_HPP
