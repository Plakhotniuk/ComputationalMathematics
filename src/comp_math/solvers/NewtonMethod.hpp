//
// Created by Арсений Плахотнюк on 04.12.2022.
//

#ifndef COMPUTATIONALMATHEMATICS_NEWTONMETHOD_HPP
#define COMPUTATIONALMATHEMATICS_NEWTONMETHOD_HPP
#include "comp_math/Exceptions/SlaeBaseException.hpp"
#include <vector>

namespace Slae::Solvers{
    double solverNewtonMethod(const double init_approach, const double step, std::function<double(double)>& f,
                              const double tolerance = 1.e-10){
        double x_n = init_approach;
        std::function<double(double)> f_der = [f, step](double x) { return (f(x + step) - f(x))/step; };
        while(std::abs(1 / f_der(x_n) * f(x_n)) > tolerance){
            x_n = x_n - f(x_n)/f_der(x_n);
        }
        return x_n;
    }
}


#endif //COMPUTATIONALMATHEMATICS_NEWTONMETHOD_HPP
