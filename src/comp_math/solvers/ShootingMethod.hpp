//
// Created by Арсений Плахотнюк on 04.12.2022.
//

#ifndef COMPUTATIONALMATHEMATICS_SHOOTINGMETHOD_HPP
#define COMPUTATIONALMATHEMATICS_SHOOTINGMETHOD_HPP
#include "comp_math/Exceptions/SlaeBaseException.hpp"
#include "comp_math/solvers/NewtonMethod.hpp"
#include <vector>

namespace Slae::Solvers{
    std::vector<double> solveShootingMethod(std::function<double(double, double, double)>& f, const double h,
                                            const double x_0, const double x_f, const double y_0, const double y_f,
                                            const double tolerance = 1.e-10){
        double alpha =  (y_f - y_0) / (x_f - x_0);

        alpha = solverNewtonMethod(alpha, h, );
    };
}

#endif //COMPUTATIONALMATHEMATICS_SHOOTINGMETHOD_HPP
