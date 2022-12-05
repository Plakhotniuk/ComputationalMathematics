//
// Created by Арсений Плахотнюк on 04.12.2022.
//

#ifndef COMPUTATIONALMATHEMATICS_RUNGEKUTTA4_HPP
#define COMPUTATIONALMATHEMATICS_RUNGEKUTTA4_HPP
#include "comp_math/Exceptions/SlaeBaseException.hpp"
#include "comp_math//types/BasicTypes.hpp"
#include <vector>

std::vector<double> RungeKutta4(std::function<double(double, double)>& f, const double x_0, const double x_f,
                                const double y_0, const double num_of_steps){
    double k1, k2, k3, k4;
    std::vector<double> y = {y_0};
    const double h = (x_f - x_0)/num_of_steps;
    double x_n = x_0;
    for(int i = 0; i < num_of_steps; ++i){
        x_n += h;
        k1 = f(x_n, y[i]);
        k2 = f(x_n + h/2, y[i] + k1 * h / 2);
        k3 = f(x_n + h/2, y[i] + k2 * h / 2);
        k4 = f(x_n + h, y[i] + k3 * h);
        y.push_back(y[i] + h*(k1/6 + k2/3 + k3/3 + k4/6));
    }
    return y;
}

std::vector<VectorXd> RungeKutta4Vec(std::function<VectorXd(double, const VectorXd&)>& f, const double x_0,
                                     const double x_f, const VectorXd& y_0, const int num_of_steps){
    VectorXd k1, k2, k3, k4;
    const double h = (x_f - x_0)/num_of_steps;
    std::vector<VectorXd> y = {y_0};
    double x_n = x_0;
    for(int i = 0; i < num_of_steps; ++i){
        x_n += h;
        k1 = f(x_n, y[i]);
        k2 = f(x_n + h/2, y[i] + k1 * h / 2);
        k3 = f(x_n + h/2, y[i] + k2 * h / 2);
        k4 = f(x_n + h, y[i] + k3 * h);
        y.emplace_back(y[i] + h*(k1/6 + k2/3 + k3/3 + k4/6));
    }
    return y;
}

#endif //COMPUTATIONALMATHEMATICS_RUNGEKUTTA4_HPP
