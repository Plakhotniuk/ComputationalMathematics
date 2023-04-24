//
// Created by Арсений Плахотнюк on 04.12.2022.
//

#ifndef COMPUTATIONALMATHEMATICS_SHOOTINGMETHOD_HPP
#define COMPUTATIONALMATHEMATICS_SHOOTINGMETHOD_HPP
#include "compMath/Exceptions/SlaeBaseException.hpp"
#include "compMath/solvers/NewtonMethod.hpp"
#include "compMath/integrators/RungeKutta4.hpp"
#include <vector>

namespace Slae::Solvers{
    std::vector<VectorXd> solveShootingMethod(std::function<double(double, const VectorXd&)>& f, const int num_of_steps,
                                            const double x_0, const double x_f, const double y_0, const double y_f,
                                            const double tolerance = 1.e-10){
        // f - скалярная функция f(x, y, y')
        // y_ = [y, y']
        double alpha =  (y_f - y_0) / (x_f - x_0); // начальный пристрелочный параметр
        double p_0 = alpha;
        VectorXd y_vec_0(2); // вектор начальных условий задачи Коши для системы уравнений
        y_vec_0 << y_0, p_0;

        std::function<VectorXd(double, const VectorXd&)> f_vec = [f](double x, const VectorXd& y_){
            VectorXd res(2);
            res << y_[1], f(x, y_);
            return res;
        };
        double h = (x_f - x_0) / num_of_steps;

        double F_alpha_n = RungeKutta4Vec(f_vec, x_0, x_f, y_vec_0, num_of_steps)[num_of_steps](0) - y_f;
        p_0 = alpha + h;
        y_vec_0(0) = y_0;
        y_vec_0(1) = p_0;

        double F_alpha_n_h = RungeKutta4Vec(f_vec, x_0, x_f, y_vec_0, num_of_steps)[num_of_steps](0) - y_f;

        while(h * F_alpha_n / (F_alpha_n_h - F_alpha_n) > tolerance){
            alpha = alpha - h * F_alpha_n / (F_alpha_n_h - F_alpha_n);

            F_alpha_n = RungeKutta4Vec(f_vec, x_0, x_f, y_vec_0, num_of_steps)[num_of_steps](0) - y_f;
            p_0 = alpha + h;
            y_vec_0(0) = y_0;
            y_vec_0(1) = p_0;
            F_alpha_n_h = RungeKutta4Vec(f_vec, x_0, x_f, y_vec_0, num_of_steps)[num_of_steps](0) - y_f;
        }

        y_vec_0(0) = y_0;
        y_vec_0(1) = alpha;
        std::vector<VectorXd> answer = RungeKutta4Vec(f_vec, x_0, x_f, y_vec_0, num_of_steps);
        return answer;
    };
}

#endif //COMPUTATIONALMATHEMATICS_SHOOTINGMETHOD_HPP
