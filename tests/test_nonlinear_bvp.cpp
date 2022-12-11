//
// Created by Арсений Плахотнюк on 11.12.2022.
//
#include "gtest/gtest.h"
#include "comp_math/types/BasicTypes.hpp"
#include "comp_math/solvers/ShootingMethod.hpp"
#include "comp_math/matrix/NonLinearBoundaryValueProblemMatrix3.hpp"
#include <iostream>
#include <fstream>


TEST(SHOOTING4, METHOD4){
    //Задача 4

    // скалярная функция правой части уравнения
    std::function<double(double, const VectorXd&)> f = [](double x, const VectorXd& y_){
        return sqrt(-y_(0)*exp(y_(1)) + exp(1)/log(x)*y_(0)*y_(0) + 1/(x*x));
    };
    int num_of_steps = 10000;

    double x_0 = exp(1); // граничные условия
    double y_0 = exp(1);
    double x_f = exp(1)*exp(1);
    double y_f = 2.* exp(1)*exp(1);
    double tolerance = 1.e-10;
    double h = (x_f - x_0) / num_of_steps;
    std::vector<VectorXd> answer = Slae::Solvers::solveShootingMethod(f, num_of_steps, x_0, x_f, y_0, y_f, tolerance);
    std::fstream file;
    file.open("test_nonlin_shooting_method.txt", std::fstream::out);

    for(int i = 0; i < answer.size(); ++i){

        file << x_0 + h * i << " "; // x
    }
    file<<"\n";

    for(int i = 0; i < answer.size(); ++i){
        file<< answer[i](0)<< " "; // y
    }
    file<<"\n";

    // Вторая ветвь решения
//    std::function<double(double, const VectorXd&)> f2 = [](double x, const VectorXd& y_){
//        return -sqrt(fabs(-y_(0)*exp(y_(1)) + exp(1)/log(x)*y_(0)*y_(0) + 1/(x*x)));
//    };
//
//    std::vector<VectorXd> answer2 = Slae::Solvers::solveShootingMethod(f2, h, x_0, x_f, y_0, y_f, tolerance);
//    for(int i = 0; i < answer2.size(); ++i){
//        file<< answer[i](0)<< " "; // y
//    }
//    file<<"\n";

    file.close();
}
