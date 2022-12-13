//
// Created by Арсений Плахотнюк on 11.12.2022.
//
#include "gtest/gtest.h"
#include "comp_math/types/BasicTypes.hpp"
#include "comp_math/solvers/ShootingMethod.hpp"
#include "comp_math/matrix/NonLinearBoundaryValueProblemMatrix3.hpp"
#include "comp_math/interpolation/Linear.hpp"
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
    std::vector<double> x(num_of_steps + 1);
    for(int i = 0; i < x.size(); ++i){
        x[i] = x_0 + h * i;
    }

    std::vector<VectorXd> answer = Slae::Solvers::solveShootingMethod(f, num_of_steps, x_0, x_f, y_0, y_f, tolerance);
//    std::vector<double> solution(answer.size());
//    for(int i = 0; i < answer.size(); ++i){
//        solution[i] = answer[i](0);
//    }
//
//    std::vector<double> answer_x = {3.5, 4, 4.5, 5, 5.5};
//    for(double i : answer_x){
//        std::cout<< linear_interpolation(x, solution, i) << std::endl;
//    }

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


TEST(NONLINEARBVP, NONLINEARBVP) {
    std::function<double(double, double)> f = [](double x, double y) { return 0.; };

    std::function<double(double, double, double)> a = [](double x, double y, double y_der) {
        return 0.;
    };

    std::function<double(double, double, double)> b = [](double x, double y, double y_der) {
        return -sqrt(abs(-y*exp(y_der) + exp(1) / log(x) * y * y + 1 / (x * x)));
    };

    std::function<double(double)> y_func = [](double x) { return 2 * x; };

    double left_bound_x = exp(1);
    double right_bound_x = exp(1)*exp(1);
    double left_bound_y = exp(1);
    double right_bound_y = 2. * exp(1)*exp(1);

    int n_iterations = 7;
    int max_number_of_splits = 1000;
    std::vector<double> y(max_number_of_splits + 1);
    double h = (right_bound_x - left_bound_x)/max_number_of_splits;

    std::fstream file;
    file.open("test_nonlin_bvp.txt", std::fstream::out);

    std::pair<std::vector<double>, std::vector<double>> solution;

    std::pair<std::vector<double>, std::vector<double>> y_0_y_0_der = InitialApproach_y3(left_bound_x, right_bound_x,
                                                                                         left_bound_y, right_bound_y,
                                                                                         max_number_of_splits, y_func);

//approach
    solution = y_0_y_0_der;

    for(int i = 1; i < n_iterations; ++i){
        std::pair<Slae::Matrix::ThreeDiagonalMatrix, std::vector<double>> matrix =
                ExpMatrixNonLinearBVP3(left_bound_x, right_bound_x,
                                       left_bound_y, right_bound_y,
                                       max_number_of_splits, a, b, f, solution);
        solution.first = Slae::Solvers::solveThreeDiagonal(matrix.first, matrix.second);
    }
    for(int i = 0; i < solution.first.size(); ++i){
        file<< solution.first[i]<< " ";

    }
    file<<"\n";
    for(int i = 0; i < solution.first.size(); ++i){

        file << left_bound_x + h * i << " ";
    }
    file<<"\n";


    file.close();
}
