//
// Created by Арсений Плахотнюк on 11.12.2022.
//
#include "gtest/gtest.h"
#include "comp_math/types/BasicTypes.hpp"
#include "comp_math/solvers/ShootingMethod.hpp"
#include "comp_math/matrix/NonLinearBoundaryValueProblemMatrix3.hpp"
#include <iostream>
#include <fstream>

TEST(SHOOTING, METHOD){
    //XI.9.3 a

    // скалярная функция правой части уравнения
    std::function<double(double, const VectorXd&)> f = [](double x, const VectorXd& y_){
        return x * sqrt(y_(0));
    };
    int num_of_steps = 1000;

    double x_0 = 0; // граничные условия
    double y_0 = 0;
    double x_f = 1;
    double y_f = 2;
    double tolerance = 1.e-10;
    std::vector<VectorXd> answer = Slae::Solvers::solveShootingMethod(f, num_of_steps, x_0, x_f, y_0, y_f, tolerance);
    std::fstream file;
    file.open("test_shooting_method.txt", std::fstream::out);
    double h = (x_f - x_0) / num_of_steps;
    for(int i = 0; i < answer.size(); ++i){

        file << x_0 + h * i << " "; // x
    }
    file<<"\n";

    for(int i = 0; i < answer.size(); ++i){
        file<< answer[i](0)<< " "; // y
    }
    file<<"\n";

    file.close();
}

TEST(THREE_DIAG, CHECK){
    std::function<double(double, double)> f = [](double x, double y) { return 0.; };

    std::function<double(double, double, double)> a = [](double x, double y, double y_dev) { return 0.; };

    std::function<double(double, double, double)> b = [](double x, double y, double y_dev) { return -x/ sqrt(y); };

    std::function<double(double)> y_func = [](double x) { return 2 * x; };

    double left_bound_x = 0.;
    double right_bound_x = 1.;
    double left_bound_y = 0.;
    double right_bound_y = 2.;

    int n_iterations = 7;
    int max_number_of_splits = 100;
    std::vector<double> y(max_number_of_splits + 1);
    double h = (right_bound_x - left_bound_x)/max_number_of_splits;

    std::fstream file;
    file.open("test_check_shooting_method.txt", std::fstream::out);


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

        file << left_bound_x + h * i << " "; // x
    }
    file<<"\n";

    for(int i = 0; i < solution.first.size(); ++i){
        file<< solution.first[i]<< " "; // y

    }
    file<<"\n";

    file.close();
}
