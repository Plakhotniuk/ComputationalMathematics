//
// Created by Арсений Плахотнюк on 04.12.2022.
//
#include "gtest/gtest.h"
#include "comp_math/matrix/LinearBoundaryValueProblemMatrix3.hpp"
#include "comp_math/utility/Overloads.hpp"
#include "comp_math/interpolation/Linear.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
// Задача 3


TEST(LINEARBOUNDARYVALUEPROBLEM1, TEST1) {
    std::function<double(double)> f = [](double x) { return 2 - 6*x + 2*pow(x, 3) + (x*x - 3)*exp(x)*sin(x)*(1+cos(x)) +
                                                     cos(x)*(exp(x) + (x*x - 1) + pow(x, 4) - 3*x*x); };

    std::function<double(double)> a = [](double x) { return x*x - 3; };

    std::function<double(double)> b = [](double x) { return (x*x - 3)*cos(x); };


    double left_bound_x = 0.;
    double right_bound_x = M_PI;
    double left_bound_y = 0.;
    double right_bound_y = M_PI*M_PI;

    int max_number_of_splits = 100;

    std::vector<double> x(max_number_of_splits + 1);
    double h = (right_bound_x - left_bound_x)/max_number_of_splits;
    for(int i = 0; i < x.size(); ++i){
        x[i] = left_bound_x + h * i;
    }


    std::pair<Slae::Matrix::ThreeDiagonalMatrix, std::vector<double>> matrix =
            ExpandedMatrixForLinearBoundaryValueProblem3(left_bound_x, right_bound_x,
                                                         left_bound_y, right_bound_y,
                                                         max_number_of_splits, a, b, f);

    std::vector<double> solution = Slae::Solvers::solveThreeDiagonal(matrix.first, matrix.second);


    std::vector<double> answer_x = {0.5, 1, 1.5, 2, 2.5, 3};
    for(double i : answer_x){
        std::cout<< linear_interpolation(x, solution, i) << std::endl;
    }

    std::fstream file;
    file.open("test_boundary_problem.txt", std::fstream::out);

    for(double i : solution){
        file<< i<< " ";

    }
    file<<"\n";

    for(double i : x){
        file<< i<< " ";

    }
    file<<"\n";
    file.close();

}