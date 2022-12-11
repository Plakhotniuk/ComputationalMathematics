//
// Created by Арсений Плахотнюк on 13.04.2022.
//

#ifndef MY_PROJECT_NONLINEARBOUNDARYVALUEPROBLEMMATRIX3_H
#define MY_PROJECT_NONLINEARBOUNDARYVALUEPROBLEMMATRIX3_H
#include "comp_math/utility/Overloads.hpp"
#include "comp_math/sparse/CSR.hpp"
#include <sstream>
#include "comp_math/Exceptions/SlaeBaseException.hpp"
#include <functional>
#include "ThreeDiagonalMatrix.hpp"
#include "cmath"
#include "comp_math/solvers/ThreeDiadonalSolver.hpp"
#include "ostream"

enum COLUMNINDEX {
    C_FIRST = 1,
    C_SECOND = 2,
    C_THIRD = 3,
};

template<int C_Index>
struct Calc{
    static double calc(std::function<double(double, double, double)>& a,
                       std::function<double(double, double, double)>& b, double h, double x, double y, double y_der);
};

template<>
double Calc<COLUMNINDEX::C_FIRST>::calc(std::function<double(double, double, double)>& a,
                                        std::function<double(double, double, double)>& b,
                                        double h, double x, double y, double y_der){
    return 1/(h*h) - a(x, y, y_der)/(2*h);
}

template<>
double Calc<COLUMNINDEX::C_SECOND>::calc(std::function<double(double, double, double)>& a,
                                         std::function<double(double, double, double)>& b,
                                         double h, double x, double y, double y_der){
    return -2/(h*h) + b(x, y, y_der);
}

template<>
double Calc<COLUMNINDEX::C_THIRD>::calc(std::function<double(double, double, double)>& a,
                                        std::function<double(double, double, double)>& b,
                                        double h, double x, double y, double y_der){
    return 1/(h*h) + a(x, y, y_der)/(2*h);
}

std::vector<double> Approach_y_der3(std::vector<double> y_prev, double h){
    std::vector<double> y_der(y_prev.size());
    y_der[0] = (y_prev[1] - y_prev[0]) / h;
    y_der[1] = (y_prev[2] - y_prev[0]) / (2*h);

    for(int i = 2; i < y_prev.size() - 3; ++i)
        y_der[i] = (-y_prev[i+2] + 8*y_prev[i+1] - 8*y_prev[i-1] + y_prev[i-2]) / (12*h);

    y_der[y_prev.size() - 2] = (y_prev.back() - y_prev[y_prev.size() - 3]) / (2*h);
    y_der.back() = (y_prev.back() - y_prev[y_prev.size() - 2]) / h;

    return y_der;
}

std::pair<std::vector<double>, std::vector<double>> InitialApproach_y3(double left_bound_x, double right_bound_x,
                                                                      double left_bound_y, double right_bound_y,
                                                                      unsigned int number_of_splits,
                                                                      std::function<double(double)>& func){
    // шаг разбиения
    auto h = (right_bound_x - left_bound_x) / number_of_splits;

    // начальное приближение
    std::vector<double> y_0(number_of_splits + 1);
    for(unsigned int i = 0; i < y_0.size(); ++i)
    {
//        y_0[i] = func(left_bound_x + h * i);
        y_0[i] = left_bound_y + (left_bound_x + h * i - left_bound_x) /
                                (right_bound_x - left_bound_x) * (right_bound_y - left_bound_y);
    }
    std::vector<double> y_0_der(number_of_splits + 1, (right_bound_y - left_bound_y) / (right_bound_x - left_bound_x));
//    std::vector<double> y_0_der = Approach_y_der(y_0, h);
    return {y_0, y_0_der};
}



std::pair<Slae::Matrix::ThreeDiagonalMatrix, std::vector<double>>
ExpMatrixNonLinearBVP3(double left_bound_x, double right_bound_x,
                       double left_bound_y, double right_bound_y,
                       unsigned int number_of_splits,
                       std::function<double(double, double, double)>& a,
                       std::function<double(double, double, double)>& b,
                       std::function<double(double, double)>& f, std::pair<std::vector<double>, std::vector<double>> y_prev){

    // шаг разбиения
    auto h = (right_bound_x - left_bound_x) / number_of_splits;

    // Матрица коэффициентов
    Slae::Matrix::ThreeDiagonalMatrix data = Slae::Matrix::ThreeDiagonalMatrix(number_of_splits + 1);
    // Коэф в начале
    data(0, 1) = 1;
    //Коэф в конце
    data(data.rows() - 1, 1) = 1;

    //Заполение остальных
    for(unsigned int i = 1; i < data.rows() - 1; ++i){
        data.fill_row(i, Calc<COLUMNINDEX::C_FIRST>::calc(a, b, h, left_bound_x + h * i, y_prev.first[i], y_prev.second[i]),
                      Calc<COLUMNINDEX::C_SECOND>::calc(a, b, h, left_bound_x + h * i, y_prev.first[i], y_prev.second[i]),
                      Calc<COLUMNINDEX::C_THIRD>::calc(a, b, h, left_bound_x + h * i, y_prev.first[i], y_prev.second[i]));
    }

    std::vector<double> y(number_of_splits + 1);
    y[0] = left_bound_y;
    for(unsigned int i = 1; i < y.size() - 1; ++i){
        y[i] = f(left_bound_x + h * i, y_prev.first[i]);
    }
    y.back() = right_bound_y;
    return {data, y};
}

#endif //MY_PROJECT_NONLINEARBOUNDARYVALUEPROBLEMMATRIX3_H
