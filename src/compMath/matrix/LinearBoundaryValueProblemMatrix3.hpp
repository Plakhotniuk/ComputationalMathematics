//
// Created by Арсений Плахотнюк on 11.03.2022.
//

#ifndef MY_PROJECT_LINEARBOUNDARYVALUEPROBLEMMATRIX3_HPP
#define MY_PROJECT_LINEARBOUNDARYVALUEPROBLEMMATRIX3_HPP
#include "compMath/utility/Overloads.hpp"
#include "compMath/sparse/CSR.hpp"
#include <sstream>
#include "compMath/Exceptions/SlaeBaseException.hpp"
#include <functional>
#include "ThreeDiagonalMatrix.hpp"
#include "cmath"
#include "compMath/solvers/ThreeDiadonalSolver.hpp"
#include "ostream"

enum COLUMNINDEX {
    C_FIRST = 1,
    C_SECOND = 2,
    C_THIRD = 3,
};

template<int C_Index>
struct Calc{
    static double calc(std::function<double(double)>& a,
                       std::function<double(double)>& b, double h, double x);
};

template<>
double Calc<COLUMNINDEX::C_FIRST>::calc(std::function<double(double)>& a,
                                                                std::function<double(double)>& b, double h, double x){
    return 1/(h*h) - a(x)/(2*h);
}

template<>
double Calc<COLUMNINDEX::C_SECOND>::calc(std::function<double(double)>& a,
                                                                 std::function<double(double)>& b, double h, double x){
    return -2/(h*h) + b(x);
}

template<>
double Calc<COLUMNINDEX::C_THIRD>::calc(std::function<double(double)>& a,
                                                                std::function<double(double)>& b, double h, double x){
    return 1/(h*h) + a(x)/(2*h);
}

std::pair<Slae::Matrix::ThreeDiagonalMatrix, std::vector<double>> ExpandedMatrixForLinearBoundaryValueProblem3(double left_bound_x, double right_bound_x, double left_bound_y,
                                                                                                               double right_bound_y, int number_of_splits,
                                                                                                               std::function<double(double)>& a, std::function<double(double)>& b,
                                                                                                               std::function<double(double)>& f) {
    // шаг разбиения
    auto h = (right_bound_x - left_bound_x) / number_of_splits;

    // Матрица коэффициентов
    Slae::Matrix::ThreeDiagonalMatrix data = Slae::Matrix::ThreeDiagonalMatrix(number_of_splits + 1);

    // Коэф в начале
    data(0, 1) = 1;
    //Коэф в конце
    data(data.rows() - 1, 1) = 1;
    //Заполение остальных
    for(int i = 1; i < data.rows() - 1; ++i){
            data.fill_row(i, Calc<COLUMNINDEX::C_FIRST>::calc(a, b, h, left_bound_x + h * i),
                          Calc<COLUMNINDEX::C_SECOND>::calc(a, b, h, left_bound_x + h * i),
                          Calc<COLUMNINDEX::C_THIRD>::calc(a, b, h, left_bound_x + h * i));
    }
    std::vector<double> y(number_of_splits + 1);
    y[0] = left_bound_y;
    for(int i = 1; i < y.size() - 1; ++i){
        y[i] = f(left_bound_x + h * i);
    }
    y.back() = right_bound_y;
    return {data, y};
}
#endif //MY_PROJECT_LINEARBOUNDARYVALUEPROBLEMMATRIX3_HPP
