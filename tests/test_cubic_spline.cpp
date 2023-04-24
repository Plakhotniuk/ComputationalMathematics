//
// Created by Арсений Плахотнюк on 27.09.2022.
//
#include "gtest/gtest.h"
#include "compMath/solvers/ThreeDiadonalSolver.hpp"
#include <iostream>
#include "functional"
#include <fstream>
#include <vector>

double get_dif(const std::vector<double> &x, const std::vector<double> &f)
{
    if (f.size() > 2)
    {
        std::vector<double> f_left(f.size() - 1);
        std::vector<double> x_left(x.size() - 1);
        std::vector<double> f_right(f.size() - 1);
        std::vector<double> x_right(x.size() - 1);
        for(int i = 0; i < f.size() - 1; ++i){
            f_left[i] = f[i];
            x_left[i] = x[i];
            f_right[i] = f[i+1];
            x_right[i] = x[i+1];
        }
        return (get_dif(x_left, f_left) - get_dif(x_right, f_right)) / (x[x.size()-1] - x[0]);
    }
    else if (f.size() == 2) {
        return (f[1] - f[0]) / (x[1] - x[0]);
    }
}

TEST(THREECUBICSPLINE, TEST) {

    std::vector<double> x = {1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000};
    std::vector<double> y = {92228496, 106021537, 123202624, 132164569, 151325798, 179323175, 203211926, 226545805,
                             248709873, 281421906};
    int n = x.size();
    std::vector<double> h(n - 1);

    for (int i = 0; i < n - 1; ++i) {
        h[i] = x[i + 1] - x[i];
    }

    Slae::Matrix::ThreeDiagonalMatrix matrix = Slae::Matrix::ThreeDiagonalMatrix::Zero(n - 2);
    std::vector<double> u(n - 2);
    for (int i = 1; i < n - 1; ++i) {
        u[i - 1] = 6 * get_dif(std::vector({x[i - 1], x[i], x[i + 1]}), std::vector({y[i - 1], y[i], y[i + 1]}));
    }
    matrix.fill_row(0, 0, 2, h[1] / (h[0] + h[1]));

    for (int i = 1; i < n - 3; ++i) {
        matrix.fill_row(i,
                        h[i - 1] / (h[i - 1] + h[i]),
                        2,
                        h[i - 1] / (h[i - 1] + h[i]));
    }
    matrix.fill_row(n - 3, h[n - 3] / (h[n - 3] + h[n - 2]), 2, 0);

    std::vector<double> c = Slae::Solvers::solveThreeDiagonal(matrix, u);
    c.insert(c.begin(), 0);
    c.emplace_back(0);

    std::vector<double> a = y;
//    a.erase(a.begin());

    std::vector<double> b(n - 1);

    for (int i = 0; i < n; ++i) {
        b[i] = (c[i + 1] * h[i + 1] / 3) + (c[i] * h[i + 1] / 6) + get_dif({x[i], x[i + 1]}, {y[i], y[i + 1]});
    }

    std::vector<double> d(n - 1);
    d[0] = c[1] / h[0];
    for (int i = 1; i < n; ++i) {
        d[i] = (c[i + 1] - c[i]) / h[i];
    }
    double var = 2010;
    int i = n-2;
    for(int j = 0; j < n - 2; ++j){
        if(x[j] <= var && var < x[j+1]){
            i = j;
        }
    }

    double res = a[i] + b[i] * (var - x[i]) +
                 c[i] * (var - x[i]) * (var - x[i]) / 2 + d[i] * (var - x[i]) * (var - x[i]) * (var - x[i]) / 6;
    std::cout << std::setprecision(15) << res << std::endl;


}

TEST(THREECUBICSPLINE_LABA, TEST) {

    std::vector<double> x = {2134, 1866, 1828, 1816, 2150, 2184, 2196, 2224, 2244, 2256, 2582, 2552, 2488, 2474, 2450, 2424, 2414, 2374, 2366, 2350, 2338, 2322, 2302, 2280, 2270, 2254, 2548, 2316, 2102, 2092, 1918, 1478, 808, 250};
    std::vector<double> y = {5852, 5401, 5341, 5331, 5882, 5945, 5976, 6030, 6074, 6096, 7032, 6929, 6717, 6678, 6599, 6533, 6507, 6402, 6383, 6334, 6305, 6267, 6217, 6164, 6143, 6096, 6907, 6234, 5791, 5770, 5461, 4916, 4358, 4047};
    int n = x.size();
    std::vector<double> h(n - 1);

    std::sort(x.begin(), x.end());
    std::sort(y.begin(), y.end());

    for (int i = 0; i < n - 1; ++i) {
        h[i] = x[i + 1] - x[i];
    }

    Slae::Matrix::ThreeDiagonalMatrix matrix = Slae::Matrix::ThreeDiagonalMatrix::Zero(n - 2);
    std::vector<double> u(n - 2);
    for (int i = 1; i < n - 1; ++i) {
        u[i - 1] = 6 * get_dif(std::vector({x[i - 1], x[i], x[i + 1]}), std::vector({y[i - 1], y[i], y[i + 1]}));
    }
    matrix.fill_row(0, 0, 2, h[1] / (h[0] + h[1]));

    for (int i = 1; i < n - 3; ++i) {
        matrix.fill_row(i,
                        h[i - 1] / (h[i - 1] + h[i]),
                        2,
                        h[i - 1] / (h[i - 1] + h[i]));
    }
    matrix.fill_row(n - 3, h[n - 3] / (h[n - 3] + h[n - 2]), 2, 0);

    std::vector<double> c = Slae::Solvers::solveThreeDiagonal(matrix, u);
    c.insert(c.begin(), 0);
    c.emplace_back(0);

    std::vector<double> a = y;
//    a.erase(a.begin());

    std::vector<double> b(n - 1);

    for (int i = 0; i < n; ++i) {
        b[i] = (c[i + 1] * h[i + 1] / 3) + (c[i] * h[i + 1] / 6) + get_dif({x[i], x[i + 1]}, {y[i], y[i + 1]});
    }

    std::vector<double> d(n - 1);
    d[0] = c[1] / h[0];
    for (int i = 1; i < n; ++i) {
        d[i] = (c[i + 1] - c[i]) / h[i];
    }
    double var = 782;
    int i = n-2;
    for(int j = 0; j < n - 2; ++j){
        if(x[j] <= var && var < x[j+1]){
            i = j;
        }
    }

    double res = a[i] + b[i] * (var - x[i]) +
                 c[i] * (var - x[i]) * (var - x[i]) / 2 + d[i] * (var - x[i]) * (var - x[i]) * (var - x[i]) / 6;
    std::cout << std::setprecision(15) << res << std::endl;


}



