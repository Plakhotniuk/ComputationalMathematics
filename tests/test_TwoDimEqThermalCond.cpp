//
// Created by Арсений Плахотнюк on 01.05.2023.
//

#include "gtest/gtest.h"
#include "compMath/solvers/ThreeDiadonalSolver.hpp"
#include "compMath/types/BasicTypes.hpp"
#include <iostream>
#include <fstream>
#include <vector>

double initCondition(double x, double y){
    /**
     * Начальное условие
     */
    return std::cos(M_PI * x) * std::sin(5 * M_PI * y);
}

double leftBoundCondX(double y, double t, double lambda_ = 1.e-4){
    return std::sin(5 * M_PI * y) * std::exp(-50*M_PI*M_PI*lambda_*t);
}

double rightBoundCondX(double y, double t, double lambda_ = 1.e-4){
    return -std::sin(5 * M_PI * y) * std::exp(-50*M_PI*M_PI*lambda_*t);
}

double leftBoundCondY(double x, double t, double lambda_ = 1.e-4){
    return 0;
}

double rightBoundCondY(double x, double t, double lambda_ = 1.e-4){
    return 0;
}

double analyticalSolution(double x, double y, double t, double lambda_){
    return std::cos(M_PI * x) * std::sin(5 * M_PI * y) * std::exp(-50*M_PI*M_PI*lambda_*t);
}

VectorXd getLinspace(double leftBound, double rightBound, uint N){
    VectorXd result = VectorXd::Zero(N);
    double h = (rightBound - leftBound) / N;
    for(int i = 0; i < N; ++i){
        result(i) = i*h;
    }
    return result;
}

double coefA(double tau, double hx, double hy, double lambda_ = 1.e-4){
    return 1. - 2. * tau * lambda_ * (25. / (hx*hx) + 1. / (hy*hy));
}

double coefB(double tau, double hx, double lambda_ = 1.e-4){
    return -tau * 25. * lambda_ / (hx*hx);
}

double coefF(double tau, double hy, double lambda_ = 1.e-4){
    return -tau * lambda_ / (hy*hy);
}

double coefC(double tau, double hx, double lambda_ = 1.e-4){
    return -tau * 25. * lambda_ / (hx*hx);
}

double coefG(double tau, double hy, double lambda_ = 1.e-4){
    return -tau * lambda_ / (hy*hy);
}


TEST(THERMAL, CONDUCTIVITY){

    const double xLeftBound = 0;
    const double xRightBound = 1;
    const double yLeftBound = 0;
    const double yRightBound = 1;
    const double tStart = 0;
    const double tEnd = 1;
    const double tStep = 0.001;
    double t = tStart;

    const uint NX = 10;
    const uint NY = NX;
    VectorXd x = getLinspace(xLeftBound, xRightBound, NX);
    double hx = (xRightBound - xLeftBound) / NX;
    VectorXd y = getLinspace(yLeftBound, yRightBound, NY);
    double hy = (yRightBound - yLeftBound) / NY;



    const double lambda_ = 1.e-4;


    Slae::Matrix::ThreeDiagonalMatrix matrix = Slae::Matrix::ThreeDiagonalMatrix::Zero(NX);
    matrix.fill_row(0, 0, 1, 0);
    matrix.fill_row(NX - 1, 0, 1, 0);

    std::vector<double> a(NX, coefA(tStep, hx, hy));
    std::vector<double> b(NX, coefB(tStep, hx));
    std::vector<double> c(NX, coefC(tStep, hx));
    std::vector<double> d(NX, 0.);
    std::vector<double> f(NX, coefF(tStep, hy));
    std::vector<double> g(NX, coefG(tStep, hy));

    std::vector<double> dTilda(NX, 0.); // нулевые начальные условия

    // Генерация сетки и начальные условия
    std::vector<std::vector<double>> solution(NY, std::vector<double>(NX));

    for(int j = 0; j < NY; ++j){
        for(int i = 0; i < NX; ++i){
            solution[j][i] = initCondition(x(i), y(j));
        }
    }


    while(t < tEnd){
        for(int j = 1; j < NY - 1; ++j){
            d[0] = solution[j][0];
            d[d.size() - 1] = solution[j][NX - 1];
            for(int i = 1; i < NX - 1; ++i){
                d[i] = solution[j][i];
                dTilda[i] = d[i] - g[i] * solution[j][i-1] - f[i] * solution[j][i+1];
                matrix.fill_row(i, c[i], a[i], b[i]);
            }
            solution[j] = Slae::Solvers::solveThreeDiagonal(matrix, d);

        }

        t += tStep;
    }


}