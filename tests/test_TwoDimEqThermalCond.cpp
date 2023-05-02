//
// Created by Арсений Плахотнюк on 01.05.2023.
//

#include "gtest/gtest.h"
#include "compMath/solvers/ThreeDiadonalSolver.hpp"
#include "compMath/types/BasicTypes.hpp"
#include "compMath/utility/Norm.hpp"
#include "compMath/utility/Overloads.hpp"
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

double analyticalSolution(double x, double y, double t, double lambda_ = 1.e-4){
    return std::cos(M_PI * x) * std::sin(5 * M_PI * y) * std::exp(-50*M_PI*M_PI*lambda_*t);
}

std::vector<double> getLinspace(double leftBound, double rightBound, uint N){
    std::vector<double> result(N, 0.);
    double h = (rightBound - leftBound) / (N-1);
    for(int i = 0; i < N; ++i){
        result[i] = leftBound + i*h;
    }
    return result;
}

double coefA(double tau, double hx, double hy, double lambda_ = 1.e-4){
    return 1. + 2. * tau * lambda_ * (25. / (hx*hx) + 1. / (hy*hy));
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

    const uint NX = 150;
    const uint NY = NX;
    std::vector<double> x = getLinspace(xLeftBound, xRightBound, NX);
    double hx = (xRightBound - xLeftBound) / (NX-1);
    std::vector<double> y = getLinspace(yLeftBound, yRightBound, NY);
    double hy = (yRightBound - yLeftBound) / (NY-1);



    const double lambda_ = 1.e-4;


    Slae::Matrix::ThreeDiagonalMatrix matrix = Slae::Matrix::ThreeDiagonalMatrix::Zero(NY);

    std::vector<double> a(NX, coefA(tStep, hx, hy));
    std::vector<double> b(NX, coefB(tStep, hx));
    std::vector<double> c(NX, coefC(tStep, hx));
    std::vector<double> d(NX, 0.); // нулевые начальные условия
    std::vector<double> f(NX, coefF(tStep, hy));
    std::vector<double> g(NX, coefG(tStep, hy));

    std::vector<double> dTilda(NX, 0.);

    // Генерация сетки и начальные условия
    std::vector<std::vector<double>> solution(NY, std::vector<double>(NX));
    for(int j = 0; j < NY; ++j){
        for(int i = 0; i < NX; ++i){
            solution[j][i] = initCondition(x[i], y[j]);
        }
    }

    int inv = 1;
    while(t < tEnd){

        // Граничные условия
        for(int i = 0; i < NX; ++i){
            solution[0][i] = leftBoundCondY(x[i], t);
            solution[NY-1][i] = rightBoundCondY(x[i], t);
            solution[i][0] = leftBoundCondX(y[i], t);
            solution[i][NX-1] = rightBoundCondX(y[i], t);
        }

        if(inv % 2 == 1){
            for(int j = 1; j < NY - 1; ++j){
                d[0] = leftBoundCondY(x[j], t);
                dTilda[0] = d[0];
                matrix.fill_row(0, 0, 1, 0);
                for(int i = 1; i < NX - 1; ++i){
                    d[i] = solution[j][i];
                    dTilda[i] = d[i] - g[i] * solution[j-1][i] - f[i] * solution[j+1][i];
                    matrix.fill_row(i, c[i], a[i], b[i]);
                }
                d[NY - 1] = rightBoundCondY(x[j], t);
                dTilda[NY - 1] = d[NY - 1];
                matrix.fill_row(NY - 1, 0, 1, 0);
                solution[j] = Slae::Solvers::solveThreeDiagonal(matrix, dTilda);
            }
        }
        else{
            for(int j = 1; j < NY - 1; ++j){
                d[0] = leftBoundCondX(y[j], t);
                dTilda[0] = d[0];
                matrix.fill_row(0, 0, 1, 0);
                for(int i = 1; i < NX - 1; ++i){
                    d[i] = solution[i][j];
                    dTilda[i] = d[i] - c[i] * solution[j-1][i] - b[i] * solution[i][j+1];
                    matrix.fill_row(i, g[i], a[i], f[i]);
                }
                d[NX - 1] = rightBoundCondX(y[j], t);
                dTilda[NY - 1] = d[NY - 1];
                matrix.fill_row(NY - 1, 0, 1, 0);
                solution[j] = Slae::Solvers::solveThreeDiagonal(matrix, dTilda);
            }
        }
        t += tStep;
        ++inv;
    }

    std::vector<std::vector<double>> anSol(NY, std::vector<double>(NX));
    for(int j = 0; j < NY; ++j){
        for(int i = 0; i < NX; ++i){
            anSol[j][i] = analyticalSolution(x[i], y[j], t);
        }
    }

    std::vector<double> diff(NY);
    for(int i = 0; i < NY; ++i){
        diff[i] = Norm<double, NormType::FirstNorm>::get_norm(anSol[i] - solution[i]);
    }

    std::cout << Norm<double, NormType::FirstNorm>::get_norm(diff) << std::endl;

}
