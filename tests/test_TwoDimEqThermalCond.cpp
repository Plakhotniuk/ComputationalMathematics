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

    const std::string FILE_PATH = __FILE__;
    const std::string DIR_PATH = FILE_PATH.substr(0, FILE_PATH.size() - 28);

    const double xLeftBound = 0;
    const double xRightBound = 1;
    const double yLeftBound = 0;
    const double yRightBound = 1;
    const double tStart = 0;
    const double tEnd = 1;
    const double tStep = 0.0001;
    double t = tStart + tStep;

    const uint NX = 100;
    const uint NY = NX;
    std::vector<double> x = getLinspace(xLeftBound, xRightBound, NX);
    double hx = (xRightBound - xLeftBound) / (NX-1);
    std::vector<double> y = getLinspace(yLeftBound, yRightBound, NY);
    double hy = (yRightBound - yLeftBound) / (NY-1);


    const double lambda_ = 1.e-4;


    Slae::Matrix::ThreeDiagonalMatrix matrix = Slae::Matrix::ThreeDiagonalMatrix::Zero(NY);

    const double a =  coefA(tStep, hx, hy);
    const double b = coefB(tStep, hx);
    const double c = coefC(tStep, hx);
    const double f = coefF(tStep, hy);
    const double g = coefG(tStep, hy);

    std::vector<double> dTilda(NX, 0.);

    std::fstream file;
    file.open(DIR_PATH + "data_files/twoDimEqThermalCond.txt", std::ios::out);

    // Генерация сетки и начальные условия
    std::vector<std::vector<double>> solution(NY, std::vector<double>(NX));
    for(int j = 0; j < NY; ++j){
        for(int i = 0; i < NX; ++i){
            solution[j][i] = initCondition(x[i], y[j]);
            file << solution[j][i] << " ";
        }
        file<<"\n";
    }

    std::vector<double> buff(NY-2, 0.);


    //region Цикл по времени

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
            for(int k = 1; k < NY - 1; ++k){
//                dTilda[0] = solution[k][0];
//                matrix.fill_row(0, 0, 1, 0);

                dTilda[1] = solution[k][1] - g * solution[k-1][1] - f * solution[k+1][1] - c * solution[k][0];
                matrix.fill_row(1, 0, a, b);

                for(int m = 2; m < NX - 2; ++m){
                    dTilda[m] = solution[k][m] - g * solution[k - 1][m] - f * solution[k + 1][m];
                    matrix.fill_row(m, c, a, b);
                }
                dTilda[NX - 2] = solution[k][NX - 2] - g * solution[k-1][NX - 2] - f * solution[k+1][NX - 2] - b*solution[k][NX-1];
                matrix.fill_row(NY - 2, c, a, 0);

                dTilda[NX - 1] = solution[k][NX - 1];
                matrix.fill_row(NY - 1, 0, 1, 0);

                solution[k] = Slae::Solvers::solveThreeDiagonal(matrix, dTilda);
            }
        }
        else{
            for(int m = 1; m < NX - 1; ++m){
                dTilda[0] = solution[0][m];
                matrix.fill_row(0, 0, 1, 0);

                dTilda[1] = solution[1][m] - c * solution[1][m-1] - b * solution[1][m+1] - g * solution[0][m];
                matrix.fill_row(0, 0, a, f);

                for(int k = 2; k < NY - 2; ++k){
                    dTilda[k] = solution[m][k] - c * solution[m][k - 1] - b * solution[m][k + 1];
                    matrix.fill_row(k, g, a, f);
                }
                dTilda[NY - 2] = solution[NY - 2][m] - c * solution[NY - 2][m-1] - b * solution[NY - 2][m+1] - f*solution[NY-1][m];
                matrix.fill_row(NY - 2, 0, g, a);

                dTilda[NY - 1] = solution[NY - 1][m];
                matrix.fill_row(NY - 1, 0, 1, 0);

                buff = Slae::Solvers::solveThreeDiagonal(matrix, dTilda);

                for(int k = 0; k < NY; ++k){
                    solution[k][m] = buff[k];
                }
            }
        }

        if(inv % 100 == 0) {
            for (int j = 0; j < NY; ++j) {
                for (int i = 0; i < NX; ++i) {
                    solution[j][i] = initCondition(x[i], y[j]);
                    file << solution[j][i] << " ";
                }
                file << "\n";
            }
        }

        t += tStep;
        ++inv;
    }

//    file.close();



    //region Аналитическое решение
    std::vector<std::vector<double>> anSol(NY, std::vector<double>(NX));
    for(int j = 0; j < NY; ++j){
        for(int i = 0; i < NX; ++i){
            anSol[j][i] = analyticalSolution(x[i], y[j], t);
        }
    }
    //endregion


    //region Сравнение по норме
    std::vector<double> diff(NY);
    for(int i = 0; i < NY; ++i){
        diff[i] = Norm<double, NormType::FirstNorm>::get_norm(anSol[i] - solution[i]);
    }

    std::cout << Norm<double, NormType::FirstNorm>::get_norm(diff) << std::endl;

    //endregion

}
