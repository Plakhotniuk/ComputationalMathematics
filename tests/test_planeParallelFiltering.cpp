//
// Created by Арсений Плахотнюк on 11.12.2022.
//
#include "gtest/gtest.h"
#include "compMath/solvers/ThreeDiadonalSolver.hpp"
#include <iostream>
#include <fstream>
#include <vector>

double getRho(double rho0, double cf, double p, double p0){
    return rho0 * (1. + cf * (p - p0));
}

double getRhoPlus(double p, double pPlus, double rho0, double cf, double p0){
    return p < pPlus ? getRho(rho0, cf, pPlus, p0) : getRho(rho0, cf, p, p0);
}

double getRhoMinus(double p, double pMinus, double rho0, double cf, double p0){
    return pMinus < p ? getRho(rho0, cf, p, p0) : getRho(rho0, cf, pMinus, p0);
}


double getCoefC(double k, double rhoMinus, double mu, double h){
    return k * rhoMinus / (mu * h * h);
}

double getCoefB(double k, double rhoPlus, double mu, double h){
    return k * rhoPlus / (mu * h * h);
}

double getCoefA(double c, double b, double phi, double cf, double rho0, double tau){
    return -c - b - (phi * cf * rho0) / tau;
}

double getCoefD(double phi, double rho0, double cf, double p, double tau){
    return -phi * cf * rho0 * p / tau;
}



TEST(PLANEPARALLELFIKTERING, TASK3){
    double Tmax = 10 * 24 * 60 * 60; // время расчета в секундах
    double Tstop = Tmax;
    double TStep = 60 * 60; // шаг по времени в секундах
    double t = 0;
    uint Nt = std::ceil(Tstop / TStep); // количество узлов по времени

    // Начальные данные
    double dZ = 10; // m
    double L = 500; // m
    double h = 1; // m - шаг сетки по пространству
    uint Nx = std::ceil(L / h); // количество узлов по пространству
    std::vector<double> X(Nx, h);
    double P0 = 100; // atm
    double Pinj = 150; // atm
    double Pprod = 50; // atm
    double phi = 0.2; // пористость пласта
    double rho0 = 1000; // kg/m^3
    double cf = 1.e-4; // atm^(-1)
    double mu = 1.e-3; // Pa * s
    double k = 1.e-14; // m^2


    std::vector<double> d(Nx); // вектор правой части
    std::vector<double> P(Nx, P0); // вектор давлений
    P[0] = Pinj;
    P[P.size() - 1] = Pprod;
    d[0] = Pinj;
    d[d.size() - 1] = Pprod;


    Slae::Matrix::ThreeDiagonalMatrix matrix = Slae::Matrix::ThreeDiagonalMatrix::Zero(Nx); // трех диагональная матрица системы
    matrix.fill_row(0, 0, 1, 0); // фиктивный узел в начале
    matrix.fill_row(Nx-1, 0, 1, 0); // фиктивный узел в конце

    double rhoMinus;
    double rhoPlus;
    double c, b, a;

    while(t < Tstop){
        // Расчет коэффициентов матрицы и вектора правой части
        for(int i = 1; i < Nx - 1; ++i){
            d[i] = getCoefD(phi, rho0, cf, P[i], TStep);
            rhoMinus = getRhoMinus(P[i], P[i-1], rho0, cf, P0);
            rhoPlus = getRhoPlus(P[i], P[i+1], rho0, cf, P0);
            c = getCoefC(k ,rhoMinus, mu, h);
            b = getCoefB(k, rhoPlus, mu, h);
            a = getCoefA(c, b, phi, cf, rho0, TStep);
            matrix.fill_row(i, c, a, b);
        }

        P = Slae::Solvers::solveThreeDiagonal(matrix, d);
        t = t + TStep;
    }

    const std::string FILE_PATH = __FILE__;
    const std::string DIR_PATH = FILE_PATH.substr(0, FILE_PATH.size() - 31);

    std::fstream file;
    file.open(DIR_PATH + "data_files/planeParallelFiltering.txt", std::ios::out);

    for(double p : P){
        file<< p << " "; // y
    }
    file<<std::endl;

    for(double dx : X){
        file<< dx << " "; // y
    }
    file<<std::endl;

    file.close();


}

