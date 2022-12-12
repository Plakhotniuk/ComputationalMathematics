//
// Created by Арсений Плахотнюк on 11.12.2022.
//
#include "gtest/gtest.h"
#include "comp_math/integrators/DormandPrince.hpp"
#include "comp_math/integrators/RungeKutta4.hpp"
#include "comp_math/types/BasicTypes.hpp"
#include <iostream>
#include <fstream>

TEST(THREEBODIES, SOLVEDP7){
// Задача 5
// y_ = [x, u, y, nu]
    std::function<VectorXd(double, const VectorXd&)> f = [](double t, const VectorXd& y_){
        VectorXd res(4);
        double mu = 0.012277471;
        double eta = 1 - mu;
        double x = y_(0);
        double u = y_(1);
        double y = y_(2);
        double nu = y_(3);
        double A = pow((x + mu)*(x + mu) + y*y, 1.5);
        double B = pow((x - eta)*(x - eta) + y*y, 1.5);
        res << u,
        x + 2*nu - eta*((x + mu) / A) - mu*((x - eta) / B),
        nu,
        y - 2*u - eta*(y/A) - mu*(y/B);

        return res;
    };
    double tolerance = 1.e-12;
    double T = 17.0652165601579625588917206249;
    double t_0 = 0.;
    double t_f = 5*T;
    const int n = 100000;
    double h = (t_f - t_0)/n;
    VectorXd u_0(4);
    u_0 << 0.994, 0., 0., -2.00158510637908252240537862224;

    std::pair<std::vector<double>, std::vector<VectorXd>> answer = DormandPrince(f, t_0, t_f, u_0, h, tolerance);

    std::fstream file;
    file.open("test_three_bodies_dp7.txt", std::fstream::out);


    for(int i = 0; i < answer.first.size(); i+=100){
    // x
    file<< answer.second[i](0)<< " ";

    }
    file<<"\n";

    for(int i = 0; i < answer.first.size(); i+=100){
    // y
    file<< answer.second[i](2)<< " ";

    }
    file<<"\n";


    file.close();

}

