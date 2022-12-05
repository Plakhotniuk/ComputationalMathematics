//
// Created by Арсений Плахотнюк on 05.12.2022.
//
#include "gtest/gtest.h"
#include "comp_math/integrators/DormandPrince.hpp"
#include "comp_math/integrators/RungeKutta4.hpp"
#include "comp_math/types/BasicTypes.hpp"
#include <iostream>
#include <fstream>

TEST(RAYLEIGH1, RK4){
    // X.9.3
    // y_ = [y, y']
    // p = y'
    std::function<double(double, double)> f = [](double x, double p){ return 1000*(1. - p*p)*p - x; };
    double t_0 = 0.;
    double t_f = 1000.;
    const int n = 1000000;
    const double p_0 = 0.001;
    const double x_0 = 0.;

    std::vector<double> p = RungeKutta4(f, t_0, t_f, p_0, n);

    std::vector<double> t(n);
    double h = (t_f - t_0)/n;
    for(int i = 0; i < n; ++i){
        t[i] = t_0 + i*h;
    }


    std::function<double(double, double)> g = [t_0, h, p](double x, double y){ return p[static_cast<int>((x - t_0) / h)]; };

    std::vector<double> answer = RungeKutta4(g, t_0, t_f, x_0, n);

    std::fstream file;
    file.open("test_rayleigh_rk.txt", std::fstream::out);

    for(double i : answer){
        file<< i<< " ";

    }
    file<<"\n";

    for(int i = 0; i < answer.size(); ++i){

        file << t_0 + h * i << " ";
    }
    file<<"\n";
    file.close();

}


TEST(RAYLEIGH2, RK4V){
    // X.9.3
    // y_ = [y, y']
    std::function<VectorXd(double, const VectorXd&)> f = [](double x, const VectorXd& y_){
        VectorXd res(2);
        res << y_[1], 1000*(1. - y_[1]*y_[1])*y_[1] - x;
        return res;
    };
    double t_0 = 0.;
    double t_f = 1000.;
    const int n = 1000000;
    double h = (t_f - t_0)/n;
    VectorXd u_0(2);
    u_0 << 0., 0.001;

    std::vector<VectorXd> answer = RungeKutta4Vec(f, t_0, t_f, u_0, n);

    std::fstream file;
    file.open("test_rayleigh_rk.txt", std::fstream::out);

    for(int i = 0; i < answer.size(); i+=1000){
        file<< answer[i](0)<< " ";

    }
    file<<"\n";

    for(int i = 0; i < answer.size(); i+=1000){

        file << t_0 + h * i << " ";
    }
    file<<"\n";
    file.close();

}


TEST(RAYLEIGH3, DP7){
    // X.9.3
    // y_ = [y, y']
    std::function<VectorXd(double, const VectorXd&)> f = [](double x, const VectorXd& y_){
        VectorXd res(2);
        res << y_[1], 1000*(1. - y_[1]*y_[1])*y_[1] - x;
        return res;
    };
    double t_0 = 0.;
    double t_f = 1000.;
    const int n = 10000;
    double h = (t_f - t_0)/n;
    VectorXd u_0(2);
    u_0 << 0., 0.001;

    std::pair<std::vector<double>, std::vector<VectorXd>> answer = DormandPrince(f, t_0, t_f, u_0, h);

    std::fstream file;
    file.open("test_rayleigh_dp.txt", std::fstream::out);

    for(int i = 0; i < answer.first.size(); i+=1000){
        file<< answer.second[i](0)<< " ";

    }
    file<<"\n";

    for(int i = 0; i < answer.first.size(); i+=1000){
        file << answer.first[i] << " ";
    }
    file<<"\n";
    file.close();

}
