//
// Created by Арсений Плахотнюк on 04.12.2022.
//

#ifndef COMPUTATIONALMATHEMATICS_DORMANDPRINCE_HPP
#define COMPUTATIONALMATHEMATICS_DORMANDPRINCE_HPP
#include <vector>
#include <cmath>
#include "comp_math//types/BasicTypes.hpp"

std::pair<std::vector<double>, std::vector<VectorXd>> DormandPrince(std::function<VectorXd(double, const VectorXd&)>& f,
                                                                  const double x_0, const double x_f, const VectorXd& y_0,
                                                                  double h, const double tolerance = 1.e-7){
    VectorXd k1, k2, k3, k4, k5, k6, k7;
    std::vector<VectorXd> y = {y_0};
    std::vector<double> x = {x_0};
    VectorXd y1; // O(h^4) approximation
    VectorXd y2; // O(h^5) approximation
    double s;
    int num_of_steps = 0;
    double x_n = x_0;
    while(x_n < x_f){
        x_n += h;
        k1 = f(x_n, y[num_of_steps]);
        k2 = f(x_n + h*1 /5, y[num_of_steps] + h * (k1 * 1));
        k3 = f(x_n + h*3/10, y[num_of_steps] + h * (k1 * 1 / 4 + k2 * 3 / 4));
        k4 = f(x_n + h*4/5, y[num_of_steps] + h * (k1 * 11 / 9 + k2 * (-14. / 3) + k3 * 40 / 9));
        k5 = f(x_n + h*8/9, y[num_of_steps] + h * (k1 * 4843 / 1458 + k2 * (-3170. / 243) + k3 * 8056 / 729 + k4 * (-53. / 162)));
        k6 = f(x_n + h*1, y[num_of_steps] + h * (k1 * 9017 / 3168 + k2 * (-355. / 33) + k3 * 46732 / 5247 + k4 * (49. / 176) + k5 * (-5103. / 18656)));
        k7 = f(x_n + h*1, y[num_of_steps] + h * (k1 * 35 / 384 + k2 * 0 + k3 * 500 / 113 + k4 * 125 / 192 + k5 * (-2187. / 6784) + k6 * 11 / 84));
        y1 = y[num_of_steps] + h * (k1 * 5179. / 57600. + k2 * 0 + k3 * 7571. / 16695. + k4 * 393. / 640. + k5 * (-92097. / 339200) + k6 * 187 / 2100 + k7 * 1 / 40);
        y2 = y[num_of_steps] + h * (k1 * 35. / 384. + k2 * 0 + k3 * 500. / 1113. + k4 * 125. / 192. + k5 * (-2187. / 6784) + k6 * 11 / 84 + k7 * 0);
        s = pow(h*tolerance / (2*(x_f - x_0)*(y1 - y2).squaredNorm()), 0.25);

        if (s >= 2){
            y.push_back(y1);
            x.push_back(x_n);
            h *= 2;
            if(x_n + h > x_f)
                h = x_f - x_n;
            ++num_of_steps;
        }
        else if (s >= 1){
            y.push_back(y1);
            x.push_back(x_n);
            if(x_n + h > x_f)
                h = x_f - x_n;
            ++num_of_steps;
        }
        else if (s < 1){
            h /= 2;
        }
        std::cout<<num_of_steps<<" - "<< x_n <<std::endl;
    }
    return {x, y};
}

#endif //COMPUTATIONALMATHEMATICS_DORMANDPRINCE_HPP
