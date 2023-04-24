//
// Created by Арсений Плахотнюк on 04.12.2022.
//

#ifndef HOMETASK2_NEWTON_HPP
#define HOMETASK2_NEWTON_HPP
#include <vector>
namespace NewtonInterpolation {

    double get_dif(const std::vector<double> &x, const std::vector<double> &f) {
        if (f.size() > 2) {
            std::vector<double> f_left(f.size() - 1);
            std::vector<double> x_left(x.size() - 1);
            std::vector<double> f_right(f.size() - 1);
            std::vector<double> x_right(x.size() - 1);
            for (int i = 0; i < f.size() - 1; ++i) {
                f_left[i] = f[i];
                x_left[i] = x[i];
                f_right[i] = f[i + 1];
                x_right[i] = x[i + 1];
            }
            return (get_dif(x_left, f_left) - get_dif(x_right, f_right)) / (x[x.size() - 1] - x[0]);
        } else if (f.size() == 2) {
            return (f[1] - f[0]) / (x[1] - x[0]);
        }
    }

    std::pair<double, std::vector<double>> get_val(std::vector<double> x_list, std::vector<double> f_list, double x) {
        double res = f_list[0];
        std::vector<double> diffs = {f_list[0]}; // Разделенные разности

        for (int i = 0; i < f_list.size(); ++i) {
            double buf = 1;
            std::vector<double> x_temp;
            std::vector<double> f_temp;
            for (int j = 0; j < i + 1; ++j) {
                x_temp.push_back(x_list[j]);
                f_temp.push_back(f_list[j]);
                if (j < i)
                    buf *= x - x_list[j];
            }
            diffs.push_back(get_dif(x_temp, f_temp));
            res += get_dif(x_temp, f_temp) * buf;
        }
        return {res, diffs};
    }
} // namespace NewtonInterpolation

#endif //HOMETASK2_NEWTON_HPP
