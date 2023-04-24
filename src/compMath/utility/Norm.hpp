//
// Created by petrov on 19.02.2022.
//

#ifndef SLAE_NORM_HPP
#define SLAE_NORM_HPP

#include "vector"
#include "cmath"

enum NormType
{
    FirstNorm = 1,
    SecondNorm = 2,
    ThirdNorm = 3,
};
/**
 * Функция, реализующая различные нормы вектора
 * @tparam T шаблонный тип
 * @param vector вектор
 * @param normType тип нормы
 * @return norm соответствующая норма вектора
 */
template <typename T, int Type>
struct Norm{
    static double get_norm(const std::vector<T>& vector);
};

template<>
double Norm<double, NormType::FirstNorm>::get_norm(const std::vector<double> &vector) {
    double norm = 0.;
    for (const auto& elm: vector)
    {
        double abs = std::abs(elm);
        if (abs > norm) norm = abs;
    }
    return norm;
}

template<>
double Norm<double, NormType::SecondNorm>::get_norm(const std::vector<double> &vector) {
    double norm = 0.;
    for (const auto& elm: vector)
        norm += std::abs(elm);
    return norm;
}

template<>
double Norm<double, NormType::ThirdNorm>::get_norm(const std::vector<double> &vector) {
    double norm = 0.;
    for(int i = 0; i < vector.size(); ++i)
        norm += vector[i] * vector[i];
    return std::sqrt(norm);
}

#endif//SLAE_NORM_HPP
