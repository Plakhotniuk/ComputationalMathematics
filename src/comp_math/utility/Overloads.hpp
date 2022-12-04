//
// Created by petrov on 12.02.2022.
//

#ifndef SLAE_OVERLOADS_HPP
#define SLAE_OVERLOADS_HPP
#include <iosfwd>
#include <vector>
#include "differential_equations/Exceptions/SlaeBaseException.hpp"
/***
 * Умножение числа на вектор
 * @tparam T Тип числа
 * @param k Число
 * @param b Вектор
 * @return Результат умножения
 */
template<typename T>
std::vector<T> operator*(const T &k, const std::vector<T> &b){
    std::vector<T> res(b.size());
    for (size_t i = 0; i < b.size(); ++i) {
        res[i] = b[i] * k;
    }
    return res;
};

template<typename T>
std::vector<T> operator/(const std::vector<T> &b, const T &k){
    std::vector<T> res(b.size());
    for (size_t i = 0; i < b.size(); ++i) {
        res[i] = b[i] / k;
    }
    return res;
};

/***
 * Умножение вектора на число
 * @tparam T Тип числа
 * @param b Вектор
 * @param k Число
 * @return Результат умножения
 */
template<typename T>
std::vector<T> operator*(const std::vector<T> &b, const T &k){
    return k * b;
};

/***
 * Скалярное перемножение векторов
 * @tparam T Тип элементов векторов
 * @param a Первый вектор
 * @param b Второй вектор
 * @return Результат скалярного перемножения
 */
template<typename T>
T operator*(const std::vector<T> &a, const std::vector<T> &b){
#ifndef NDEBUG
    if (a.size() != b.size()) {
        std::stringstream buff;
        buff << "Vectors must have same size! A.size " << a.size() << "! B.size " << b.size() << "!" << std::endl;
        throw Slae::SlaeBaseExceptionCpp(buff.str());
    }
#endif //NDEBUG
    T res = 0;
    for (int i = 0; i < a.size(); ++i) {
        res += a[i] * b[i];
    }
    return res;
};

/***
 * Вычитание двух векторов
 * @tparam T Тип элементов векторов
 * @param a Первый вектор
 * @param b Второй вектор
 * @return Результат вычитания векторов
 */
template<typename T>
std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b){
#ifndef NDEBUG
    if (a.size() != b.size()) {
        std::stringstream buff;
        buff << "Vectors must have same size! A.size " << a.size() << "! B.size " << b.size() << "!" << std::endl;
        throw Slae::SlaeBaseExceptionCpp(buff.str());
    }
#endif //NDEBUG
    std::vector<T> res(a.size());
    for(int i = 0; i < a.size(); ++i){
        res[i] = a[i] - b[i];
    }
    return res;
};

/***
 * Сложение двух векторов
 * @tparam T Тип элементов векторов
 * @param a Первый вектор
 * @param b Второй вектор
 * @return Результат сложения векторов
 */
template<typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b){
    return a - static_cast<T>(-1) * b;
};

/**
 * оператор вывода вектора в поток
 * @tparam T используемые тип данных
 * @param os поток
 * @param b вектор
 * @return поток с вектором
 */
template<typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &b){
    for (auto i: b) {
        os << i << " ";
    }
    os << "\n";
    return os;
};

#endif//SLAE_OVERLOADS_HPP
