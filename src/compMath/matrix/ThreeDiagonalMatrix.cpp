//
// Created by Арсений Плахотнюк on 05.02.2022.
//

#include "ThreeDiagonalMatrix.hpp"
#include <iostream>
using std::cout;
using std::array;
using std::vector;

namespace Slae::Matrix{

    ThreeDiagonalMatrix::ThreeDiagonalMatrix(unsigned int size): data_(size) {};

    ThreeDiagonalMatrix ThreeDiagonalMatrix::Zero(unsigned int size)
    {
        ThreeDiagonalMatrix result(size);
        for (int i = 0; i < size; ++i) {
            result.data_[i] = {0., 0., 0.};
        }
        return result;
    }

    ThreeDiagonalMatrix ThreeDiagonalMatrix::Identity(unsigned int size)
    {
        ThreeDiagonalMatrix result(size);
        for (int i = 0; i < size; ++i) {
            result.data_[i] = {0., 1., 0.};
        }
        return result;
    }

    ThreeDiagonalMatrix ThreeDiagonalMatrix::ThreeNumbers(unsigned int size, double a, double b, double c) {
        ThreeDiagonalMatrix result(size);
        for (int i = 0; i < size; ++i) {
            result.data_[i] = {a, b, c};
        }
        return result;
    }

    unsigned int ThreeDiagonalMatrix::rows() const noexcept
    {
        return data_.size();
    }

    double & ThreeDiagonalMatrix::operator()(int i, int j)
    {
#ifndef NDEBUG
        if (i >= data_.size())
        {
            std::stringstream buff;
            buff << "Index i exceeds matrix size! Received index: " << i << ". Matrix size: "
                 << data_.size() << ". File: " << __FILE__ << ". Line: " << __LINE__;

            throw SlaeBaseExceptionCpp(buff.str());
        }

        if (j > 2)
        {

            std::stringstream buff;
            buff << "Index j must belong to {0, 1, 2}! Received index: " << i << ". Matrix size: "
                 << data_.size() << data_.size() << ". File: " << __FILE__ << ". Line: " << __LINE__;

            throw SlaeBaseExceptionCpp(buff.str());
        }
#endif
        return data_[i][j];
    }

    const double & ThreeDiagonalMatrix::operator()(int i, int j) const
    {
#ifndef NDEBUG
        if (i >= data_.size())
        {
            std::stringstream buff;
            buff << "Index i exceeds matrix size! Received index: " << i << ". Matrix size: "
                 << data_.size() << ". File: " << __FILE__ << ". Line: " << __LINE__;

            throw SlaeBaseExceptionCpp(buff.str());
        }

        if (j > 2)
        {
            std::stringstream buff;
            buff << "Index j must belong to {0, 1, 2}! Received index: " << j << ". Matrix size: "
                 << data_.size() << ". File: " << __FILE__ << ". Line: " << __LINE__;

            throw SlaeBaseExceptionCpp(buff.str());
        }
#endif
        return data_[i][j];

    }
    void ThreeDiagonalMatrix::fill_row(unsigned ind, double a, double b, double c){
        data_[ind][0] = a;
        data_[ind][1] = b;
        data_[ind][2] = c;
    }
    void ThreeDiagonalMatrix::check_diagonal_domimance() const {
        bool one_strict_inequality = false;
        bool non_strict_inequality = true;
        for(int i = 0; i < this->rows(); ++i){
            if(abs(data_[i][1]) <  abs(data_[i][0]) + abs(data_[i][2])){
                cout << "The sufficient condition of diagonal dominance is not fulfilled in row: " << i << std::endl;
                non_strict_inequality =!non_strict_inequality;
                break;
            }
            if(abs(data_[i][1]) >  abs(data_[i][0]) + abs(data_[i][2])){
                one_strict_inequality = !one_strict_inequality;
            }
        }
        if(one_strict_inequality && non_strict_inequality){
            cout << "The sufficient condition of diagonal dominance is fulfilled!"
                    "Strict inequality is satisfied in at least one row! " << std::endl;
        }
    }

}
