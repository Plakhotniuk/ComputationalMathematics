//
// Created by Арсений Плахотнюк on 05.02.2022.
//

#ifndef MY_PROJECT_THREEDIAGONALMATRIX_HPP
#define MY_PROJECT_THREEDIAGONALMATRIX_HPP
#include "compMath/Exceptions/SlaeBaseException.hpp"
#include <vector>
#include <array>
#include <sstream>
#include <string>

namespace Slae::Matrix{
    class ThreeDiagonalMatrix{

    private:
        std::vector<std::array<double, 3>> data_;

    public:

        /** @brief ThreeDiagonalMatrix class constructor
         * Creates three-diagonal matrix with size 'size'
         *
         * @param size -- matrix size
        */
        explicit ThreeDiagonalMatrix(unsigned int size);

        /** @brief ThreeDiagonalMatrix class static constructor
         * Creates three-diagonal zero matrix with size 'size'
         *
         * @param size -- matrix size
        */
        static ThreeDiagonalMatrix Zero(unsigned int size);

        /** @brief ThreeDiagonalMatrix class static constructor
         * Creates three-diagonal identity matrix with size 'size'
         *
         * @param size -- matrix size
        */
        static ThreeDiagonalMatrix Identity(unsigned int size);

        /** @brief ThreeDiagonalMatrix class static constructor
         * Creates three-diagonal matrix with size 'size' and filled with three numbers
         * @param size -- matrix size
         * @param a -- value on lower diagonal
         * @param b -- value on main diagonal
         * @param c -- value on upper diagonal
        */
        static ThreeDiagonalMatrix ThreeNumbers(unsigned int size, double a, double b, double c);

        /** @brief ThreeDiagonalMatrix class method
         *
         * @param i -- index of row
         * @param j -- index of column
         * @return matrix element with i,j indexes
        */
        double & operator()(int i, int j);

        /** @brief ThreeDiagonalMatrix class method
        *
        * @param i -- index of row
        * @param j -- index of column
        * @return matrix element with i,j indexes
        */
        [[nodiscard]] const double & operator()(int i, int j) const;

        /** @brief ThreeDiagonalMatrix class method
         * @return matrix row-size
        */
        [[nodiscard]] unsigned int rows() const noexcept;
        /**
         * @brief ThreeDiagonalMatrix class method
         * @param ind -- index of row
         * @param  a -- value on lower diagonal
         * @param b -- value on main diagonal
         * @param c -- value on upper diagonal
         */
        void fill_row(unsigned ind, double a, double b, double c);

        void check_diagonal_domimance() const;

    };

}   //namespace Slae::Matrix

#endif //MY_PROJECT_TREEDIAGONALMATRIX_HPP
