//
// Created by Арсений Плахотнюк on 07.02.2022.
//

#include "ThreeDiadonalSolver.hpp"
#include <iostream>

namespace Slae::Solvers
{
    std::vector<double> solveThreeDiagonal(const Slae::Matrix::ThreeDiagonalMatrix &matrix, const std::vector<double> &col)
    {
        if (matrix.rows() != col.size()) {
            std::stringstream buff;
            buff << "Matrix and column dimensions are not equal! Matrix size: " << matrix.rows() << ". Column size: "
                 << col.size() << ". File: " << __FILE__ << ". Line: " << __LINE__;
            throw SlaeBaseExceptionCpp(buff.str());
        }

        matrix.check_diagonal_domimance();

        int n = col.size();

        std::vector<double> result(n);

        std::vector<std::array<double, 2>> pr(n - 1);

        pr[0][0] = matrix(0,2) / matrix(0,1);

        pr[0][1] = col[0] /  matrix(0, 1);

        for (int i = 1; i < n - 1; ++i){

            pr[i][0] = matrix(i,2) / (matrix(i,1) - matrix(i,0) * pr[i - 1][0]);

            pr[i][1] = (col[i] - matrix(i, 0) * pr[i-1][1]) / (matrix(i,1) - matrix(i,0) * pr[i - 1][0]);

        }

        result[n-1] = (col[n-1] - matrix(n-1, 0) * pr[n-2][1]) / (matrix(n-1,1) - matrix(n-1,0) * pr[n - 2][0]);;

        for (int i = n-2; i >= 0; --i){
            result[i] = pr[i][1] - pr[i][0] * result[i+1];
        }

        return result;

    }
}