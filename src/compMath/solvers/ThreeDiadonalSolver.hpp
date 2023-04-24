//
// Created by Арсений Плахотнюк on 07.02.2022.
//

#ifndef MY_PROJECT_THREEDIADONALSOLVER_HPP
#define MY_PROJECT_THREEDIADONALSOLVER_HPP

#include "compMath/matrix/ThreeDiagonalMatrix.hpp"
#include "compMath/Exceptions/SlaeBaseException.hpp"
#include <sstream>
#include <vector>

namespace Slae::Solvers
{
    /** @brief Метод решает систему уравнений при помощи метода прогонки
* Решает систему линейных алгебраических уравнений при помощи метода прогонки. О методе прогонки можно узнать из
* https://mipt.ru/education/chair/computational_mathematics/study/materials/compMath/other/Aristova_Zavyalova_Lobanov_2014.pdf
*
* @param matrix трехдиагональная матрица
* @param col столбец правой части
* @return решение СЛАУ
*
* @throw SlaeBaseExceptionCpp выбрасывается, если количество строк матрицы и высота столбца не совпадают
*/
    [[nodiscard]] std::vector<double> solveThreeDiagonal(const Matrix::ThreeDiagonalMatrix &matrix, const std::vector<double> &col);

}  // namespace Slae::Solvers

#endif //MY_PROJECT_THREEDIADONALSOLVER_HPP
