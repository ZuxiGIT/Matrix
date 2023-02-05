#pragma once

#include "../matrix/matrix.hpp"

namespace Generators
{
    class MatrixGenerator
    {


    public:

        template <typename T>
        Linal::Matrix<T> GenerateUpperTriang(int n);
    };

    template <>
    Linal::Matrix<double> GenerateUpperTriang(int n)
    {

    }
}
