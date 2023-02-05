#pragma once

#include <random>
#include <type_traits>
#include "../matrix/matrix.hpp"

namespace Generators
{
    class MatrixGenerator
    {

        std::random_device m_rd;
        std::mt19937 m_generator;

    public:
        MatrixGenerator()
        :
        m_rd(),
        m_generator(m_rd())
        {}

        template <typename T>
        Linal::Matrix<T> generateUpperTriang(int n)
        {
            static_assert(std::is_same<T, double>::value ||
                          std::is_same<T, int>::value,
                           "Only double and int Matrix generation is supported");
            return Linal::Matrix<T>(1,1);
        }

        template <typename T>
        Linal::Matrix<T> generateWithExpectedDeterminant(int n, T determinant);
    };

    template <>
    Linal::Matrix<double> MatrixGenerator::generateUpperTriang(int n)
    {
        std::uniform_real_distribution distr {-0.5, 0.5};
        Linal::dMatrix res (n, n);

        for(int i = 0; i < n; ++i)
            for(int j = i; j < n; ++j)
                res[i][j] = distr(m_generator);

        return res;
    }

    template <>
    Linal::Matrix<int> MatrixGenerator::generateUpperTriang(int n)
    {
        std::uniform_int_distribution distr {-1000, 1000};
        Linal::iMatrix res (n, n);

        for(int i = 0; i < n; ++i)
            for(int j = i; j < n; ++j)
                res[i][j] = distr(m_generator);

        return res;
    }

    template<typename T>
    Linal::Matrix<T> MatrixGenerator::generateWithExpectedDeterminant(int n, T determinant)
    {
        Linal::Matrix<T> m1 = generateUpperTriang<T>(n);
        Linal::Matrix<T> m2 = generateUpperTriang<T>(n);


        for(int i = 0; i < n; ++i)
            m2[i][i] = m1[i][i] = 1;

        m2[0][0] = determinant;

        m2.transpose();
        return m1 * m2;
    }
}
