#pragma once

#include <vector>
#include <stdexcept>
#include <iostream>

// Matrix is stored as a contiguous array of rows
// a_11 a_12 ... a_1n a_21 a_22 ... a_2n ...       a_mn
// |______1 row_____| |______2 row_____| ...|_.. m_row_|

namespace Linal
{
    template <typename T>
    class Matrix
    {
        size_t m_rows = 0;
        size_t m_cols = 0;
        T* m_data = nullptr;

        struct ProxyRow
        {
            T* m_row = nullptr;
            size_t m_size = 0;
            T& operator[](size_t index);
            const T& operator[](size_t index) const;
        };

        static T* data_safe_copy(const T* src, size_t srcsize);

    public:

        // ctors
        Matrix() = default;
        explicit Matrix(size_t N);
        Matrix(size_t rows, size_t cols);
        Matrix(size_t rows, size_t cols, const std::vector<T>& data);
        Matrix(const Matrix& rhs);
        Matrix(Matrix&& rhs) noexcept;

        //assignment
        Matrix& operator=(const Matrix& rhs);
        Matrix& operator=(Matrix&& rhs) noexcept;

        //dtor
        ~Matrix() { delete[] m_data; }

        // selectors
        int nCols() const { return m_cols; }
        int nRows() const { return m_rows; }

        ProxyRow operator[] (size_t index);
        const ProxyRow operator[] (size_t index) const;

        //arithmetic operators
        Matrix& operator+= (const Matrix& rhs);
        Matrix& operator-= (const Matrix& rhs);
        Matrix& operator*= (const Matrix& rhs);
        Matrix& operator*= (const T& num);

        operator std::vector<T>();

        static T DeterminantNaive(const Matrix& matrix);

        Matrix& transpose();
    };

    template <typename T>
    T* Matrix<T>::data_safe_copy(const T* src, size_t srcsize)
    {
        T* dest = new T[srcsize];

        try
        {
            for(size_t i = 0; i < srcsize; ++i)
                dest[i] = src[i];
        }
        catch(...)
        {
            delete [] dest;
            throw;
        }

        return dest;
    }

    using iMatrix = Matrix<int>;
    using fMatrix = Matrix<float>;
    using dMatrix = Matrix<double>;

    template <typename T>
    T& Matrix<T>::ProxyRow::operator[](size_t index)
    {
        if(index >= m_size)
            throw std::invalid_argument("Out of bounds in column indexing");
        return m_row[index];
    }

    template <typename T>
    const T& Matrix<T>::ProxyRow::operator[](size_t index) const
    {
        if(index >= m_size)
            throw std::invalid_argument("Out of bounds in column indexing");
        return m_row[index];
    }

    template <typename T>
    Matrix<T>::Matrix(size_t N)
    :
    m_rows(N),
    m_cols(N),
    m_data(new T[m_rows*m_cols])
    {}

    template <typename T>
    Matrix<T>::Matrix(size_t rows, size_t cols)
    :
    m_rows(rows),
    m_cols(cols),
    m_data(new T[m_rows*m_cols]{})
    {}

    template <typename T>
    Matrix<T>::Matrix(size_t rows, size_t cols, const std::vector<T>& vec)
    :
    m_rows(rows),
    m_cols(cols),
    m_data(data_safe_copy(vec.data(), vec.size()))
    {
        if(vec.size() != m_rows * m_cols)
        {
            delete [] m_data;
            throw std::invalid_argument("Initializer vector dimension mismatch");
        }
    }

    template <typename T>
    Matrix<T>::Matrix(const Matrix& rhs)
    :
    m_rows(rhs.m_rows),
    m_cols(rhs.m_cols),
    m_data(data_safe_copy(rhs.m_data, rhs.m_rows * rhs.m_cols))
    {}

    template <typename T>
    Matrix<T>::Matrix(Matrix&& rhs) noexcept
    {
        std::cout << __PRETTY_FUNCTION__ << " was called" << std::endl;
        std::swap(m_cols, rhs.m_cols);
        std::swap(m_rows, rhs.m_rows);
        std::swap(m_data, rhs.m_data);
    }

    template <typename T>
    Matrix<T>& Matrix<T>::operator=(const Matrix& rhs)
    {
        if(this == &rhs)
            return *this;

        Matrix temp {rhs};

        *this = std::move(temp);

        return *this;
    }

    template <typename T>
    Matrix<T>& Matrix<T>::operator=(Matrix&& rhs) noexcept
    {
        if(this == &rhs)
            return *this;

        std::swap(m_cols, rhs.m_cols);
        std::swap(m_rows, rhs.m_rows);
        std::swap(m_data, rhs.m_data);

        return *this;
    }

    template <typename T>
    typename Matrix<T>::ProxyRow Matrix<T>::operator[] (size_t index)
    {
        if(index >= m_rows)
            throw std::invalid_argument("Out of bounds in row indexing");

        return ProxyRow {m_data + index * m_cols, m_cols};
    }


    template <typename T>
    typename Matrix<T>::ProxyRow const Matrix<T>::operator[] (size_t index) const
    {
        if(index >= m_rows)
            throw std::invalid_argument("Out of bounds in row indexing");

        return ProxyRow {m_data + index * m_cols, m_cols};
    }

    template <typename T>
    Matrix<T>& Matrix<T>::operator+=(const Matrix& rhs)
    {
        if(m_rows != rhs.m_rows || m_cols != rhs.m_cols)
            throw std::invalid_argument("Dimension of matrices in \"+=\" mismatch");

        for(size_t i = 0; i < m_rows * m_cols; ++i)
            m_data[i] += rhs.m_data[i];

        return *this;
    }

    template <typename T>
    Matrix<T> operator+(const Matrix<T>& lhs, const Matrix<T>& rhs)
    {
        Matrix<T> temp {lhs};
        temp += rhs;
        return temp;
    }

    template <typename T>
    Matrix<T>& Matrix<T>::operator-=(const Matrix& rhs)
    {
        if(m_rows != rhs.m_rows || m_cols != rhs.m_cols)
            throw std::invalid_argument("Dimension of matrices in \"-=\" mismatch");

        for(size_t i = 0; i < m_rows * m_cols; ++i)
            m_data[i] -= rhs.m_data[i];

        return *this;
    }

    template <typename T>
    Matrix<T> operator-(const Matrix<T>& lhs, const Matrix<T>& rhs)
    {
        Matrix<T> temp {lhs};
        temp -= rhs;
        return temp;
    }

    template <typename T>
    Matrix<T>& Matrix<T>::operator*=(const Matrix& rhs)
    {
        if(m_cols != rhs.m_rows)
            throw std::invalid_argument("Dimension of matrices in \"*=\" mismatch");

        Matrix<T> res {m_rows, rhs.m_cols};

        for(size_t i = 0; i < m_rows; ++i)
            for(size_t j = 0; j < rhs.m_cols; ++j)
                for(size_t k = 0; k < m_cols; ++k)
                    res[i][j] += (*this)[i][k] * rhs[k][j];

        *this = std::move(res);
        return *this;
    }

    template <typename T>
    Matrix<T> operator*(const Matrix<T>& lhs, const Matrix<T>& rhs)
    {
        Matrix<T> temp {lhs};
        temp *= rhs;
        return temp;
    }

    template <typename T>
    Matrix<T>& Matrix<T>::operator*=(const T& num)
    {
        for(size_t i = 0; i < m_rows * m_cols; ++i)
            m_data[i] *= num;

        return *this;
    }

    template <typename T>
    Linal::Matrix<T>::operator std::vector<T>()
    {
        std::vector<T> res = {};
        res.reserve(m_rows*m_cols);

        for(size_t i = 0; i < m_rows; ++i)
            for(size_t j = 0; j < m_cols; ++j)
                //res.push_back((*this)[i][j]);
                res.push_back(*(m_data + i * m_cols + j));

        return res;
    }

    template <typename T>
    T Matrix<T>::DeterminantNaive(const Matrix& matrix)
    {
        if(matrix.m_cols != matrix.m_rows)
            throw std::invalid_argument("DeterminantFailure: "
                                        "Matrix is not square");

        if(matrix.m_rows == 1)
            return matrix[0][0];

        T summ = {};

        size_t n_rows = matrix.m_rows;
        size_t n_cols = matrix.m_cols;

        for(size_t i = 0; i < n_cols; ++i)
        {
            // minor
            Matrix<T> temp {n_rows - 1, n_cols - 1};

            // minor filling
            for(size_t j = 0; j < n_rows - 1; j++)
                for(size_t k = 0; k < n_cols; k++)
                {
                    if(k == i)
                        continue;

                    if(k < i)
                        temp[j][k] = matrix[j+1][k];
                    else
                        temp[j][k - 1] = matrix[j+1][k];
                }

            summ += matrix[0][i] * DeterminantNaive(temp) * (i % 2 == 0 ? 1 : -1);
        }

        return summ;
    }


    template <typename T>
    Matrix<T>& Matrix<T>::transpose()
    {
        Matrix<T> temp {m_cols, m_rows};

        for(size_t i = 0; i < m_rows; ++i)
            for(size_t j = 0; j < m_cols; ++j)
                //temp[j][i] = (*this)[i][j];
                temp[j][i] = *(m_data + i * m_cols + j);

        *this = std::move(temp);

        return *this;
    }
} // namespace Linal

template <typename T>
std::ostream& operator<<(std::ostream& os, const Linal::Matrix<T>& matrix)
{
    size_t rows = matrix.nRows();
    size_t cols = matrix.nCols();

    for(size_t i = 0; i < rows; ++i)
    {
        for(size_t j = 0; j < cols; ++j)
            os << matrix[i][j] << " ";
        os << std::endl;
    }

    return os;
}
