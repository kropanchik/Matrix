
#pragma once

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <tuple>
#include <ostream>
#include <cmath>
#include <string>

using namespace std;

class MyExeption {};

template <typename T>
class Matrix {
    template <class t1>
    friend class PCA;
protected:
    int m;
    int n;
    vector<T> a;
public:
    Matrix();
    Matrix(int m, int n);
    Matrix(int m, int n, vector<T> vec);

    void Create(int i, int j);
    void Create_with_data(int i, int j, vector<T> vec);
    [[nodiscard]] pair<int, int> M_Size() const;

    Matrix<T> AdamarMultiple(const Matrix<T> &obj);
    T Trace() const;
    [[nodiscard]] double FrobeniusNorm () const;
    [[nodiscard]] double EuclidNorm() const;
    T MaxNorm() const;
    T ScalarMul(const Matrix<T> &obj) const;
    Matrix<T> VectorMul(const Matrix<T> &obj);
    double Det();
    int Rank();
    double AngelofVec(const Matrix<T> &obj);
    void Transpose();
    void Transpose_vec();
    Matrix<T> Inverse();

    Matrix<T> operator*(const Matrix<T> &obj);
    Matrix<T> operator/( T r);
    Matrix<T> operator+(const Matrix<T> &obj);
    Matrix<T> operator-(const Matrix<T> &obj);
    Matrix<T>& operator=(const Matrix<T> &obj);



    template <typename T1>
    friend Matrix<T1> operator*(Matrix<T1> &obj, T1 r);
    template <typename T1>
    friend Matrix<T1> operator*(T1 r, Matrix<T1> &obj);
    template <typename T1>
    friend ostream &operator<<(ostream &out, const Matrix<T1> &matrix);
    template <typename T1>
    friend istream &operator>>(istream &in, Matrix<T1> &matrix);

    void write_to_bin_file(ostream&) const;
    void read_from_bin_file(istream&);
    void operator>> (ofstream &file) const;
    void operator<< (ifstream &file);

};

template <typename T>
class PCA
{
public:
    Matrix<T> mat;
    void Centering();
    void Normalize();
    double TRV();
    double ERV();
    double Leverage(int num);
    explicit PCA(Matrix<T> &m)
    {
        mat = m;
    }
    vector<Matrix<T>> NIPALS(int num_PC);

};



template <typename T>
class UnitMatrix: public Matrix<T>
{
public:
    using Matrix<T>::Matrix;
    using Matrix<T>::operator=;
    explicit UnitMatrix(int dim);

};

template <typename T>
class DiagonalMatrix :public Matrix<T>
{
public:
    using Matrix<T>::Matrix;
    using Matrix<T>::operator=;
    explicit DiagonalMatrix(int dim, vector<T> vec);

};

template <typename T>
class UpperTriangleMatrix :public Matrix<T>
{
public:
    using Matrix<T>::Matrix;
    using Matrix<T>::operator=;
};

template <typename T>
class LowerTriangleMatrix :public Matrix<T>
{
public:
    using Matrix<T>::Matrix;
    using Matrix<T>::operator=;
};

template <typename T>
class SymmetricMatrix :public Matrix<T>
{
public:
    using Matrix<T>::Matrix;
    using Matrix<T>::operator=;
};
