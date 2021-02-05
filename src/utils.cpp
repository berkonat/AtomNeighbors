#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <functional>
#include <numeric>
#include <algorithm>

#include "utils.hpp"

std::vector<double> read_radius(std::string fname)
{
    std::ifstream infile(fname);
    double content;
    std::string str;
    std::vector<double> X;
    int c=0;
    while (std::getline(infile, str)){
        std::cout << str << std::endl;
        if (c==1) break;
        c++;
    }
    while (infile >> content)
    {
        X.push_back(content);
        //std::cout << content << std::endl;
    }
    return X;
}

std::vector<double> read_xyz(std::string fname,
                             std::vector<double>* box_lengths)
{
    std::ifstream infile(fname);
    double content;
    std::string str;
    std::vector<double> X;
    int i;
    std::getline(infile, str);
    std::getline(infile, str);
    std::cout << str << std::endl;
    for(i=0;i<3;i++) {
        infile >> content;
        box_lengths->push_back(content);
    }
    while (infile >> content)
    {
        std::vector<double> row;
        row.push_back(content);
        for(i=1;i<3;i++){
            if(!(infile >> content)){
                break;
            } else {
                row.push_back(content);
            }
        }
        X.insert(X.end(), row.begin(), row.end());
        //std::cout << row[0] << ' ' << row[1] << ' ' << row[2] << std::endl;
    }
    return X;
}

template<typename T>
std::vector<T> vecslice(std::vector<T> const &A, int start, int step)
{
    std::vector<T> vec(A.cbegin() + start, A.cbegin() + start + step + 1);
    return vec;
}


template<typename T>
std::vector<T> vecslice(std::vector<T> const &A, int start, int step, int stride)
{
    std::vector<T> vec(step);
    for(int i=0;i<step;i++)
        vec[i] = A[i*stride + start];
    return vec;
}


template<typename T1, typename T2>
std::vector<T1> gemv(std::vector<T1> const &A, std::vector<T2> const &x)
{
    int i;
    int d = x.size();
    int nA = A.size()/d;
    std::vector<T1> u(nA, 0.0);
    for(i=0;i<nA;i++){
        u[i] = std::inner_product(&A[i*d], &A[i*d+d], x.begin(), 0);
    }
    return u;
}

template<typename T1, typename T2>
std::vector<T1> gemv(std::vector<T1> const &A, std::vector<T2> const &B,
                     int ncolA, int* axis)
{
    int i,j;
    int ax = 0;
    if(axis != NULL) ax = *axis;
    int nrowA = A.size()/ncolA;
    int nB = B.size()/ncolA;
    std::vector<T1> M(nB * nrowA, 0.0);
    if(ax == 0){
        for(j=0;j<nB;j++)
            for(i=0;i<nrowA;i++)
                M[j*ncolA+i] += std::inner_product(&A[i*ncolA], &A[i*ncolA+ncolA], &B[j*ncolA], 0);
    } else {
        std::vector<T2> u(ncolA);
        for(j=0;j<nB;j++){
            for(i=0;i<ncolA;i++) u[i] = B[i*nB + j];
            for(i=0;i<nrowA;i++)
                M[j+i*nB] += std::inner_product(&A[i*ncolA], &A[i*ncolA+ncolA], &B[j*ncolA], 0);
        }
    }
    return M;
}

template<typename T>
double norm(std::vector<T> const &u)
{
    return sqrt(std::inner_product(u.begin(), u.end(), u.begin(), 0));
}

template<typename T>
std::vector<double> norm(std::vector<T> const &A, int ncolA, int* axis)
{
    if (axis != NULL){
        int nrowA = A.size()/ncolA;
        std::vector<double> rtn(nrowA, 0.);
        int i;
        if ((*axis) == 1){
            std::vector<double> u;
            for(i=0; i<nrowA; i++){
                u = vecslice(A,i,ncolA,nrowA);
                rtn[i] = sqrt(std::inner_product(u.begin(), u.end(), u.begin(), 0));
            }
        } else {
            for(i=0; i<nrowA; i++){
                rtn[i] = sqrt(std::inner_product(&A[i*ncolA], &A[i*ncolA+ncolA], &A[i*ncolA], 0));
            }
        }
        return rtn;
    } else {
        double sm = sqrt(std::inner_product(A.begin(), A.end(), A.begin(), 0));
        std::vector<double> rtn = {sqrt(sm)};
        return rtn;
    }
}

template<typename T1, typename T2>
std::vector<T1> vsub(std::vector<T1> const &u, T2 const &v)
{
    std::vector<T1> rtn(u.size());
    for(int i=0; i<u.size(); i++) rtn[i] = u[i] - *v;
    return rtn;
}

template<typename T1, typename T2>
std::vector<T1> vsub(std::vector<T1> const &u, std::vector<T2> const &v)
{
    std::vector<T1> rtn(u.size());
    for(int i=0; i<u.size(); i++) rtn[i] = u[i] - v[i];
    return rtn;
}

template<typename T1, typename T2>
std::vector<T1> vsub(std::vector<T1> const &A, std::vector<T2> const &v, int ncolA)
{
    std::vector<T1> rtn(A.size());
    for(int i=0; i<A.size(); i++) rtn[i] = A[i%ncolA] - v[i];
    return rtn;
}

template<typename T1, typename T2>
std::vector<T1> vadd(std::vector<T1> const &u, T2 const &v)
{
    std::vector<T1> rtn(u.size());
    for(int i=0; i<u.size(); i++) rtn[i] = u[i] + *v;
    return rtn;
}

template<typename T1, typename T2>
std::vector<T1> vadd(std::vector<T1> const &u, std::vector<T2> const &v)
{
    std::vector<T1> rtn(u.size());
    for(int i=0; i<u.size(); i++) rtn[i] = u[i] + v[i];
    return rtn;
}

template<typename T1, typename T2>
std::vector<T1> vadd(std::vector<T1> const &A, std::vector<T2> const &v, int ncolA)
{
    std::vector<T1> rtn(A.size());
    for(int i=0; i<A.size(); i++) rtn[i] = A[i%ncolA] + v[i];
    return rtn;
}

template<typename T1, typename T2>
std::vector<T1> vmul(std::vector<T1> const &u, std::vector<T2> const &v)
{
    std::vector<T1> rtn(u.size());
    for(int i=0; i<u.size(); i++) rtn[i] = u[i] * v[i];
    return rtn;
}

template<typename T1, typename T2>
std::vector<T1> vmul(std::vector<T1> const &A, std::vector<T2> const &v, int ncolA)
{
    std::vector<T1> rtn(A.size());
    for(int i=0; i<A.size(); i++) rtn[i] = A[i%ncolA] * v[i];
    return rtn;
}

template<typename T1, typename T2>
std::vector<T1> vmul(std::vector<T1> const &u, T2 const &v)
{
    std::vector<T1> rtn(u.size());
    for(int i=0; i<u.size(); i++) rtn[i] = u[i] + *v;
    return rtn;
}

template<typename T1, typename T2>
std::vector<T1> vdiv(std::vector<T1> const &u, std::vector<T2> const &v)
{
    std::vector<T1> rtn(u.size());
    for(int i=0; i<u.size(); i++) rtn[i] = u[i] / v[i];
    return rtn;
}

template<typename T1, typename T2>
std::vector<T1> vdiv(std::vector<T1> const &A, std::vector<T2> const &v, int ncolA)
{
    std::vector<T1> rtn(A.size());
    for(int i=0; i<A.size(); i++) rtn[i] = A[i%ncolA] / v[i];
    return rtn;
}

template<typename T1, typename T2>
std::vector<T1> vdiv(std::vector<T1> const &u, T2 const v)
{
    std::vector<T1> rtn(u.size());
    for(int i=0; i<u.size(); i++) rtn[i] = u[i] / *v;
    return rtn;
}

//std::vector<double> vadd(std::vector<double> &u,
//                         double &v)
//{
//    std::vector<double> rtn;
//    std::transform(u.begin(), u.end(),
//                   std::back_inserter(rtn),
//                   [v](auto& c){return c + v;});
//                   //std::plus<double>());
//                   //std::minus<T1>());
//                   //std::multiplies<double>());
//                   //std::divides<double>());
//    //for(auto& item : u){
//    //    rtn.push_back(item + *v);
//    return rtn;
//}

// MAcros for instantiate functions for each type

#define local _VTYPES = "double float" \
        local _VTYPES2 = "int" \
foreach VT of _VTYPES { \
    template std::vector<VT> vecslice(std::vector<VT> const &A, int start, int step); \
    template std::vector<VT> vecslice(std::vector<VT> const &A, int start, int step, int stride); \
    template double norm(std::vector<VT> const &u); \
    template std::vector<double> norm(std::vector<VT> const &A, int ncolA, int* axis); \
    foreach VT2 of _VTYPES { \
        template std::vector<VT> gemv(std::vector<VT> const &A, std::vector<VT2> const &x); \
        template std::vector<VT> gemv(std::vector<VT> const &A, std::vector<VT2> const &B, int ncolA, int* axis); \
        template std::vector<VT> vadd(std::vector<VT> const &u, VT2 const &v); \
        template std::vector<VT> vsub(std::vector<VT> const &u, VT2 const &v); \
        template std::vector<VT> vmul(std::vector<VT> const &u, VT2 const &v); \
        template std::vector<VT> vdiv(std::vector<VT> const &u, VT2 const &v); \
        template std::vector<VT> vadd(std::vector<VT> const &u, std::vector<VT2> const &v); \
        template std::vector<VT> vsub(std::vector<VT> const &u, std::vector<VT2> const &v); \
        template std::vector<VT> vmul(std::vector<VT> const &u, std::vector<VT2> const &v); \
        template std::vector<VT> vdiv(std::vector<VT> const &u, std::vector<VT2> const &v); \
        template std::vector<VT> vadd(std::vector<VT> const &A, std::vector<VT2> const &v, int ncolA); \
        template std::vector<VT> vsub(std::vector<VT> const &A, std::vector<VT2> const &v, int ncolA); \
        template std::vector<VT> vmul(std::vector<VT> const &A, std::vector<VT2> const &v, int ncolA); \
        template std::vector<VT> vdiv(std::vector<VT> const &A, std::vector<VT2> const &v, int ncolA); \
    } \
    foreach VT2 of _VTYPES2 { \
        template std::vector<VT> gemv(std::vector<VT> const &A, std::vector<VT2> const &x); \
        template std::vector<VT> gemv(std::vector<VT> const &A, std::vector<VT2> const &B, int ncolA, int* axis); \
        template std::vector<VT> vadd(std::vector<VT> const &u, VT2 const &v); \
        template std::vector<VT> vsub(std::vector<VT> const &u, VT2 const &v); \
        template std::vector<VT> vmul(std::vector<VT> const &u, VT2 const &v); \
        template std::vector<VT> vdiv(std::vector<VT> const &u, VT2 const &v); \
        template std::vector<VT> vadd(std::vector<VT> const &u, std::vector<VT2> const &v); \
        template std::vector<VT> vsub(std::vector<VT> const &u, std::vector<VT2> const &v); \
        template std::vector<VT> vmul(std::vector<VT> const &u, std::vector<VT2> const &v); \
        template std::vector<VT> vdiv(std::vector<VT> const &u, std::vector<VT2> const &v); \
        template std::vector<VT> vadd(std::vector<VT> const &A, std::vector<VT2> const &v, int ncolA); \
        template std::vector<VT> vsub(std::vector<VT> const &A, std::vector<VT2> const &v, int ncolA); \
        template std::vector<VT> vmul(std::vector<VT> const &A, std::vector<VT2> const &v, int ncolA); \
        template std::vector<VT> vdiv(std::vector<VT> const &A, std::vector<VT2> const &v, int ncolA); \
    } \
} \
foreach VT of _VTYPES { \
    foreach VT2 of _VTYPES2 { \
        template std::vector<VT> gemv(std::vector<VT2> const &A, std::vector<VT> const &x); \
        template std::vector<VT> gemv(std::vector<VT2> const &A, std::vector<VT> const &B, int ncolA, int* axis); \
        template std::vector<VT> vadd(std::vector<VT2> const &u, VT const &v); \
        template std::vector<VT> vsub(std::vector<VT2> const &u, VT const &v); \
        template std::vector<VT> vmul(std::vector<VT2> const &u, VT const &v); \
        template std::vector<VT> vdiv(std::vector<VT2> const &u, VT const &v); \
        template std::vector<VT> vadd(std::vector<VT2> const &u, std::vector<VT> const &v); \
        template std::vector<VT> vsub(std::vector<VT2> const &u, std::vector<VT> const &v); \
        template std::vector<VT> vmul(std::vector<VT2> const &u, std::vector<VT> const &v); \
        template std::vector<VT> vdiv(std::vector<VT2> const &u, std::vector<VT> const &v); \
        template std::vector<VT> vadd(std::vector<VT2> const &A, std::vector<VT> const &v, int ncolA); \
        template std::vector<VT> vsub(std::vector<VT2> const &A, std::vector<VT> const &v, int ncolA); \
        template std::vector<VT> vmul(std::vector<VT2> const &A, std::vector<VT> const &v, int ncolA); \
        template std::vector<VT> vdiv(std::vector<VT2> const &A, std::vector<VT> const &v, int ncolA); \
    } \
}


