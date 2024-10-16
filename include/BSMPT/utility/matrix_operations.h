#pragma once

#include <cassert>
#include <iostream>
#include <vector>

typedef std::vector<int> VecInt;

typedef std::vector<double> VecDoub;
typedef std::vector<VecDoub> MatDoub;
typedef std::vector<MatDoub> Mat3DDoub;

VecDoub operator*(const VecDoub &a, const double b);
VecDoub operator*(const double a, const VecDoub &b);
VecDoub operator+(const VecDoub &a, const VecDoub &b);
VecDoub operator*(const MatDoub &a, const VecDoub &b);

MatDoub operator+(const MatDoub &a, const MatDoub &b);
MatDoub operator-(const MatDoub &a, const MatDoub &b);
MatDoub operator*(const MatDoub &a, const MatDoub &b);

void set_zero(MatDoub &a);
void printvec(const VecDoub &a);
void printmat(const MatDoub &a);