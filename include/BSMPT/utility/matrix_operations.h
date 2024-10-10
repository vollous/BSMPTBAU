#pragma once

#include <iostream>
#include <vector>
#include <cassert>

typedef std::vector<double> VecDoub;
typedef std::vector<VecDoub> MatDoub;

VecDoub operator*(const VecDoub &a, const double b);
VecDoub operator*(const double a, const VecDoub &b);

VecDoub operator*(const MatDoub &a, const VecDoub &b);

MatDoub operator+(const MatDoub &a, const MatDoub &b);
MatDoub operator-(const MatDoub &a, const MatDoub &b);
MatDoub operator*(const MatDoub &a, const MatDoub &b);

void printvec(const VecDoub &a);
void printmat(const MatDoub &a);