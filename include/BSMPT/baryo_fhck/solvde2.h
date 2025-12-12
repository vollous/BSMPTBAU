#pragma once

#include <BSMPT/utility/data_structures.h>
#include <fstream>

class LUdcmp {
    private:
        size_t n;
        MatDoub lu, &aref;
        VecInt indx;
        double d;

    public:
        LUdcmp(MatDoub &a);
        void solve(VecDoub &b, VecDoub &x);
        void solve(MatDoub &b, MatDoub &x);
        void inverse(MatDoub &ainv);
        double det();
        ~LUdcmp() {};
};

struct VecMat{
    VecDoub vec;
    MatDoub mat;
    VecMat(VecDoub &v, MatDoub &m) : vec(v), mat(m) {
        if((v.size() != m.rows()) || (v.size() != m.cols())) {
            std::cout << "Error in VecMat: Vector and Matrix are not of same dimension\n";
            exit(1); 
        }
    };
    ~VecMat() {};
};

// based on the Thomas algorithm for matrices
class Lin1ODEBVPSolver {
    private:
        const size_t N, dim;
        MatDoub Mtemp;
        VecDoub Vtemp;
        std::vector<MatDoub> C;
    public:
        Lin1ODEBVPSolver(std::vector<VecMat> &ode, const double &x1, const double &x2, VecDoub &y1, VecDoub &y2);
        void SetCMatrices();
        ~Lin1ODEBVPSolver() {};
};