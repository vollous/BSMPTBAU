#include <BSMPT/baryo_fhck/solvde2.h>

LUdcmp::LUdcmp(MatDoub &a) : n(a.rows()), lu(a), aref(a), indx(n) {
    const double TINY = 1.0e-40;
    double big, temp;
    VecDoub vv(n);  // pivot coeffs for each row
    d = 1.0;        // sign
    // find the pivot for each row
    for (size_t i = 0; i < n; i++) {
        big = 0.0;
        for (size_t j = 0; j < n; j++){
            temp = std::abs(lu[i][j]);
            if (temp > big) big = temp;
        }
        if (big == 0.0) throw("Singular matrix in LUdcmp");
        vv[i] = 1.0 / big;
    }
    for (size_t k = 0; k < n; k++) {
        big = 0.0;
        size_t imax = k;
        for (size_t i = k; i < n; i++) {
            temp = vv[i] * std::abs(lu[i][k]);
            if (temp > big) {
                big = temp;
                imax = i;
            }
        }
        if (k != imax) {
            for (size_t j = 0; j < n; j++) {
                SWAP(lu[imax][j], lu[k][j]);
            }
            d = -d;
            vv[imax] = vv[k];
        }
        indx[k] = imax;
        if (lu[k][k] == 0.0) lu[k][k] = TINY;
        for (size_t i = k + 1; i < n; i++) {
            temp = lu[i][k] /= lu[k][k];
            for (size_t j = k + 1; j < n; j++)
                lu[i][j] -= temp * lu[k][j];
        }
    }
}

void LUdcmp::solve(VecDoub &b, VecDoub &x) {
    size_t ii = 0, ip;
    double sum = 0;
    if (b.size() != n || x.size() != n) throw("LUdcmp::solve bad sizes");
    for (size_t i = 0; i < n; i++) x[i] = b[i];
    for (size_t i = 0; i < n; i++) {
        ip = indx[i];
        sum = x[ip];
        x[ip] = x[i];
        if (ii != 0)
            for (size_t j = ii - 1; j < i; j++) sum -= lu[i][j] * x[j];
        else if (sum != 0.0)
            ii = i + 1;
        x[i] = sum;
    }
    for (size_t i = n; i-- > 0;) {
        sum = x[i];
        for (size_t j = i + 1; j < n; j++) sum -= lu[i][j] * x[j];
        x[i] = sum / lu[i][i];
    }
}

void LUdcmp::solve(MatDoub &b, MatDoub &x) {
    size_t m = b.cols();
    if (b.rows() != n || x.rows() != n || b.cols() != x.cols()) throw("LUdcmp::solve bad sizes");
    VecDoub xx(n);
    for (size_t j = 0; j < m; j++) {
        for (size_t i = 0; i < n; i++) xx[i] = b[i][j];
        solve(xx, xx);
        for (size_t i = 0; i < n; i++) x[i][j] = xx[i];
    }
}

void LUdcmp::inverse(MatDoub &ainv) {
    ainv.resize(n, n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) ainv[i][j] = 0.;
        ainv[i][i] = 1.;
    }
    solve(ainv, ainv);
}

double LUdcmp::det() {
    double dd = d;
    for (size_t i = 0; i < n; i++) dd *= lu[i][i];
    return dd;
}


Lin1ODEBVPSolver::Lin1ODEBVPSolver(std::vector<VecMat> &ode, const double &x1, const double &x2, VecDoub &y1, VecDoub &y2) 
: N(ode.size()), dim(ode[0].vec.size()) 
{
    const double h = (x2 - x1) / (double) (N + 1);

    // set C matrices to diagonal unit matrices
    SetCMatrices();

    // rescale ode matrices and vectors by 2h
    for (size_t i = 0; i < N; i++){
        ode[i].vec = ode[i].vec * (2. * h);
        ode[i].mat = ode[i].mat * (2. * h);
    }

    // set first boundary condition
    ode[0].vec = ode[0].vec + y1;

    // set second boundary condition
    ode[N - 1].vec = ode[N - 1].vec - y2;

    // first step:
    Mtemp = ode[0].mat;
    LUdcmp lu(Mtemp);
    lu.inverse(C[0]);
    ode[0].vec = C[0] * ode[0].vec;
    
    // i'th step
    for(size_t i = 1; i < N; i++){
        Mtemp = ode[i].mat;
        Mtemp = Mtemp + C[i - 1];
        LUdcmp dcmp(Mtemp);
        dcmp.inverse(Mtemp);
        if (i != N - 1)
            C[i] = Mtemp;
        ode[i].vec = ode[i].vec + ode[i - 1].vec;
        ode[i].vec = Mtemp * ode[i].vec;
    }

    // Back substitution step
    for (int i = (int)(N - 2); i >= 0; i--) {
        Vtemp = C[i] * ode[i + 1].vec;
        ode[i].vec = ode[i].vec - Vtemp;
    }

    std::ofstream file("test1.dat");
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < dim; j++)
            file << ode[i].vec[j] << "\t";
        file << "\n";
    }
    file.close();
}

void Lin1ODEBVPSolver::SetCMatrices() {
    for (size_t i = 0; i < N; i++) {
        C.push_back(MatDoub(dim, dim, 0.));
        for (size_t j = 0; j < dim; j++) 
            C[i][j][j] = 1.;
    }
}