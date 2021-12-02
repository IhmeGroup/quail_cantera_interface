#include "helpers.h"

extern "C" {

void get_massfractions(const int nsp, const double rho, const double *rhoY, double *Y)
{
    double Y_last = 1.0;
    for (int is = 0; is < nsp - 1; is++){
        Y[is] = rhoY[is] / rho;
        Y_last += -Y[is];
    }
    Y[nsp - 1] = Y_last;
}

} // end extern "C"