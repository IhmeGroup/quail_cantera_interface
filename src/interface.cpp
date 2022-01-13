// #include "cantera/thermo.h"
#include <iostream>
#include <omp.h>
#include "thermo.h"
#include "helpers.h"
#include <limits>
#include "cantera/kinetics.h"

extern "C" {

using namespace Cantera;


void get_pressure_interface(const double* Uq, double* P, int ne, int nq, 
    int ns, int nsp, int dim, char* filename)
{

// Create a new phase
// std::unique_ptr<IdealGasMix> gas(new IdealGasMix(filename));
    std::unique_ptr<ThermoPhase> gas(newPhase(filename));
    for (int ie = 0; ie < ne; ie++){
        for (int iq = 0; iq < nq; iq++){

            auto start = ie * (nq*ns) + iq * ns;
            const double *U = Uq + start;
            const double rho = U[0];
            const double rhoE = U[dim + 1];
            double rhoKE = 0;
            for (int id = 0; id < dim; id++){
                rhoKE = rhoKE + U[id + 1] * U[id + 1] / rho;
            }
            auto e = (rhoE - 0.5 * rhoKE) / rho;

            double Y[nsp];
            get_massfractions(nsp, rho, U + dim + 2, Y);
            P[ie * (nq) + iq] = get_P(gas, rho, e, Y);
        }
    }
} // end get_pressure_interface

void get_wavespeed_interface(const double* Uq, double* gamma, double* P, int ne, int nq, 
    int ns, int nsp, int dim, char* filename)
{

// Create a new phase
std::unique_ptr<ThermoPhase> gas(newPhase(filename));

// #pragma omp parallel for schedule(static, 1)
    for (int ie = 0; ie < ne; ie++){
        for (int iq = 0; iq < nq; iq++){

            auto start = ie * (nq*ns) + iq * ns;
            const double *U = Uq + start;
            const double rho = U[0];
            const double rhoE = U[dim + 1];
            double rhoKE = 0;
            for (int id = 0; id < dim; id++){
                rhoKE = rhoKE + U[id + 1] * U[id + 1] / rho;
            }
            auto e = (rhoE - 0.5 * rhoKE) / rho;

            double Y[nsp];
            get_massfractions(nsp, rho, U + dim + 2, Y);
            P[ie * (nq) + iq] = get_P(gas, rho, e, Y);
            gamma[ie * (nq) + iq] = get_gamma(gas);
        }
    }
} // end get_wavespeed_interface

void get_temperature_interface(const double* Uq, double* T, int ne, int nq, 
    int ns, int nsp, int dim, char* filename)
{

// Create a new phase
std::unique_ptr<ThermoPhase> gas(newPhase(filename));

// #pragma omp parallel for schedule(static, 1)
    for (int ie = 0; ie < ne; ie++){
        for (int iq = 0; iq < nq; iq++){
            
            auto start = ie * (nq*ns) + iq * ns;
            const double *U = Uq + start;
            const double rho = U[0];
            const double rhoE = U[dim + 1];
            double rhoKE = 0;
            for (int id = 0; id < dim; id++){
                rhoKE = rhoKE + U[id + 1] * U[id + 1] / rho;
            }
            auto e = (rhoE - 0.5 * rhoKE) / rho;

            double Y[nsp];
            get_massfractions(nsp, rho, U + dim + 2, Y);
            T[ie * (nq) + iq] = get_T(gas, rho, e, Y);

        } // end nq loop
    } // end ne loop
} // end get_temperature_interface

void get_gamma_interface(const double* Uq, double* gamma, int ne, int nq, 
    int ns, int nsp, int dim, char* filename)
{

// Create a new phase
std::unique_ptr<ThermoPhase> gas(newPhase(filename));

// #pragma omp parallel for schedule(static, 1)
    for (int ie = 0; ie < ne; ie++){
        for (int iq = 0; iq < nq; iq++){
            
            auto start = ie * (nq*ns) + iq * ns;
            const double *U = Uq + start;
            const double rho = U[0];
            const double rhoE = U[dim + 1];
            double rhoKE = 0;
            for (int id = 0; id < dim; id++){
                rhoKE = rhoKE + U[id + 1] * U[id + 1] / rho;
            }
            auto e = (rhoE - 0.5 * rhoKE) / rho;

            double Y[nsp];
            get_massfractions(nsp, rho, U + dim + 2, Y);
            double P = get_P(gas, rho, e, Y);
            gamma[ie * (nq) + iq] = get_gamma(gas);

        } // end nq loop
    } // end ne loop
} // end get_temperature_interface

void get_net_production_rates_interface(const double* Uq, double* S, int ne, int nq, 
    int ns, int nsp, int dim, char* filename)
{

// Create a new phase
std::unique_ptr<ThermoPhase> gas(newPhase(filename));

std::vector<ThermoPhase *> phases_{gas.get()};
std::unique_ptr<Kinetics> kin(newKineticsMgr(gas->xml(), phases_));
double T;
double Wk;
// #pragma omp parallel for schedule(static, 1)
    for (int ie = 0; ie < ne; ie++){
        for (int iq = 0; iq < nq; iq++){
            
            auto start = ie * (nq*ns) + iq * ns;
            const double *U = Uq + start;
            const double rho = U[0];
            const double rhoE = U[dim + 1];
            double rhoKE = 0;
            for (int id = 0; id < dim; id++){
                rhoKE = rhoKE + U[id + 1] * U[id + 1] / rho;
            }
            auto e = (rhoE - 0.5 * rhoKE) / rho;

            double Y[nsp];
            double wdot[nsp];
            get_massfractions(nsp, rho, U + dim + 2, Y);
            T = get_T(gas, rho, e, Y);

            kin->getNetProductionRates(wdot);

	       // multiply wdot by molecular weights
	       for (int is = 0; is < nsp - 1; is++){
                Wk = gas->molecularWeight(is);
                S[ie*nq*ns + iq * ns + is + dim + 2] = Wk * wdot[is];
	       }
        } // end nq loop
    } // end ne loop
} // end get_temperature_interface

} // end extern "C"
