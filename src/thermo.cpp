#include "cantera/thermo.h"
#include "cantera/IdealGasMix.h"
#include <iostream>
#include <memory>
#include <omp.h>
#include "helpers.h"

extern "C" {

using namespace Cantera;

void get_massfractions(const int nsp, const double rho, const double *rhoY, double *Y, std::string& Ys,
    const std::vector<std::string>& Ynames)
{
    double Y_last = 1.0;
    for (int is = 0; is < nsp - 1; is++){
        Y[is] = rhoY[is] / rho;
        Y_last += -Y[is];

        Ys.append(Ynames[is] + ": " + std::to_string(Y[is]) + ", ");
    }
    Y[nsp - 1] = Y_last;
    Ys.append(Ynames[nsp - 1] + ": " + std::to_string(Y_last));
}

void get_pressure_update(const double* Uq, double* P, int ne, int nq, 
    int ns, int nsp, int dim, char* filename)
{
double T_final; // converged Temperature
// Create a new phase
std::unique_ptr<ThermoPhase> gas(newPhase(filename));
// IdealGasMix gas(filename);
// #pragma omp parallel for schedule(static, 1)
    for (int ie = 0; ie < ne; ie++){
        for (int iq = 0; iq < nq; iq++){

            auto start = ie * (nq*ns) + iq * ns;
            const double *U = Uq + start;
            const double rho = U[0];
            const double rhou = U[1];
            const double rhoE = U[2];
            double T_guess = 350.0;
            auto e = rhoE / rho - 0.5 * rhou * rhou / rho;

            double Y[nsp];
            std::string Y_string;
            get_massfractions(nsp, rho, U + dim + 2, Y, Y_string, gas->speciesNames());
            set_mixture_TRY(gas, T_guess, rho, Y);

            T_final = get_temperature_from_energy(gas, e);
            gas->setTemperature(T_final);

            // gas->setMassFractionsByName(Y_string);
            // gas->setState_UV(e, nu);
            P[ie * (nq) + iq] = gas->pressure();
        }
    }
} // end get_pressure



void get_pressure(const double* Uq, double* P, int ne, int nq, 
    int ns, int nsp, int dim, char* filename)
{
double T_final; // converged Temperature
// Create a new phase
std::unique_ptr<ThermoPhase> gas(newPhase(filename));
// IdealGasMix gas(filename);
// #pragma omp parallel for schedule(static, 1)
    for (int ie = 0; ie < ne; ie++){
        for (int iq = 0; iq < nq; iq++){

            auto start = ie * (nq*ns) + iq * ns;
            const double *U = Uq + start;
            const double rho = U[0];
            const double rhou = U[1];
            const double rhoE = U[2];
            double T_guess = 350.0;
            auto e = rhoE / rho - 0.5 * rhou * rhou / rho;

            double Y[nsp];
            std::string Y_string;
            get_massfractions(nsp, rho, U + dim + 2, Y, Y_string, gas->speciesNames());
            set_mixture_TRY(gas, T_guess, rho, Y);

            T_final = get_temperature_from_energy(gas, e);
            gas->setTemperature(T_final);

            // gas->setMassFractionsByName(Y_string);
            // gas->setState_UV(e, nu);
            P[ie * (nq) + iq] = gas->pressure();
        }
    }
} // end get_pressure


void get_specificheatratio(const double* Uq, double* gam, int ne, int nq, 
    int ns, int nsp, int dim, char* filename)
{

// int nThreads = omp_get_max_threads();
// writelog("Running on {} threads\n\n", nThreads);
// Create a new phase
std::unique_ptr<ThermoPhase> gas(newPhase(filename));

// #pragma omp parallel for schedule(static, 1)
    for (int ie = 0; ie < ne; ie++){
        for (int iq = 0; iq < nq; iq++){

            auto start = ie * (nq*ns) + iq * ns;
            const double *U = Uq + start;
            const double rho = U[0];
            const double rhou = U[1];
            const double rhoE = U[2];
            const double rhoY0 = U[3];

            auto e = rhoE / rho - 0.5 * rhou * rhou / rho;
            auto nu = 1.0 / rho;
            
            double Y[nsp];
            std::string Y_string;
            get_massfractions(nsp, rho, U + dim + 2, Y, Y_string, gas->speciesNames());
            gas->setMassFractionsByName(Y_string);
            gas->setState_UV(e, nu);
            gam[ie * (nq) + iq] = gas->cp_mass() / gas->cv_mass();
        }
    }
} // end get_specificheatratio


void get_temperature(const double* Uq, double* T, int ne, int nq, 
    int ns, int nsp, int dim, char* filename)
{

// int nThreads = omp_get_max_threads();
// writelog("Running on {} threads\n\n", nThreads);
// Create a new phase
std::unique_ptr<ThermoPhase> gas(newPhase(filename));

// #pragma omp parallel for schedule(static, 1)
    for (int ie = 0; ie < ne; ie++){
        for (int iq = 0; iq < nq; iq++){

            auto start = ie * (nq*ns) + iq * ns;
            const double *U = Uq + start;
            const double rho = U[0];
            const double rhou = U[1];
            const double rhoE = U[2];
            const double rhoY0 = U[3];

            auto e = rhoE / rho - 0.5 * rhou * rhou / rho;
            auto nu = 1.0 / rho;
            
            double Y[nsp];
            std::string Y_string;
            get_massfractions(nsp, rho, U + dim + 2, Y, Y_string, gas->speciesNames());
            gas->setMassFractionsByName(Y_string);
            gas->setState_UV(e, nu);
            T[ie * (nq) + iq] = gas->temperature();
        }
    }
} // end get_temperature


// NOTE: Not currently used -> actually not faster
// void get_pressure_multithread(const double* Uq, double* P, int ne, int nq, int ns)
// {

// // int nThreads = omp_get_max_threads();
// // writelog("Running on {} threads\n\n", nThreads);
// std::vector<std::unique_ptr<IdealGasMix>> gases;//(newPhase("air_test.xml"));

// std::unique_ptr<IdealGasMix> gas;
// gas = std::make_unique<IdealGasMix>("air_test.xml");
// // gases.resize(nThreads);

// // #pragma omp parallel
// // {
// //     size_t j = omp_get_thread_num();
// //     gases[j] = std::make_unique<IdealGasMix>("air_test.xml");
// // }
// // for (int i = 0; i < nThreads; i++) {
//     // gases.emplace_back(new IdealGasMix("air_test.xml"));
// // }
// // Doesn't currently work
// #pragma omp parallel
// {
// #pragma omp for nowait //schedule(static, 1)
// for (int ie = 0; ie < ne; ie++){
//     for (int iq = 0; iq < nq; iq++){

//             // Create a new phase
//             // std::unique_ptr<ThermoPhase> gas(newPhase("air_test.xml"));
//             // size_t j = omp_get_thread_num();
//             // IdealGasMix& gas = *gases[j];
//             // std::cout << j << std::endl;
//             auto start = ie * (nq*ns) + iq * ns;
//             const double *U = Uq + start;
//             const double rho = U[0];
//             const double rhou = U[1];
//             const double rhoE = U[2];
//             const double rhoY0 = U[3];
//             const double rhoY1 = U[4];

//             auto e = rhoE / rho - 0.5 * rhou * rhou / rho;
//             auto nu = 1.0 / rho;
//             auto YO2 = rhoY0 / rho;
//             auto YN2 = rhoY1 / rho;

//             std::string name = "N2:" + std::to_string(YN2) +", O2:" + std::to_string(YO2);
//             gas->setMassFractionsByName(name);
//             gas->setState_UV(e, nu);
//             P[ie * (nq) + iq] = gas->pressure();
//         }
//     }

// }
// }


// void get_pressure_multithread2(const double* Uq, double* P, int ne, int nq, int ns)
// {

// int nThreads = omp_get_max_threads();
// // writelog("Running on {} threads\n\n", nThreads);
// std::vector<std::unique_ptr<IdealGasMix>> gases;//(newPhase("air_test.xml"));
// gases.resize(nThreads);

// #pragma omp parallel
// {
//     size_t j = omp_get_thread_num();
//     gases[j] = std::make_unique<IdealGasMix>("air_test.xml");
// }
// // for (int i = 0; i < nThreads; i++) {
//     // gases.emplace_back(new IdealGasMix("air_test.xml"));
// // }


// #pragma omp parallel for //schedule(static, 1)
// for (int ie = 0; ie < ne; ie++){
//     for (int iq = 0; iq < nq; iq++){

//             size_t j = omp_get_thread_num();
//             IdealGasMix& gas = *gases[j];
//             auto start = ie * (nq*ns) + iq * ns;
//             const double *U = Uq + start;
//             const double rho = U[0];
//             const double rhou = U[1];
//             const double rhoE = U[2];
//             const double rhoY0 = U[3];
//             const double rhoY1 = U[4];

//             auto e = rhoE / rho - 0.5 * rhou * rhou / rho;
//             auto nu = 1.0 / rho;
//             auto YO2 = rhoY0 / rho;
//             auto YN2 = rhoY1 / rho;

//             std::string name = "N2:" + std::to_string(YN2) +", O2:" + std::to_string(YO2);
//             gas.setMassFractionsByName(name);
//             gas.setState_UV(e, nu);
//             P[ie * (nq) + iq] = gas.pressure();
//         }
//     }

// }


// void get_pressure(const double* Uq, double* P, int ne, int nq, 
//     int ns, int nsp, int dim, char* filename)
// {

// // Create a new phase
// std::unique_ptr<ThermoPhase> gas(newPhase(filename));
// // IdealGasMix gas(filename);
// // #pragma omp parallel for schedule(static, 1)
//     for (int ie = 0; ie < ne; ie++){
//         for (int iq = 0; iq < nq; iq++){

//             auto start = ie * (nq*ns) + iq * ns;
//             const double *U = Uq + start;
//             const double rho = U[0];
//             const double rhou = U[1];
//             const double rhoE = U[2];

//             auto e = rhoE / rho - 0.5 * rhou * rhou / rho;
//             auto nu = 1.0 / rho;

//             double Y[nsp];
//             std::string Y_string;
//             get_massfractions(nsp, rho, U + dim + 2, Y, Y_string, gas->speciesNames());
//             gas->setMassFractionsByName(Y_string);
//             gas->setState_UV(e, nu);
//             P[ie * (nq) + iq] = gas->pressure();
//         }
//     }
// } // end get_pressure




} // end extern C
