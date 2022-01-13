// #include "cantera/thermo.h"
#include <iostream>
#include <memory>
#include <omp.h>
#include "helpers.h"
#include "thermo.h"

extern "C" {

using namespace Cantera;

void set_mixture_TRY(std::unique_ptr<ThermoPhase> & gas, const double T, 
    const double rho, const double *Y) {
    gas->setState_TRY(T, rho, Y);
}


double get_P(std::unique_ptr<ThermoPhase> & gas, const double rho, const double e, double *Y){
    double T = get_T(gas, rho, e, Y);
    return rho * GAS_CONSTANT / gas->meanMolecularWeight() * T;
} // end get_P


double get_T(std::unique_ptr<ThermoPhase> & gas, const double rho, const double e, double *Y){
    const double T_guess = 300.0;
    set_mixture_TRY(gas, T_guess, rho, Y);
    return get_temperature_from_energy(gas, e);
} // end get_T

double get_gamma(std::unique_ptr<ThermoPhase> & gas){
    return gas->cp_mass() / gas->cv_mass();
} // end get_gamma


double get_energy_from_temperature(std::unique_ptr<ThermoPhase> & gas, const double T_in){
    gas->setTemperature(T_in);
    return gas->intEnergy_mass();
} // end get_energy_from_temperature


double get_temperature_from_energy(std::unique_ptr<ThermoPhase> & gas, 
    const double e_in) {

    constexpr double reltol = std::numeric_limits<double>::epsilon() * 10.0;
    constexpr double relaxed_reltol = 1e-6;
    constexpr unsigned maxiter = 100;

    double T_n = gas->temperature(), dT = 0.0;
    double f_n; // initial guess
    int count = 0;

    while (count++ < maxiter) {
        f_n = get_energy_from_temperature(gas, T_n) - e_in;
        dT = f_n / gas->cv_mass();
        T_n -= dT;
        if (fabs(dT) < T_n * reltol) {
            return T_n;
        }
    }
    if (abs(dT) < T_n * relaxed_reltol) {

        std::cout << "MultiSpecies::GetTemperatureFromEnergy: " <<std::endl;
        std::cout << "too many iterations, " << std::endl;
        std::cout << "last dT = " << dT << ", " << std::endl;
        std::cout << "last T_n = " << T_n << "K!\n"
            << std::endl;
      }
} // end get_temperature_from_energy

} // end extern C
