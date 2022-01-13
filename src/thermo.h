#include "cantera/thermo.h"
#include "cantera/IdealGasMix.h"



#include <iostream>

#define GAS_CONSTANT 8.3144621000e3 // [J / (K kmole)]

extern "C" {

using namespace Cantera;

double get_P(std::unique_ptr<ThermoPhase> & gas, const double rho, const double e, double *Y);

double get_T(std::unique_ptr<ThermoPhase> & gas, const double rho, const double e, double *Y);

double get_gamma(std::unique_ptr<ThermoPhase> & gas);

double get_energy_from_temperature(std::unique_ptr<ThermoPhase> & gas, 
    const double T_in);

double get_temperature_from_energy(std::unique_ptr<ThermoPhase> & gas, 
    const double e_in);

void set_mixture_TRY(std::unique_ptr<ThermoPhase> & gas, const double T, 
    const double rho, const double *Y);

} // end extern "C"