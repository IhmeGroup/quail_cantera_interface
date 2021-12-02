#include "cantera/thermo.h"
#include <iostream>

extern "C" {

using namespace Cantera;

void set_mixture_TRY(std::unique_ptr<ThermoPhase> & gas, const double T, 
    const double rho, const double *Y);

double get_energy_from_temperature(std::unique_ptr<ThermoPhase> & gas, 
    const double T_in);

double get_temperature_from_energy(std::unique_ptr<ThermoPhase> & gas, 
    const double e_in);

} // end extern "C"