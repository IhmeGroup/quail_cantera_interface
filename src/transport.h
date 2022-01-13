#include "cantera/thermo.h"
#include "cantera/transport.h"
#include <iostream>

extern "C" {

using namespace Cantera;

double get_viscosity_from_temperature(std::unique_ptr<ThermoPhase> & gas, 
	std::unique_ptr<Transport> & trans, 
	const double T_in);

double get_thermalconductivity_from_temperature(std::unique_ptr<ThermoPhase> & gas, 
    std::unique_ptr<Transport> & trans, 
    const double T_in);

} // end extern "C"