#include <iostream>
#include <memory>
#include <omp.h>
// #include "helpers.h"
#include "transport.h"

extern "C" {

using namespace Cantera;

double get_viscosity_from_temperature(std::unique_ptr<ThermoPhase> & gas, 
    std::unique_ptr<Transport> & trans, 
    const double T_in){
    gas->setTemperature(T_in);
    return trans->viscosity();
} // end get_viscosity_from_temperature

double get_thermalconductivity_from_temperature(std::unique_ptr<ThermoPhase> & gas, 
    std::unique_ptr<Transport> & trans, 
    const double T_in){
    gas->setTemperature(T_in);
    return trans->thermalConductivity();
} // end get_thermalconductivity_from_temperature
} // end extern C
