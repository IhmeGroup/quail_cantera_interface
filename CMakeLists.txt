cmake_minimum_required(VERSION 3.9)
project (quail_cantera_interface)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 14)

find_package(Threads REQUIRED)
list(APPEND CMAKE_CXX_FLAGS "-fopenmp")
include_directories($ENV{CANTERA_ROOT}/include /home/ihme/bornhoft/utilities/sundials/sundials-3.1.1/include $ENV{EIGEN3_ROOT} $ENV{BOOST_ROOT}/include)
link_directories($ENV{CANTERA_ROOT}/lib /home/ihme/bornhoft/utilities/sundials/sundials-3.1.1/install/lib)

add_library(quail_cantera_interface SHARED
    src/interface.cpp
    src/thermo.cpp
    src/helpers.cpp
)

target_link_libraries(quail_cantera_interface cantera sundials_cvodes sundials_ida sundials_nvecserial)
