# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(header_path "${BSMPT_SOURCE_DIR}/include/BSMPT/baryo_fhck")
set(header
    ${header_path}/difeq_transport_equations.h
    ${header_path}/TransportModel.h
    ${header_path}/TransportEquations.h)

set(src
    TransportModel.cpp
    TransportEquations.cpp)

add_library(BaryoFHCK ${header} ${src})
target_link_libraries(
  BaryoFHCK
  PUBLIC Eigen3::Eigen
         GSL::gsl
         Minimizer
         Kfactors
         LibWallThickness
         Utility
         Baryo
         Boost::boost
         BounceSolution
         GW
         VacuumProfile)
# target_compile_options(BaryoFHCK PUBLIC ${OpenMP_CXX_FLAGS})
target_include_directories(BaryoFHCK PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(BaryoFHCK PUBLIC cxx_std_17)
