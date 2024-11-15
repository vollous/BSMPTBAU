# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(header_path "${BSMPT_SOURCE_DIR}/include/BSMPT/baryo_fhck")
set(header
    ${header_path}/TransportEquations.h
    ${header_path}/BubbleWallProfile.h)

set(src
    TransportEquations.cpp
    BubbleWallProfile.cpp)

add_library(BaryoFHCK ${header} ${src})
target_link_libraries(
  BaryoFHCK
  PUBLIC Eigen3::Eigen
         GSL::gsl
         Minimizer
         Kfactors
         LibWallThickness
         Utility
         Boost::boost
         BounceSolution
         GW)
# target_compile_options(BaryoFHCK PUBLIC ${OpenMP_CXX_FLAGS})
target_include_directories(BaryoFHCK PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(BaryoFHCK PUBLIC cxx_std_17)
