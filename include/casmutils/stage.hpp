#ifndef STAGE_HH
#define STAGE_HH

#include "casmutils/definitions.hpp"
#include "casmutils/structure.hpp"
#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Structure.hh>
#include <iostream>


std::vector<std::vector<std::pair<CASM::Site,double>>> find_nearest_neighbors( const Rewrap::Structure &struc, double min_radius, double max_radius);

#endif
