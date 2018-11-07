#include "casmutils/definitions.hpp"
#include "casmutils/handlers.hpp"
#include "casmutils/structure.hpp"
#include "casmutils/stage.hpp"
#include <boost/program_options.hpp>
#include <casm/crystallography/Structure.hh>
#include <fstream>
#include <iostream>

namespace Utilities
{

void nearest_neighbors_distance_initializer(po::options_description& nearest_neighbors_distance_desc)
{
    UtilityProgramOptions::add_help_suboption(nearest_neighbors_distance_desc);
//    UtilityProgramOptions::add_output_suboption(nearest_neighbors_distance_desc);

    nearest_neighbors_distance_desc.add_options()("structure,s", po::value<fs::path>()->required(),
                               "POS.vasp like file that you want to get the nearest neighbors for.");
    
    nearest_neighbors_distance_desc.add_options()("max_radius", po::value<double>()->required(),
                               "Maximum radius to look for neighbors in");

    nearest_neighbors_distance_desc.add_options()("min_radius", po::value<double>()->default_value(0),
                               "minimum radius to look for neighbors in");
    return;
}
} // namespace Utilities

using namespace Utilities;

int main(int argc, char* argv[])
{
    Handler nearest_neighbors_distance_launch(argc, argv, nearest_neighbors_distance_initializer);

    /* if(nearest_neighbors_distance_launch.count("help") || nearest_neighbors_distance_launch.argc()<2) */
    if (nearest_neighbors_distance_launch.count("help"))
    {
        std::cout << nearest_neighbors_distance_launch.desc() << std::endl;
        return 1;
    }

    try
    {
        nearest_neighbors_distance_launch.notify();
    }

    catch (po::required_option& e)
    {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    auto path = nearest_neighbors_distance_launch.fetch<fs::path>("structure");

    auto struc = Rewrap::Structure(path);
    auto max_radius= nearest_neighbors_distance_launch.fetch<double>("max_radius");
    auto min_radius=  nearest_neighbors_distance_launch.fetch<double>("min_radius");
//    if (nearest_neighbors_distance_launch.vm().count("output"))
//    {
//        auto out_path = nearest_neighbors_distance_launch.fetch<fs::path>("output");
//        Simplicity::write_poscar(prim_struc, out_path);
//    }

//    else
//    {
//    }

    auto neighbor_analysis=find_nearest_neighbors(struc,min_radius,max_radius);
    // access site i's neighbor by doing neighbor_analysis[i]
       for (int i =0 ; i < struc.basis.size() ; i++){
//             if struc.basis[i].occ_name()== "Sn"
       for ( auto pair : neighbor_analysis[i]){
                std::cout << pair.first.occ_name() << " is " << pair.second << " away from " << struc.basis[i].occ_name() << std::endl;
       }
       }

    return 0;
}
