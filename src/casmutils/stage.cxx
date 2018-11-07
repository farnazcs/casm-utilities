#include "casmutils/stage.hpp"

    typedef std::pair<CASM::Site,double> NeighborInfo;
    typedef std::vector<NeighborInfo> SiteNeighborInfo;

std::vector<SiteNeighborInfo> find_nearest_neighbors( const Rewrap::Structure &struc, double min_radius, double max_radius){
// need lattice for temporary point
    auto lattice = struc.lattice();
    std::vector<SiteNeighborInfo> neighbor_analysis;
    // The periodic shifts of the unit cell will be stored here
    CASM::Coordinate shift_value(lattice);

   // This should be read in from user input
   // dim stores a maximum range of the periodic boundary we
   // wish to explore
   auto dim = lattice.enclose_sphere(max_radius);
        // loop over each central site
       for(int i = 0; i < struc.basis.size(); i++) {
               CASM::Site tatom(struc.basis[i]);
         //get distance to closest basis site in the unit cell at the origin
   // this is a container holding the sites within max radius of site i
  SiteNeighborInfo neighbors_and_dists;
        //loop over all other sites to find the distance
         for(int j = 0; j < struc.basis.size(); j++) {
   // This gives all the periodic images of the unit cell
   // that we will be checking for neighbors in
   CASM::EigenCounter<Eigen::Vector3i > grid_count(-dim, dim, Eigen::Vector3i::Constant(1));
   // loop over periodic images last
     do{
             // construct the current periodic shift value
       shift_value.frac() = grid_count().cast<double>();
            // this is the jth basis site with shift = shift_value
            CASM::Site tatomj(struc.basis[j] + shift_value);
            // this is the distance from site i to the shifted site j
          double dist = tatom.dist(tatomj);
           if ((dist <= max_radius) && (dist >=min_radius) ) {
              NeighborInfo info = std::make_pair(tatomj,dist);
              neighbors_and_dists.push_back(info);
           }
          } while (++grid_count);
       }
	 //sort before push back
	auto mycomparefunction = [](const NeighborInfo& lhs, const NeighborInfo& rhs)->bool{
		if (lhs.second == rhs.second){
			return lhs.first.occ_name() < rhs.first.occ_name();
		}
		return lhs.second < rhs.second;
	};
	 std::sort (neighbors_and_dists.begin(),neighbors_and_dists.end(),mycomparefunction);
        neighbor_analysis.push_back(neighbors_and_dists);
               }
        return neighbor_analysis;
}

