#include "pch.h"
#include "CrautNew.h"
#include "Linker.h"

CrautNew::CrautNew(std::string pdb_filepath, char subunit, std::vector<Helix> helices, std::vector<Angle> angles) :
	Craut(pdb_filepath, subunit, helices, angles)
{
}

double CrautNew::length_weight(Linker linker) {
	return (linker.get_h1().length + linker.get_h2().length) / termini_dist();
}

double CrautNew::angle_weight(Linker linker) {
	double sum = 0.;
	for (Angle a : angles)
		sum += gaussian(get_linker_angle(linker), a.std_dev, a.mean);

	// this will find the minimum standard deviation for all angle sequences, and will compute only once
	static double min_std_dev = std::min_element(angles.begin(), angles.end(), [](Angle l, Angle r) {return l.std_dev < r.std_dev; })->std_dev;

	//std::cout << " SUM:" << sum << " MIN_STD_DEV:" << min_std_dev << std::endl;

	return 2. / (min_std_dev * sqrt(2. * PI)) - sum;
}

double CrautNew::dist_weight(Linker linker) {
	Eigen::Vector3d tri_point = get_linker_tri_point(linker);

	// get min distance. we have at least 2 atoms for sure, since we have 2 termini
	// this could be optimized by using an axis-aligned bounding box tree
	double min_dist = (atoms[0] - tri_point).norm();
	for (int i = 1; i < (int)atoms.size(); i++)
		min_dist = std::min(min_dist, (atoms[i] - tri_point).norm());

	return min_dist;
}

double CrautNew::get_weight(Linker linker) {
	// these linear weights of the different weighting contributions
	// were taken from the Heidelberg IGEM team's code
	return 1.85707748e-06 * length_weight(linker) +
		0.57330412e+00 * angle_weight(linker) +
		dist_weight(linker);
}

std::string CrautNew::get_weighting_info(Linker linker){
	// the weights are computed a second time which is inefficient, but it's still relatively fast
	return "   weight=" + std::to_string(get_weight(linker)) +
		"   length_weight=" + std::to_string(length_weight(linker)) +
		"   angle_weight=" + std::to_string(angle_weight(linker)) +
		"   dist_weight=" + std::to_string(dist_weight(linker));
}
