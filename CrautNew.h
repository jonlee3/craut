#pragma once
#include "Craut.h"

struct Helix;
struct Angle;

/**
Extends the <Craut> class by providing weighting functionality
for linkers. The weighting scheme is a simplified version of that
found in the Heidelberg IGEM team's code. The class is labeled as
"New", because it is an updated version of their weighting scheme.

*/
class CrautNew :
	public Craut
{
public:
	CrautNew(std::string pdb_filepath, char subunit, std::vector<Helix> helices, std::vector<Angle> angles);

protected:
	/**
	Returns the total length of the linker, divided by the termini distance.
	*/
	double length_weight(Linker linker);

	/**
	Returns a weight for a linker based on the angle in between the helices.
	The computation is based on a gaussian distribution.
	*/
	double angle_weight(Linker linker);

	/**
	Returns the distance from the linker to the protein. This is considered
	as the distance from the point of the linker other than the termini.
	*/
	double dist_weight(Linker linker);

	/**
	The overridden function from the base class <Craut>
	Returns the weight of a linker, to be used for
	ranking linkers.
	*/
	double get_weight(Linker linker);

	/**
	Returns a formatted string with info on the separate weight contributions.
	*/
	std::string get_weighting_info(Linker linker);

};

