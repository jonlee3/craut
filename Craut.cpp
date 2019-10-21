#include "pch.h"
#include "Craut.h"
#include "utils.h"
#include "Linker.h"

Craut::Craut(std::string pdb_filepath, char subunit, std::vector<Helix> helices, std::vector<Angle> angles) :
	helices(helices), angles(angles)
{
	load_pdb(pdb_filepath, subunit);
	init_termini();
	compute_protein_centre();
	init_angle_separators();
	generate_linkers();
}

void Craut::init_termini() {
	// identify the first appearance of the atom "N"
	for (n_terminus_index = 0; n_terminus_index < (int)atoms.size(); n_terminus_index++)
		if (atom_names[n_terminus_index] == "N")
			break;

	// if there isn't such an atom, n_terminus_index will be atoms.size()
	if (n_terminus_index == (int)atoms.size())
		throw "The N-Terminus was not found.";

	// identify the last appearance of the atom "C"
	for (c_terminus_index = (int)atoms.size() - 1; c_terminus_index >= 0; c_terminus_index--)
		if (atom_names[c_terminus_index] == "C")
			break;

	// if there isn't such an atom, c_terminus_index will be -1
	if (c_terminus_index == -1)
		throw "The C-Terminus was not found.";

	// we can now assume the termini exist
}

void Craut::init_angle_separators() {
	// This code was translated from the Heidelberg IGEM team's code

	angle_separators.resize(angles.size() + 1, 0.);

	for (int i = 0; i < (int)angles.size() - 1; i++) {
		double factor = (angles[i + 1].mean - angles[i].mean) / (angles[i + 1].std_dev + angles[i].std_dev);
		angle_separators[i + 1] = factor * angles[i].std_dev + angles[i].mean;
	}

	*(--angle_separators.end()) = PI;
}

void Craut::generate_linkers() {
	// enumerate over all pairs of helices
	// allow repetitions, and assume order matters
	for (int i = 0; i < (int)helices.size(); i++) {
		for (int j = 0; j < (int)helices.size(); j++) {
			// triangle side lengths
			double a = helices[i].length, b = helices[j].length, c = termini_dist();
			
			if (triangle_ineq_satisfied(a, b, c)) {
				// get the angle that must be between the helices
				double angle = get_angle(a, b, c);

				// find which angle segment this corresponds to
				int angle_index = -1;
				for (int i = 0; i < (int)angle_separators.size(); i++) {
					if ((angle >= angle_separators[i]) & (angle <= angle_separators[i + 1])) {
						angle_index = i;
						break;
					}
				}
				
				// An angle segment should be found for any angle in [0, PI]
				if (angle_index == -1)
					throw "Invalid angle found";

				linkers.push_back(Linker(helices[i], helices[j], angles[angle_index]));
			}
		}
	}
}

void Craut::rank_linkers() {
	// sort the linkers using the weight function
	std::sort(linkers.begin(), linkers.end(),
		[this](Linker l, Linker r) -> bool {
		double a = get_weight(l);
		return get_weight(l) < get_weight(r);
	});
}

void Craut::rank_and_print_linkers(int num_to_print) {
	rank_linkers();
	for (int i = 0; i != num_to_print && i < (int)linkers.size(); i++)
		std::cout << linkers[i].get_seq() << " " << get_weighting_info(linkers[i]) << std::endl;
}

void Craut::load_pdb(std::string filepath, char subunit) {
	// see https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
	// for the pdb file format

	// clear in case we want to load another pdb
	atoms.clear();

	FILE* file;
	fopen_s(&file, filepath.c_str(), "rb");
	{
		char buffer[PDB_LINE_SIZE];

		// read first line
		int chars_read = fread(buffer, sizeof(char), PDB_LINE_SIZE, file);

		// read until end of file
		while (chars_read == PDB_LINE_SIZE) {
			// this line describes an atom iff the first 4 characters of it are "ATOM"
			// also check if this is part of the right subunit (column 21)
			if (strncmp(buffer, "ATOM", 4) == 0 && buffer[21] == subunit) {
				// store the atom coordinates
				// the x,y,z coordinates are in columns 30-53 (inclusive, 0-indexed), add 1 extra char for null terminator
				char extracted_coords[ATOM_COORD_LENGTH * 3 + 1];
				strncpy_s(extracted_coords, &buffer[30], ATOM_COORD_LENGTH * 3);
				atoms.push_back(string_to_vector3d(extracted_coords));

				// store the atom name
				char curr_atom_name[5];
				strncpy_s(curr_atom_name, &buffer[13], 4);
				atom_names.push_back(remove_whitespace(std::string(curr_atom_name)));
			}

			// read next line of the file
			chars_read = fread(buffer, sizeof(char), PDB_LINE_SIZE, file);
		}
	}
	fclose(file);
}

double Craut::termini_dist() {
	//assume n_terminus_index and c_terminus_index are valid indices into <atoms>
	return (atoms[n_terminus_index] - atoms[c_terminus_index]).norm();
}

std::string Craut::get_weighting_info(Linker linker)
{
	return std::to_string(get_weight(linker));
}

double Craut::get_linker_angle(Linker linker)
{
	double a = linker.get_h1().length, b = linker.get_h2().length, c = termini_dist();
	return get_angle(a,b,c);
}

void Craut::compute_protein_centre() {
	protein_centre = Eigen::Vector3d(0., 0., 0.);
	for (auto atom : atoms)
		protein_centre += atom;
	protein_centre /= atoms.size();
}

Eigen::Vector3d Craut::get_linker_tri_point(Linker linker) {
	Eigen::Vector3d start = atoms[n_terminus_index];
	Eigen::Vector3d end = atoms[n_terminus_index];

	// we have 2 points of the triangle, we need to compute the third

	// get the angle between h1 and the connection between the termini
	double angle = get_angle(linker.get_h1().length, termini_dist(), linker.get_h2().length);

	// height is the distance between the connection between the termini and the third traingle point
	double height = sin(angle)*termini_dist() / linker.get_h1().length;

	// get the projection of the third triangle point onto the connection between the termini
	Eigen::Vector3d proj = start + linker.get_h1().length * cos(angle) * (end - start).normalized();

	// this is the third point of the triangle formed by the linker (not the termini)
	return (proj - protein_centre).normalized() * height + proj;
}