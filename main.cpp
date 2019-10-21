#include "pch.h"
#include "CrautNew.h"
#include "utils.h"
#include "Linker.h"

static std::vector<Helix> helices = { Helix("AEAAAK", 8.7), Helix("AEAAAKA", 10.), Helix("AEAAAKAA", 10.8),
									Helix("AEAAAKEAAAK", 15.6), Helix("AEAAAKEAAAKA", 16.8), Helix("AEAAAKEAAAKEAAAKA", 24.8),
									Helix("AEAAAKEAAAKEAAAKEAAAKA", 32.3), Helix("AEAAAKEAAAKEAAAKEAAAKEAAAKA", 40.5) };

static std::vector<Angle> angles = { Angle(29.7, 8.5, "NVL"), Angle(38.7, 30., "KTA"), Angle(60., 12., "AADGTL"),
									Angle(74.5, 27., "VNLTA"), Angle(117., 12., "AAAHPEA"), Angle(140., 15., "ASLPAA"),
									Angle(160., 5., "ATGDLA") };

int main(int argc, char** argv) {
	// TODO: add more general error checking for input
	std::string usage = "Invalid arguments. Usage: CRAUT.exe pdb_filepath [number of linkers to print]";
	if (argc < 2)
		std::cout << usage << std::endl;
	//CrautNew craut("C:\\Users\\Jon\\Desktop\\CRAUT\\5xjh.pdb", 'A', helices, angles);
	CrautNew craut(argv[1], 'A', helices, angles);
	if(argc < 3)
		craut.rank_and_print_linkers();
	else
		craut.rank_and_print_linkers(atoi(argv[2]));
}