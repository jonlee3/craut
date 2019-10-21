#include "pch.h"
#include "utils.h"

Eigen::Vector3d string_to_vector3d(char* xyz_coords) {
	//vector to be returned
	Eigen::Vector3d v(0., 0., 0.);

	//position of the next coord to be read in the string
	char* next_double_pos = xyz_coords;

	for (int i = 0; i < 3; i++)
		v[i] = strtod(next_double_pos, &next_double_pos);

	return v;
}

std::string remove_whitespace(std::string s){
	s.erase(std::remove_if(s.begin(), s.end(), isspace), s.end());
	return s;
}

bool triangle_ineq_satisfied(double a, double b, double c) {
	double lengths[3] = { a, b, c };

	//check if the triangle inequality is satisfied for all cases
	for (int i = 0; i < 3; i++)
		if (lengths[i] + lengths[(i + 1) % 3] < lengths[(i + 2) % 3])
			return false;

	return true;
}

double get_angle(double a, double b, double c)
{
	return acos((a * a + b * b - c * c) / (2. * a * b));
}

double gaussian(double x, double std_dev, double mean) {
	double temp = (x - mean) / std_dev;
	return 1. / (std_dev*sqrt(2. * PI)) * exp(-0.5 * temp * temp);
}
