#pragma once

/**
Given a c-string of length 3*ATOM_COORD_LENGTH containing the coordinates of a point,
returns the point as a Vector3f. Each point must be ATOM_COORD_LENGTH characters.
*/
Eigen::Vector3d string_to_vector3d(char* xyz_coords);


/**
Returns a string that is a copy of s but without any whitepsace characters
*/
std::string remove_whitespace(std::string s);


/**
Returns true iff the triangle inequality is sasified for a triangle
with side lengths a,b,c.
*/
bool triangle_ineq_satisfied(double a, double b, double c);

/**
Returns the angle C of a triangle with side lengths a,b,c.
Uses cosine law to find the angle
*/
double get_angle(double a, double b, double c);


/**
The probability density function of the gaussian distribution
*/
double gaussian(double x, double std_dev, double mean);

const double PI = 3.141592653589;