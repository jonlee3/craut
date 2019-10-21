#pragma once

#define TO_RADIANS(x) x * PI / 180.

#include "utils.h"

/**
Represents a helix building block for a linker.
<seq> must be the sequence of animo acids
*/
struct Helix {
	std::string seq;
	double length;
	Helix() {}
	Helix(std::string seq, double length) : seq(seq), length(length) {}
};

/**
Represents an angle building block for a linker.
<seq> is the sequence of animo acids, <mean> is the mean
angle found with this sequence, and <std_dev> is the standard
deviation of the angles found with this sequence.
*/
struct Angle {
	std::string seq;
	double mean;
	double std_dev;
	Angle(double mean, double std_dev, std::string seq) : seq(seq), mean(TO_RADIANS(mean)), std_dev(TO_RADIANS(std_dev)) {}
};

/**
Represents a linker made up of 2 helices, and an angle
between the two helices.
*/
 class Linker
{
public:
	Linker(Helix h1, Helix h2, Angle angle);

	Helix get_h1();
	Helix get_h2();
	Angle get_angle();

	std::string get_seq();

private:
	Helix h1, h2;
	Angle angle;
};

