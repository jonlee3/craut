#include "pch.h"
#include "Linker.h"


Linker::Linker(Helix h1, Helix h2, Angle angle) :
	h1(h1), h2(h2), angle(angle)
{
}

Helix Linker::get_h1()
{
	return h1;
}

Helix Linker::get_h2()
{
	return h2;
}

Angle Linker::get_angle()
{
	return angle;
}

std::string Linker::get_seq() {
	return h1.seq + angle.seq + h2.seq;
}