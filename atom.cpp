#include "class.h"

void atom :: get_pos0(vec & B)
{
	pos0 = B;
}

void atom :: get_pos(vec & B)
{
	pos = B;
}

void atom :: get_dipole()
{
	dipole = pos - pos0;
}

atom & atom :: operator=(const atom& B)
{
	pos = B.pos;
	pos0 = B.pos0;
	dipole = B.dipole;
	return *this;
}
