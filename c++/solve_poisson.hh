/*
	Copyright (C) 2010  Andrew Cotter

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


/**
	\file solve_poisson.hh
	\brief definition of SolvePoisson function
*/




#ifndef __SOLVE_POISSON_HH__
#define __SOLVE_POISSON_HH__




#include "shared_array.hh"




//============================================================================
//    SolvePoisson function
//============================================================================


void SolvePoisson(
	SharedArray< float > const& solution,
	SharedArray< float > const& rhs,
	unsigned int const rows,
	unsigned int const columns,
	double const epsilon,    // termination criterion is this proportion of the L2 norm of rhs
	unsigned int const innerIterations = 16
);




#endif    /* __SOLVE_POISSON_HH__ */
