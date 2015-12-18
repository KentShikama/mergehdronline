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
	\file fix_noise.hh
	\brief definition of FixNoise function
*/




#ifndef __FIX_NOISE__
#define __FIX_NOISE__




#include "shared_array.hh"




//============================================================================
//    FixNoise function
//============================================================================


SharedArray< float > const FixNoise(
	SharedArray< float > const& image,
	unsigned int const rows,
	unsigned int const columns,
	double const pixels,
	int const radius,
	int const dimension,
	double const alpha,
	double const beta,
	int const window,
	double const epsilon
);




#endif    /* __FIX_NOISE_HH__ */
