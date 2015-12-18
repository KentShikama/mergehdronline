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
	\file fix_chromatic_aberration.hh
	\brief definition of FixChromaticAberration function
*/




#ifndef __FIX_CHROMATIC_ABERRATION_HH__
#define __FIX_CHROMATIC_ABERRATION_HH__




#include "shared_array.hh"




//============================================================================
//    FixChromaticAberration function
//============================================================================


SharedArray< float > const FixChromaticAberration(
	SharedArray< float > const& image,
	unsigned int const rows,
	unsigned int const columns,
	int const degree = 3,
	double const sigma1 = 1,       // Canny smoothing for edge-detection
	double const sigma2 = 16,      // local normalization radius for edge-detection
	double const lambda = 0.1,     // Hessian regularization parameter for edge-detection
	double const alpha = 1,        // alpha parameter to Nelder-Mead
	double const gamma = 2,        // gamma parameter to Nelder-Mead
	double const rho   = 0.5,      // rho parameter to Nelder-Mead
	double const sigma = 0.5,      // sigma parameter to Nelder-Mead
	double const delta = 0.001,    // initial simplex size
	double const epsilon = 0.0001
);




#endif    /* __FIX_CHROMATIC_ABERRATION_HH__ */
