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
	\file image.hh
	\brief definition of image manipulation functions
*/




#ifndef __IMAGE_HH__
#define __IMAGE_HH__




#include "shared_array.hh"

#include <vector>

#include <stdint.h>




//============================================================================
//    SaveEXR function
//============================================================================


void SaveEXR(
	char const* const filename,
	SharedArray< float > const& image,
	unsigned int const rows,
	unsigned int const columns
);




//============================================================================
//    SaveNonlinearImage function
//============================================================================


void SaveNonlinearImage(
	char const* const filename,
	SharedArray< float > const& image,
	unsigned int const rows,
	unsigned int const columns,
	double const stops = 1
);




//============================================================================
//    LoadEXR function
//============================================================================


SharedArray< float > const LoadEXR(
	char const* const filename,
	unsigned int& rows,
	unsigned int& columns
);




//============================================================================
//    LoadLinearTIFF function
//============================================================================


SharedArray< uint16_t > const LoadLinearTIFF(
	char const* const filename,
	unsigned int& rows,
	unsigned int& columns
);




//============================================================================
//    LoadNonlinearImage function
//============================================================================


SharedArray< uint8_t > const LoadNonlinearImage(
	char const* const filename,
	unsigned int& rows,
	unsigned int& columns
);




//============================================================================
//    MergeLinearImages function
//============================================================================


SharedArray< float > MergeLinearImages(
	std::vector< SharedArray< uint16_t > > const& images,
	unsigned int const rows,
	unsigned int const columns,
	double const pixels    // actually a lower bound on the number of pixels
);




//============================================================================
//    MergeNonlinearImages function
//============================================================================


SharedArray< float > MergeNonlinearImages(
	std::vector< SharedArray< uint8_t > > const& images,
	unsigned int const rows,
	unsigned int const columns,
	double const pixels,    // actually a lower bound on the number of pixels
	double const stops,
	double const lambda,
	double const epsilon
);




#endif    /* __IMAGE_HH__ */
