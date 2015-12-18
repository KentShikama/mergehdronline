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
	\file headers.hh
	\brief Precompiled header
*/




#ifndef __HEADERS_HH__
#define __HEADERS_HH__




#include <popt.h>
#include <tiff.h>
#include <tiffio.h>
#include <ImfRgbaFile.h>
#include <Magick++.h>


#include "squish_luminance.hh"
#include "solve_poisson.hh"
#include "fix_chromatic_aberration.hh"
#include "fix_noise.hh"
#include "image.hh"
#include "shared_array.hh"
#include "helpers.hh"


#include <vector>
#include <string>

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <memory>
#include <limits>

#include <cstdlib>
#include <cstddef>
#include <cmath>


#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>




#endif    /* __HEADERS_HH__ */
