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
	\file hdr_denoise.cc
*/




#include "headers.hh"




//============================================================================
//    main function
//============================================================================


int main( int argc, char const* argv[] ) {

	int code = EXIT_SUCCESS;

	int help = 0;
	char* input  = NULL;
	char* output = NULL;
	double pixels = std::numeric_limits< double >::quiet_NaN();
	int radius = 4;
	int dimension = 32;
	double alpha = std::numeric_limits< double >::quiet_NaN();
	double beta  = std::numeric_limits< double >::quiet_NaN();
	int window = 63;
	double epsilon = std::numeric_limits< double >::quiet_NaN();

	struct poptOption options[] = {
		{ "help",      'h',  POPT_ARG_NONE,   &help,      0, "show help",                      NULL            },
		{ "input",     'i',  POPT_ARG_STRING, &input,     0, "input file (EXR)",               "INPUT"         },
		{ "output",    'o',  POPT_ARG_STRING, &output,    0, "output file (EXR)",              "OUTPUT"        },
		{ "pixels",    'p',  POPT_ARG_DOUBLE, &pixels,    0, "how many pixels should we use?", "FLOAT(=10000)" },
		{ "radius",    'r',  POPT_ARG_INT,    &radius,    0, "descriptor radius",              "INT(=4)"       },
		{ "dimension", 'd',  POPT_ARG_INT,    &dimension, 0, "descriptor dimension",           "INT(=32)"      },
		{ "alpha",     'a',  POPT_ARG_DOUBLE, &alpha,     0, "alpha parameter",                "FLOAT(=0.2)"   },
		{ "beta",      'b',  POPT_ARG_DOUBLE, &beta,      0, "beta parameter",                 "FLOAT(=0.5)"   },
		{ "window",    'w',  POPT_ARG_INT,    &window,    0, "window side length",             "INT(=63)"      },
		{ "epsilon",   'e',  POPT_ARG_DOUBLE, &epsilon,   0, "epsilon parameter",              "FLOAT(=1e-12)" },
		{ NULL,        '\0', 0,               NULL,       0, NULL,                             NULL            }
	};

	poptContext context = poptGetContext( NULL, argc, argv, options, 0 );

	try {

		for ( ; ; ) {

			int const code = poptGetNextOpt( context );
			if ( code == -1 )
				break;
			else if ( code < -1 )
				throw std::runtime_error( poptStrerror( code ) );
		}
		if ( poptGetArg( context ) != NULL )    // check for positional arguments
			throw std::runtime_error( "unknown option" );

		if ( help ) {

			std::cout <<
				"Fixes CCD noise in an HDR input file, in EXR format, and saves the result to a" << std::endl <<
				"new HDR file, also in EXR format." << std::endl <<
				std::endl <<
				"The denoising is performed using an algorithm related to:" << std::endl <<
				"\tJeff Orchard, Mehran Ebrahimi and Alexander Wong. \"Efficient" << std::endl <<
				"\tNonlocal-Means Denoising Using the SVD\". Proceedings of the IEEE" << std::endl <<
				"\tConference on Image Processing. Pages 1732-1735. 2008." << std::endl <<
				"This paper accelerates the original non-local means algorithm by restricting" << std::endl <<
				"the search for similar pixels to only a neighborhood of the pixel of interest," << std::endl <<
				"rather than the entire image. The window parameter is the size of this" << std::endl <<
				"neighborhood. The weight on each pixel is the Gaussian-like function" << std::endl <<
				"exp(-gamma*(d^2)), where d is the Tukey-window weighted L2 distance between" << std::endl <<
				"windows of size 2*radius+1 around the pixel of interest, and the candidate" << std::endl <<
				"pixel, both projected onto the PCA subspace of the specified dimension.  The" << std::endl <<
				"gamma parameter to this weight function is determined by first finding the" << std::endl <<
				"descriptor distances between many adjacent pixels (the number of samples is" << std::endl <<
				"determined by the pixels parameter), and then choosing gamma such that" << std::endl <<
				"differences of the alpha*100th percentile have the weight given by the beta" << std::endl <<
				"parameter. The PCA subspace is found using the subspace iteration algorithm," << std::endl <<
				"applied to a covariance matrix calculated from the specified number of pixels," << std::endl <<
				"with termination threshold epsilon." << std::endl <<
				std::endl <<
				"This program is extremely slow! Be patient." << std::endl <<
				std::endl
			;
			poptPrintHelp( context, stdout, 0 );
		}
		else {

			if ( input == NULL )
				throw std::runtime_error( "input file not specified" );
			if ( output == NULL )
				throw std::runtime_error( "output file not specified" );
			if ( std::isnan( pixels ) )
				pixels = 10000;
			if ( std::isnan( alpha ) )
				alpha = 0.2;
			if ( std::isnan( beta ) )
				beta = 0.5;
			if ( std::isnan( epsilon ) )
				epsilon = 1e-12;

			std::cout << "Loading HDR image" << std::endl;
			unsigned int rows, columns;
			SharedArray< float > image = LoadEXR( input, rows, columns );

			std::cout << "Denoising image" << std::endl;
			SharedArray< float > result = FixNoise(
				image,
				rows,
				columns,
				pixels,
				radius,
				dimension,
				alpha,
				beta,
				window,
				epsilon
			);

			std::cout << "Saving HDR image to " << output << std::endl;
			SaveEXR( output, result, rows, columns );
		}
	}
	catch( std::exception& exception ) {

		std::cerr << "Error: " << exception.what() << std::endl << std::endl;
		poptPrintHelp( context, stderr, 0 );
		code = EXIT_FAILURE;
	}

	return code;
}
