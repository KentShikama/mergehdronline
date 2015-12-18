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
	\file hdr_fix_ca.cc
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
	int degree = 1;
	double epsilon = std::numeric_limits< double >::quiet_NaN();

	struct poptOption options[] = {
		{ "help",    'h',  POPT_ARG_NONE,   &help,    0, "show help",         NULL             },
		{ "input",   'i',  POPT_ARG_STRING, &input,   0, "input file (EXR)",  "INPUT"          },
		{ "output",  'o',  POPT_ARG_STRING, &output,  0, "output file (EXR)", "OUTPUT"         },
		{ "degree",  'd',  POPT_ARG_INT,    &degree,  0, "polynomial degree", "INT(=1)"        },
		{ "epsilon", 'e',  POPT_ARG_DOUBLE, &epsilon, 0, "epsilon parameter", "FLOAT(=0.0001)" },
		{ NULL,      '\0', 0,               NULL,     0, NULL,                NULL             }
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
				"Fixes chromatic aberration in an HDR input file, in EXR format, and saves the" << std::endl <<
				"result to a new HDR file, also in EXR format." << std::endl <<
				std::endl <<
				"The algorithm is somewhat different from, but strongly inspired by:" << std::endl <<
				"\tSing Bing Kang. \"Automatic Removal of Chromatic Aberration from a" << std::endl <<
				"\tSingle Image\". Proceedings of the IEEE Conference on Computer Vision" << std::endl <<
				"\tand Pattern Recognition. Pages 1-8. 2007." << std::endl <<
				"The cited paper corrects for a number of different causes of chromatic" << std::endl <<
				"aberration, but this implementation only attempts to scale the color channels," << std::endl <<
				"in proportion to a polynomial (in the radius from the center of the image)," << std::endl <<
				"where the degree parameter is the degree of this polynomial. The coefficients" << std::endl <<
				"are chosen to minimize the L2 distance between locally-normalized Canny edge" << std::endl <<
				"strengths in each channel. Following Kang, green is the reference channel, and" << std::endl <<
				"the Nelder-Mead simplex algorithm is used for optimization. The problem which" << std::endl <<
				"we are attmpeting to minimize is highly non-convex, with a large number of" << std::endl <<
				"local minima. Hence, this program can often take an extremely long time to" << std::endl <<
				"converge, and is not guaranteed to find a good solution. As usual, the epsilon" << std::endl <<
				"parameter is the termination threshold, and it may be increased to improve" << std::endl <<
				"runtime, or decreased to improve quality. In practice, the default value" << std::endl <<
				"consistently improves the chromatic aberration in the image, although perhaps" << std::endl <<
				"not quite as much as would be desired." << std::endl <<
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
			if ( std::isnan( epsilon ) )
				epsilon = 0.0001;

			std::cout << "Loading HDR image" << std::endl;
			unsigned int rows, columns;
			SharedArray< float > image = LoadEXR( input, rows, columns );

			std::cout << "Fixing chromatic aberration" << std::endl;
			SharedArray< float > result = FixChromaticAberration(
				image,
				rows,
				columns,
				degree,
				1,
				16,
				0.1,
				1,
				2,
				0.5,
				0.5,
				0.001,
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
