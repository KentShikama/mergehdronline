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
	\file hdr_squish.cc
*/




#include "headers.hh"




//============================================================================
//    FindLuminance function
//============================================================================


SharedArray< float > FindLuminance(
	SharedArray< float > const image,
	unsigned int const rows,
	unsigned int const columns
)
{
	SharedArray< float > result( rows * columns );

	#pragma omp parallel for schedule( static )
	for ( int ii = 0; ii < static_cast< int >( rows * columns ); ++ii ) {

		float const red   = image[ ii * 3 + 0 ];
		float const green = image[ ii * 3 + 1 ];
		float const blue  = image[ ii * 3 + 2 ];

		double const intensity = std::max( red + green + blue, 3.0f / 65536 ) / 3.0;

		result[ ii ] = std::log( intensity );
	}

	return result;
}




//============================================================================
//    ScaleLuminance function
//============================================================================


SharedArray< float > ScaleLuminance(
	SharedArray< float > const image,
	SharedArray< float > const luminance,
	unsigned int const rows,
	unsigned int const columns,
	double const saturation
)
{
	if ( saturation <= 0 )
		throw std::runtime_error( "saturation must be positive" );

	SharedArray< float > result( rows * columns * 3 );

	float maximum = -std::numeric_limits< float >::infinity();
	#pragma omp parallel
	{	float localMaximum = -std::numeric_limits< float >::infinity();

		#pragma omp for schedule( static )
		for ( int ii = 0; ii < static_cast< int >( rows * columns ); ++ii ) {

			float const red   = image[ ii * 3 + 0 ];
			float const green = image[ ii * 3 + 1 ];
			float const blue  = image[ ii * 3 + 2 ];

			double const intensity = std::max( red + green + blue, 3.0f / 65536 ) / 3.0;

			result[ ii * 3 + 0 ] = saturation * std::log( std::max( red,   1.0f / 65536 ) / intensity ) + luminance[ ii ];
			result[ ii * 3 + 1 ] = saturation * std::log( std::max( green, 1.0f / 65536 ) / intensity ) + luminance[ ii ];
			result[ ii * 3 + 2 ] = saturation * std::log( std::max( blue,  1.0f / 65536 ) / intensity ) + luminance[ ii ];

			localMaximum = std::max( localMaximum, result[ ii * 3 + 0 ] );
			localMaximum = std::max( localMaximum, result[ ii * 3 + 1 ] );
			localMaximum = std::max( localMaximum, result[ ii * 3 + 2 ] );
		}

		#pragma omp critical
		maximum = std::max( maximum, localMaximum );
	}

	#pragma omp parallel for schedule( static )
	for ( int ii = 0; ii < static_cast< int >( rows * columns * 3 ); ++ii )
		result[ ii ] = std::exp( result[ ii ] - maximum );

	return result;
}




//============================================================================
//    main function
//============================================================================


int main( int argc, char const* argv[] ) {

	int code = EXIT_SUCCESS;

	int help = 0;
	char* input  = NULL;
	char* output = NULL;
	double alpha      = std::numeric_limits< double >::quiet_NaN();
	double beta       = std::numeric_limits< double >::quiet_NaN();
	double delta      = std::numeric_limits< double >::quiet_NaN();
	double theta      = std::numeric_limits< double >::quiet_NaN();
	double epsilon    = std::numeric_limits< double >::quiet_NaN();
	double saturation = std::numeric_limits< double >::quiet_NaN();

	struct poptOption options[] = {
		{ "help",       'h',  POPT_ARG_NONE,   &help,       0, "show help",            NULL             },
		{ "input",      'i',  POPT_ARG_STRING, &input,      0, "input file (EXR)",     "INPUT"          },
		{ "output",     'o',  POPT_ARG_STRING, &output,     0, "output file (EXR)",    "OUTPUT"         },
		{ "alpha",      'a',  POPT_ARG_DOUBLE, &alpha,      0, "alpha parameter",      "FLOAT(=0.1)"    },
		{ "beta",       'b',  POPT_ARG_DOUBLE, &beta,       0, "beta parameter",       "FLOAT(=0.1)"    },
		{ "delta",      'd',  POPT_ARG_DOUBLE, &delta,      0, "delta parameter",      "FLOAT(=1.1)"    },
		{ "theta",      't',  POPT_ARG_DOUBLE, &theta,      0, "theta parameter",      "FLOAT(=0)"      },
		{ "epsilon",    'e',  POPT_ARG_DOUBLE, &epsilon,    0, "epsilon parameter",    "FLOAT(=0.0001)" },
		{ "saturation", 's',  POPT_ARG_DOUBLE, &saturation, 0, "saturation parameter", "FLOAT(=1)"      },
		{ NULL,         '\0', 0,               NULL,        0, NULL,                   NULL             }
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
				"Squishes the contrast in an input HDR image, and saves the result. Both the" << std::endl <<
				"input and the output images are in EXR format. A great deal of tweaking of the" << std::endl <<
				"\"alpha\", \"beta\", \"delta\" and \"theta\" arguments is usually required before one" << std::endl <<
				"gets pleasing results." << std::endl <<
				std::endl <<
				"The HDR contrast processing is performed using an implementation of:" << std::endl <<
				"\tRaanan Fattal, Dani Lischinski and Michael Werman. \"Gradient Domain" << std::endl <<
				"\tHigh Dynamic Range Compression\". Proceedings of the 29th annual" << std::endl <<
				"\tconference on computer graphics and interactive techniques. Pages" << std::endl <<
				"\t249-256. 2002." << std::endl <<
				"The alpha, beta and saturation parameters are analagous to those described in" << std::endl <<
				"the paper, except that in this implementation, alpha is a proportion of the" << std::endl <<
				"average gradient magnitude, and beta is what would be 1-beta in the paper." << std::endl <<
				"Gamma is the maximum gradient attenuation factor at a given level, which is" << std::endl <<
				"needed in order to cope with small gradients. The paper isn't clear on whether" << std::endl <<
				"the average gradient magnitudes used for the alpha scaling should be the" << std::endl <<
				"overall averages, or the average only over the current level in the pyramid. We" << std::endl <<
				"use a convex combination of these two averages, mediated by the theta" << std::endl <<
				"parameter, where 0 chooses the overall average, and 1 the average on the" << std::endl <<
				"current pyramid level." << std::endl <<
				std::endl <<
				"Essentially, alpha controls the point at which contrast-magnification gives way" << std::endl <<
				"to contrast-attenuation, beta controls the overall extent to which" << std::endl <<
				"contrast-squishing is performed, and increasing theta causes large-scale" << std::endl <<
				"(spatially) contrast differences to be given greater importance, relative to" << std::endl <<
				"small-scale differences. Changing delta seldom has any great effect. Epsilon is" << std::endl <<
				"the termination threshold--increasing it will reduce the quality of the" << std::endl <<
				"solution, but improve runtime, while decreasing it will have the opposite" << std::endl <<
				"effect. One warning: if epsilon is too small, then the solution may never" << std::endl <<
				"converge, and you'll need to abort this program." << std::endl <<
				std::endl
			;
			poptPrintHelp( context, stdout, 0 );
		}
		else {

			if ( input == NULL )
				throw std::runtime_error( "input file not specified" );
			if ( output == NULL )
				throw std::runtime_error( "output file not specified" );
			if ( std::isnan( alpha ) )
				alpha = 0.1;
			if ( std::isnan( beta ) )
				beta = 0.1;
			if ( std::isnan( delta ) )
				delta = 1.1;
			if ( std::isnan( theta ) )
				theta = 0;
			if ( std::isnan( epsilon ) )
				epsilon = 0.0001;
			if ( std::isnan( saturation ) )
				saturation = 1;

			std::cout << "Loading HDR image" << std::endl;
			unsigned int rows, columns;
			SharedArray< float > image = LoadEXR( input, rows, columns );

			std::cout << "Finding luminance" << std::endl;
			SharedArray< float > const luminance = FindLuminance( image, rows, columns );

			std::cout << "Performing HDR luminance squishing" << std::endl;
			SquishLuminance( luminance, rows, columns, epsilon, alpha, beta, delta, theta );

			std::cout << "Scaling luminance" << std::endl;
			SharedArray< float > const result = ScaleLuminance( image, luminance, rows, columns, saturation );

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
