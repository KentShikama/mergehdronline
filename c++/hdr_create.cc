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
	\file hdr_create.cc
*/




#include "headers.hh"




//============================================================================
//    main function
//============================================================================


int main( int argc, char const* argv[] ) {

	int code = EXIT_SUCCESS;

	int help = 0;
	int linear = 0;
	char* output = NULL;
	double pixels  = std::numeric_limits< double >::quiet_NaN();
	double stops   = std::numeric_limits< double >::quiet_NaN();
	double lambda  = std::numeric_limits< double >::quiet_NaN();
	double epsilon = std::numeric_limits< double >::quiet_NaN();

	struct poptOption options[] = {
		{ "help",   'h',  POPT_ARG_NONE,   &help,    0, "show help",                                NULL             },
		{ "linear", '4',  POPT_ARG_NONE,   &linear,  0, "are the input files 16-bit linear TIFFs?", NULL             },
		{ "output", 'o',  POPT_ARG_STRING, &output,  0, "output file (EXR)",                        "OUTPUT"         },
		{ "pixels", 'p',  POPT_ARG_DOUBLE, &pixels,  0, "how many pixels should we use?",           "FLOAT(=10000)"  },
		{ "stops",  's',  POPT_ARG_DOUBLE, &stops,   0, "stops spanned by images",                  "FLOAT"          },
		{ "lambda", 'l',  POPT_ARG_DOUBLE, &lambda,  0, "lambda parameter",                         "FLOAT(=1)"      },
		{ "epsilon", 'e', POPT_ARG_DOUBLE, &epsilon, 0, "epsilon parameter",                        "FLOAT(=0.0001)" },
		{ NULL,     '\0', 0,               NULL,     0, NULL,                                       NULL             }
	};

	poptContext context = poptGetContext( NULL, argc, argv, options, 0 );
	poptSetOtherOptionHelp( context, "[OPTION...] [INPUTS...]" );

	try {

		for ( ; ; ) {

			int const code = poptGetNextOpt( context );
			if ( code == -1 )
				break;
			else if ( code < -1 )
				throw std::runtime_error( poptStrerror( code ) );
		}

		if ( help ) {

			std::cout <<
				"Creates and saves an HDR image, in EXR format, from a collection of images. If" << std::endl <<
				"the input images are 16-bit linear TIFFs (as created by dcraw), use the" << std::endl <<
				"\"linear\" option." << std::endl <<
				std::endl <<
				"Nonlinear image files are merged using a modified version of:" << std::endl <<
				"\tPaul Debevec and Jitendra Malik. \"Recovering High Dynamic Range" << std::endl <<
				"\tRadiance Maps from Photographs\". Proceedings of the 24th annual" << std::endl <<
				"\tconference on computer graphics and interactive techniques. Pages" << std::endl <<
				"\t369-378. 1997." << std::endl <<
				"Linear TIFFs are merged using a different (but related) algorithm, for which" << std::endl <<
				"the \"stops\", \"lambda\" and \"epsilon\" arguments are ignored." << std::endl <<
				std::endl <<
				"The stops parameter is self-explanatory. Lambda controls the smoothness of the" << std::endl <<
				"response curves which are recovered. This term is scaled in proportion to the" << std::endl <<
				"number of pixels which are used in determining these response curves and" << std::endl <<
				"exposure settings (this is the \"pixels\" parameter), so when using more pixels," << std::endl <<
				"it may be advisable to decrease lambda. Finally, epsilon is the termination" << std::endl <<
				"threshold. If you're getting unsatisfactory results, and suspect that the" << std::endl <<
				"reason for this is that a poor-quality solution is being found, then decreasing" << std::endl <<
				"epsilon (or increasing the number of pixels) might help. Conversely, if you are" << std::endl <<
				"satisfied with the solution, but want the program to run more quickly, then" << std::endl <<
				"increasing epsilon (or decreasing the number of pixels) could be a good idea." << std::endl <<
				"One warning: if epsilon is too small, then the solution may never converge, and" << std::endl <<
				"you'll need to abort this program." << std::endl <<
				std::endl
			;
			poptPrintHelp( context, stdout, 0 );
		}
		else {

			std::vector< std::string > inputs;
			for ( ; ; ) {

				char const* const input = poptGetArg( context );
				if ( input == NULL )
					break;
				inputs.push_back( input );
			}
			if ( inputs.empty() )
				throw std::runtime_error( "input file(s) not specified" );
			if ( output == NULL )
				throw std::runtime_error( "output file not specified" );

			unsigned int rows, columns;
			SharedArray< float > image;
			if ( linear ) {

				if ( inputs.size() > 1 ) {

					if ( std::isnan( pixels ) )
						pixels = 10000;
					if ( ! std::isnan( stops ) )
						std::cerr << "Warning: stops is ignored when merging 16-bit linear TIFFs" << std::endl;
					if ( ! std::isnan( lambda ) )
						std::cerr << "Warning: lambda is ignored when merging 16-bit linear TIFFs" << std::endl;
					if ( ! std::isnan( epsilon ) )
						std::cerr << "Warning: epsilon is ignored when merging 16-bit linear TIFFs" << std::endl;

					std::cout << "Loading images" << std::endl;
					std::vector< SharedArray< uint16_t > > inputImages;
					inputImages.push_back( LoadLinearTIFF( inputs[ 0 ].c_str(), rows, columns ) );
					for ( unsigned int ii = 1; ii < inputs.size(); ++ii ) {

						unsigned int newRows, newColumns;
						inputImages.push_back( LoadLinearTIFF( inputs[ ii ].c_str(), newRows, newColumns ) );
						if ( ( newRows != rows ) || ( newColumns != columns ) )
							throw std::runtime_error( "all input images must have the same dimensions" );
					}

					std::cout << "Merging multiple exposures" << std::endl;
					image = MergeLinearImages( inputImages, rows, columns, pixels );
				}
				else {

					if ( ! std::isnan( pixels ) )
						std::cerr << "Warning: pixels is ignored when \"merging\" only one image" << std::endl;
					if ( ! std::isnan( stops ) )
						std::cerr << "Warning: stops is ignored when \"merging\" only one image" << std::endl;
					if ( ! std::isnan( lambda ) )
						std::cerr << "Warning: lambda is ignored when \"merging\" only one image" << std::endl;
					if ( ! std::isnan( epsilon ) )
						std::cerr << "Warning: epsilon is ignored when \"merging\" only one image" << std::endl;

					std::cout << "Loading image" << std::endl;
					SharedArray< uint16_t > inputImage = LoadLinearTIFF( inputs[ 0 ].c_str(), rows, columns );
					image = SharedArray< float >( rows * columns * 3 );
					#pragma omp parallel for schedule( static )
					for ( int ii = 0; ii < static_cast< int >( rows * columns * 3 ); ++ii )
						image[ ii ] = static_cast< double >( inputImage[ ii ] ) / 0xffff;
				}
			}
			else {

				if ( inputs.size() > 1 ) {

					if ( std::isnan( pixels ) )
						pixels = 10000;
					if ( std::isnan( stops ) )
						throw std::runtime_error( "stops must be specified whem merging (except for 16-bit linear TIFFs)" );
					if ( std::isnan( lambda ) )
						lambda = 1;
					if ( std::isnan( epsilon ) )
						epsilon = 0.0001;

					std::vector< SharedArray< uint8_t > > inputImages;
					inputImages.push_back( LoadNonlinearImage( inputs[ 0 ].c_str(), rows, columns ) );
					for ( unsigned int ii = 1; ii < inputs.size(); ++ii ) {

						unsigned int newRows, newColumns;
						inputImages.push_back( LoadNonlinearImage( inputs[ ii ].c_str(), newRows, newColumns ) );
						if ( ( newRows != rows ) || ( newColumns != columns ) )
							throw std::runtime_error( "all input images must have the same dimensions" );
					}

					image = MergeNonlinearImages( inputImages, rows, columns, pixels, stops, lambda, epsilon );
				}
				else {

					if ( ! std::isnan( pixels ) )
						std::cerr << "Warning: pixels is ignored when \"merging\" only one image" << std::endl;
					if ( ! std::isnan( stops ) )
						std::cerr << "Warning: stops is ignored when \"merging\" only one image" << std::endl;
					if ( ! std::isnan( lambda ) )
						std::cerr << "Warning: lambda is ignored when \"merging\" only one image" << std::endl;
					if ( ! std::isnan( epsilon ) )
						std::cerr << "Warning: epsilon is ignored when \"merging\" only one image" << std::endl;

					std::cout << "Loading image" << std::endl;
					SharedArray< uint8_t > inputImage = LoadNonlinearImage( inputs[ 0 ].c_str(), rows, columns );
					image = SharedArray< float >( rows * columns * 3 );
					#pragma omp parallel for schedule( static )
					for ( int ii = 0; ii < static_cast< int >( rows * columns * 3 ); ++ii ) {

						// attempt to undo gamma correction
						double const gamma = static_cast< double >( inputImage[ ii ] ) / 0xff;
						if ( gamma < 4.5 * 0.018 )
							image[ ii ] = gamma / 4.5;
						else
							image[ ii ] = std::pow( ( gamma + 0.099 ) / 1.099, 1.0 / 0.45 );
					}
				}
			}

			SaveEXR( output, image, rows, columns );
		}
	}
	catch( std::exception& exception ) {

		std::cerr << "Error: " << exception.what() << std::endl << std::endl;
		poptPrintHelp( context, stderr, 0 );
		code = EXIT_FAILURE;
	}

	poptFreeContext( context );

	return code;
}
