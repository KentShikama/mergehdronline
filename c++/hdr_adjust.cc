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
	\file hdr_adjust.cc
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
	double redScale   = std::numeric_limits< double >::quiet_NaN();
	double redGamma   = std::numeric_limits< double >::quiet_NaN();
	double greenScale = std::numeric_limits< double >::quiet_NaN();
	double greenGamma = std::numeric_limits< double >::quiet_NaN();
	double blueScale  = std::numeric_limits< double >::quiet_NaN();
	double blueGamma  = std::numeric_limits< double >::quiet_NaN();
	char* crop = NULL;

	struct poptOption options[] = {
		{ "help",        'h',  POPT_ARG_NONE,   &help,       0, "show help",                         NULL        },
		{ "input",       'i',  POPT_ARG_STRING, &input,      0, "input file (EXR)",                  "INPUT"     },
		{ "output",      'o',  POPT_ARG_STRING, &output,     0, "output file (EXR)",                 "OUTPUT"    },
		{ "red_scale",   'r',  POPT_ARG_DOUBLE, &redScale,   0, "scale parameter for red channel",   "FLOAT(=1)" },
		{ "red_gamma",   'R',  POPT_ARG_DOUBLE, &redGamma,   0, "gamma parameter for red channel",   "FLOAT(=1)" },
		{ "green_scale", 'g',  POPT_ARG_DOUBLE, &greenScale, 0, "scale parameter for green channel", "FLOAT(=1)" },
		{ "green_gamma", 'G',  POPT_ARG_DOUBLE, &greenGamma, 0, "gamma parameter for green channel", "FLOAT(=1)" },
		{ "blue_scale",  'b',  POPT_ARG_DOUBLE, &blueScale,  0, "scale parameter for blue channel",  "FLOAT(=1)" },
		{ "blue_gamma",  'B',  POPT_ARG_DOUBLE, &blueGamma,  0, "gamma parameter for blue channel",  "FLOAT(=1)" },
		{ "crop",        'c',  POPT_ARG_STRING, &crop,       0, "cropping",                          "WxH+X+Y"   },
		{ NULL,          '\0', 0,               NULL,        0, NULL,                                NULL        }
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
				"Performs cropping and/or simple color adjustments on an HDR input file, in EXR" << std::endl <<
				"format, and saves the result to a new HDR file, also in EXR format." << std::endl <<
				std::endl <<
				"For each of the three color channels (red, green and blue) there are scale and" << std::endl <<
				"gamma parameters. At each pixel, the value of each color channel is replaced by" << std::endl <<
				"scale*(value^gamma). One can (very roughly) think of the scale parameter as" << std::endl <<
				"controlling the strength of a color channel, and the gamma parameter its" << std::endl <<
				"saturation, although adjusting the gamma parmeters of the different color" << std::endl <<
				"channels independently often leads to unpredictable results." << std::endl <<
				std::endl
			;
			poptPrintHelp( context, stdout, 0 );
		}
		else {

			if ( input == NULL )
				throw std::runtime_error( "input file not specified" );
			if ( output == NULL )
				throw std::runtime_error( "output file not specified" );
			if ( std::isnan( redScale ) )
				redScale = 1;
			if ( std::isnan( redGamma ) )
				redGamma = 1;
			if ( std::isnan( greenScale ) )
				greenScale = 1;
			if ( std::isnan( greenGamma ) )
				greenGamma = 1;
			if ( std::isnan( blueScale ) )
				blueScale = 1;
			if ( std::isnan( blueGamma ) )
				blueGamma = 1;

			std::cout << "Loading HDR image" << std::endl;
			unsigned int rows, columns;
			SharedArray< float > image = LoadEXR( input, rows, columns );

			std::cout << "Adjusting image" << std::endl;
			if ( crop != NULL ) {

				unsigned int cropRows    = 0;
				unsigned int cropColumns = 0;
				unsigned int cropRow     = 0;
				unsigned int cropColumn  = 0;
				{	char const* ii = crop;

					cropColumns = 0;
					for ( ; *ii != 'x'; ++ii ) {

						cropColumns *= 10;
						if ( ( *ii < '0' ) || ( *ii > '9' ) )
							throw std::runtime_error( "cropping window must be specified as WxH+X+Y" );
						cropColumns += *ii - '0';
					}
					++ii;

					cropRows = 0;
					for ( ; *ii != '+'; ++ii ) {

						cropRows *= 10;
						if ( ( *ii < '0' ) || ( *ii > '9' ) )
							throw std::runtime_error( "cropping window must be specified as WxH+X+Y" );
						cropRows += *ii - '0';
					}
					++ii;

					cropColumn = 0;
					for ( ; *ii != '+'; ++ii ) {

						cropColumn *= 10;
						if ( ( *ii < '0' ) || ( *ii > '9' ) )
							throw std::runtime_error( "cropping window must be specified as WxH+X+Y" );
						cropColumn += *ii - '0';
					}
					++ii;

					cropRow = 0;
					for ( ; *ii != '\0'; ++ii ) {

						cropRow *= 10;
						if ( ( *ii < '0' ) || ( *ii > '9' ) )
							throw std::runtime_error( "cropping window must be specified as WxH+X+Y" );
						cropRow += *ii - '0';
					}
					++ii;
				}

				if ( ( cropColumn == 0 ) || ( cropRow == 0 ) )
					throw std::runtime_error( "cropping window must have nonzero dimensions" );
				if ( ( cropColumns + cropColumn > columns ) || ( cropRows + cropRow > rows ) )
					throw std::runtime_error( "cropping window must fit within the image" );

				SharedArray< float > oldImage = image;
				image = SharedArray< float >( cropRows * cropColumns * 3 );

				#pragma omp parallel for schedule( static )
				for ( int row = 0; row < static_cast< int >( cropRows ); ++row ) {

					float const* const pSource = &oldImage[ ( ( row + cropRow ) * columns + cropColumn ) * 3 ];
					float* const pDestination = &image[ row * cropColumns * 3 ];
					std::copy( pSource, pSource + cropColumns * 3, pDestination );
				}

				rows    = cropRows;
				columns = cropColumns;
			}
			{	if ( redScale <= 0 )
					throw std::runtime_error( "red_scale must be positive" );
				if ( greenScale <= 0 )
					throw std::runtime_error( "green_scale must be positive" );
				if ( blueScale <= 0 )
					throw std::runtime_error( "blue_scale must be positive" );

				float oldMaximum = -std::numeric_limits< float >::infinity();
				float newMaximum = -std::numeric_limits< float >::infinity();

				#pragma omp parallel
				{	float localOldMaximum = -std::numeric_limits< float >::infinity();
					float localNewMaximum = -std::numeric_limits< float >::infinity();
					#pragma omp for schedule( static )
					for ( int ii = 0; ii < static_cast< int >( rows * columns ); ++ii ) {

						localOldMaximum = std::max( localOldMaximum, image[ ii * 3 + 0 ] );
						localOldMaximum = std::max( localOldMaximum, image[ ii * 3 + 1 ] );
						localOldMaximum = std::max( localOldMaximum, image[ ii * 3 + 2 ] );
						image[ ii * 3 + 0 ] = redScale   * std::pow( static_cast< double >( image[ ii * 3 + 0 ] ), redGamma   );
						image[ ii * 3 + 1 ] = greenScale * std::pow( static_cast< double >( image[ ii * 3 + 1 ] ), greenGamma );
						image[ ii * 3 + 2 ] = blueScale  * std::pow( static_cast< double >( image[ ii * 3 + 2 ] ), blueGamma  );
						localNewMaximum = std::max( localNewMaximum, image[ ii * 3 + 0 ] );
						localNewMaximum = std::max( localNewMaximum, image[ ii * 3 + 1 ] );
						localNewMaximum = std::max( localNewMaximum, image[ ii * 3 + 2 ] );
					}

					#pragma omp critical
					{	oldMaximum = std::max( oldMaximum, localOldMaximum );
						newMaximum = std::max( newMaximum, localNewMaximum );
					}
				}

				if ( std::fabs( newMaximum ) <= std::numeric_limits< float >::epsilon() )
					throw std::runtime_error( "image was scaled to black" );

				double const scale = static_cast< double >( oldMaximum ) / static_cast< double >( newMaximum );
				#pragma omp parallel for schedule( static )
				for ( int ii = 0; ii < static_cast< int >( rows * columns * 3 ); ++ii )
					image[ ii ] *= scale;
			}

			std::cout << "Saving HDR image to " << output << std::endl;
			SaveEXR( output, image, rows, columns );
		}
	}
	catch( std::exception& exception ) {

		std::cerr << "Error: " << exception.what() << std::endl << std::endl;
		poptPrintHelp( context, stderr, 0 );
		code = EXIT_FAILURE;
	}

	return code;
}
