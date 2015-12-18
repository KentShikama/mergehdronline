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
	\file image.cc
	\brief implementation of image manipulation functions
*/




#include "headers.hh"




namespace {




//============================================================================
//    MergeWeight function
//============================================================================


template< typename t_Type >
inline double MergeWeightFunction( t_Type const& value, bool const low = false, bool const high = false ) {

	double const maximum = static_cast< t_Type >( -1 );
	double const lower = maximum * 0.2;
	double const upper = maximum * 0.8;

	double weight = 1;
	if ( value < lower ) {

		if ( ! low )
			weight = value / lower;
	}
	else if ( value > upper ) {

		if ( ! high )
			weight = ( maximum - value ) / ( maximum - upper );
	}

	return weight;
}




}    // anonymous namespace




//============================================================================
//    SaveEXR function
//============================================================================


void SaveEXR(
	char const* const filename,
	SharedArray< float > const& image,
	unsigned int const rows,
	unsigned int const columns
)
{
	SharedArray< Imf::Rgba > pixels( rows * columns );

	#pragma omp parallel for schedule( static )
	for ( int ii = 0; ii < static_cast< int >( rows * columns ); ++ii ) {

		pixels[ ii ].r = half( image[ ii * 3 + 0 ] );
		pixels[ ii ].g = half( image[ ii * 3 + 1 ] );
		pixels[ ii ].b = half( image[ ii * 3 + 2 ] );
		pixels[ ii ].a = half( 1 );
	}

	Imf::RgbaOutputFile file( filename, Imf::Header( columns, rows ) );
	file.setFrameBuffer( &pixels[ 0 ], 1, columns );
	file.writePixels( rows );
}




//============================================================================
//    SaveNonlinearImage function
//============================================================================


void SaveNonlinearImage(
	char const* const filename,
	SharedArray< float > const& image,
	unsigned int const rows,
	unsigned int const columns,
	double const stops
)
{
	double const maximum = std::pow( 2.0, -stops );

	Magick::Image file( Magick::Geometry( columns, rows ), "black" );
	file.type( Magick::TrueColorType );
	assert( file.colorSpace() == Magick::RGBColorspace );

#ifdef MAGICKCORE_HDRI_SUPPORT
	float const scale = static_cast< float >( QuantumRange );
#else    // MAGICKCORE_HDRI_SUPPORT
	float const scale = ( 1u << QuantumDepth ) - 1;
#endif    // MAGICKCORE_HDRI_SUPPORT

	Magick::PixelPacket* pPixelPacket = file.getPixels( 0, 0, columns, rows );

	#pragma omp parallel for schedule( static )
	for ( int row = 0; row < static_cast< int >( rows ); ++row ) {

		Magick::PixelPacket* ii = pPixelPacket + row * columns;
		float const* jj = &image[ 0 ] + row * columns * 3;
		for ( unsigned int column = 0; column < columns; ++column ) {

			double red   = *( jj + 0 ) / maximum;
			double green = *( jj + 1 ) / maximum;
			double blue  = *( jj + 2 ) / maximum;

			// gamma correction
			red   = ( red   < 0.018 ) ? ( 4.5 * red   ) : ( 1.099 * std::pow( red,   0.45 ) - 0.099 );
			green = ( green < 0.018 ) ? ( 4.5 * green ) : ( 1.099 * std::pow( green, 0.45 ) - 0.099 );
			blue  = ( blue  < 0.018 ) ? ( 4.5 * blue  ) : ( 1.099 * std::pow( blue,  0.45 ) - 0.099 );

			ii->red   = static_cast< MagickCore::Quantum >( std::max( 0.0, std::min( 1.0, red   ) ) * scale );
			ii->green = static_cast< MagickCore::Quantum >( std::max( 0.0, std::min( 1.0, green ) ) * scale );
			ii->blue  = static_cast< MagickCore::Quantum >( std::max( 0.0, std::min( 1.0, blue  ) ) * scale );
			ii->opacity = 0;

			++ii;
			jj += 3;
		}
	}

	file.syncPixels();
	file.write( filename );
}




//============================================================================
//    LoadEXR function
//============================================================================


SharedArray< float > const LoadEXR(
	char const* const filename,
	unsigned int& rows,
	unsigned int& columns
)
{
	Imf::RgbaInputFile file( filename );

	Imath::Box2i box = file.dataWindow();
	columns = box.max.x + 1 - box.min.x;
	rows    = box.max.y + 1 - box.min.y;

	SharedArray< Imf::Rgba > pixels( rows * columns );
	file.setFrameBuffer( &pixels[ 0 ] - ( box.min.y * columns + box.min.x ), 1, columns );
	file.readPixels( box.min.y, box.max.y );

	SharedArray< float > result( rows * columns * 3 );

	#pragma omp parallel for schedule( static )
	for ( int ii = 0; ii < static_cast< int >( rows * columns ); ++ii ) {

		result[ ii * 3 + 0 ] = static_cast< float >( pixels[ ii ].r );
		result[ ii * 3 + 1 ] = static_cast< float >( pixels[ ii ].g );
		result[ ii * 3 + 2 ] = static_cast< float >( pixels[ ii ].b );
	}

	return result;
}




//============================================================================
//    LoadLinearTIFF function
//============================================================================


SharedArray< uint16_t > const LoadLinearTIFF(
	char const* const filename,
	unsigned int& rows,
	unsigned int& columns
)
{
	TIFF* file = TIFFOpen( filename, "r" );
	if ( file == NULL )
		throw std::runtime_error( "unable to open image file" );

	if ( ! TIFFGetField( file, TIFFTAG_IMAGELENGTH, &rows ) )
		throw std::runtime_error( "unable to read rows" );
	if ( ! TIFFGetField( file, TIFFTAG_IMAGEWIDTH, &columns ) )
		throw std::runtime_error( "unable to read columns" );

	unsigned int bits = 0;;
	if ( ! TIFFGetField( file, TIFFTAG_BITSPERSAMPLE, &bits ) )
		throw std::runtime_error( "unable to read bits" );
	if ( bits != 16 )
		throw std::runtime_error( "expected 16-bit linear color channels" );

	unsigned int const bytes = TIFFScanlineSize( file );
	if ( bytes != columns * 6 )
		throw std::runtime_error( "wrong number of bytes per scanline" );

	SharedArray< uint16_t > result( rows * columns * 3 );

	tdata_t buffer = _TIFFmalloc( bytes );
	for ( unsigned int row = 0; row < rows; ++row ) {

		TIFFReadScanline( file, buffer, row );

		uint16_t const* ii = static_cast< uint16_t const* >( buffer );
		uint16_t* jj = &result[ 0 ] + row * columns * 3;
		for ( unsigned int column = 0; column < columns * 3; ++column )
			*( jj++ ) = *( ii++ );
	}
	_TIFFfree( buffer );

	TIFFClose( file );

	return result;
}




//============================================================================
//    LoadNonlinearImage function
//============================================================================


SharedArray< uint8_t > const LoadNonlinearImage(
	char const* const filename,
	unsigned int& rows,
	unsigned int& columns
)
{
	Magick::Image file( filename );
	rows    = file.rows();
	columns = file.columns();

	SharedArray< uint8_t > result( rows * columns * 3 );

#ifdef MAGICKCORE_HDRI_SUPPORT
	double const scale = 256 / ( static_cast< double >( QuantumRange ) + 1 );
#else    // MAGICKCORE_HDRI_SUPPORT
	double const scale = 256 / static_cast< double >( 1u << QuantumDepth );
#endif    // MAGICKCORE_HDRI_SUPPORT

	Magick::PixelPacket const* pPixelPacket( file.getConstPixels( 0, 0, columns, rows ) );

	#pragma omp parallel for schedule( static )
	for ( int row = 0; row < static_cast< int >( rows ); ++row ) {

		Magick::PixelPacket const* ii = pPixelPacket + row * columns;
		uint8_t* jj = &result[ 0 ] + row * columns * 3;
		for ( unsigned int column = 0; column < columns; ++column ) {

			*( jj + 0 ) = std::max( 0, std::min( 255, static_cast< int >( ii->red   * scale ) ) );
			*( jj + 1 ) = std::max( 0, std::min( 255, static_cast< int >( ii->green * scale ) ) );
			*( jj + 2 ) = std::max( 0, std::min( 255, static_cast< int >( ii->blue  * scale ) ) );

			++ii;
			jj += 3;
		}
	}

	return result;
}




//============================================================================
//    MergeLinearImages function
//============================================================================


/*
	This is basically the algorithm of Debevec & Malik, where a linear response
	curve is substituted into the equations.
*/
SharedArray< float > MergeLinearImages(
	std::vector< SharedArray< uint16_t > > const& images,
	unsigned int const rows,
	unsigned int const columns,
	double const pixels
)
{
	if ( pixels <= 0 )
		throw std::runtime_error( "pixels must be positive" );

	unsigned int const size = images.size();

	SharedArray< double > matrix( Square( size ) );
	SharedArray< double > vector(         size   );
	std::fill( &matrix[ 0 ], &matrix[ 0 ] + Square( size ), 0 );
	std::fill( &vector[ 0 ], &vector[ 0 ] +         size,   0 );

	// calculate the optimization problem
	{	unsigned int const sampleRows    = std::max( 1u, std::min( rows,    static_cast< unsigned int >( std::ceil( std::sqrt( pixels * rows    / columns ) ) ) ) );
		unsigned int const sampleColumns = std::max( 1u, std::min( columns, static_cast< unsigned int >( std::ceil( std::sqrt( pixels * columns / rows    ) ) ) ) );

		SharedArray< double > coefficients( images.size() );
		SharedArray< double > logValues(    images.size() );
		SharedArray< double > weights(      images.size() );

		for ( unsigned int sampleRow = 0; sampleRow < sampleRows; ++sampleRow ) {

			unsigned int const row = ( rows * ( 2 * sampleRow + 1 ) ) / ( 2 * sampleRows );

			for ( unsigned int sampleColumn = 0; sampleColumn < sampleColumns; ++sampleColumn ) {

				unsigned int const column = ( columns * ( 2 * sampleColumn + 1 ) ) / ( 2 * sampleColumns );

				unsigned int index = ( row * columns + column ) * 3;
				for ( unsigned int channel = 0; channel < 3; ++channel ) {

					double denominator = 0;
					for ( unsigned int ii = 0; ii < images.size(); ++ii ) {

						uint16_t const value = images[ ii ][ index ];
						logValues[ ii ] = ( ( value == 0 ) ? 0 : std::log( value ) );
						weights[   ii ] = MergeWeightFunction( value );

						denominator += weights[ ii ];
					}

					if ( denominator > 0 ) {

						for ( unsigned int ii = 0; ii < images.size(); ++ii ) {

							if ( weights[ ii ] > 0 ) {

								double constant = 0;
								std::fill( &coefficients[ 0 ], &coefficients[ 0 ] + images.size(), 0 );

								constant += logValues[ ii ];
								coefficients[ ii ] -= 1;
								for ( unsigned int jj = 0; jj < images.size(); ++jj ) {

									if ( weights[ jj ] > 0 ) {

										double const weight = weights[ jj ] / denominator;
										constant -= weight * logValues[ jj ];
										coefficients[ jj ] += weight;
									}
								}

								unsigned int matrixIndex = 0;
								for ( unsigned int jj = 0; jj < images.size(); ++jj ) {

									vector[ jj ] -= weights[ ii ] * constant * coefficients[ jj ];

									for ( unsigned int kk = 0; kk < images.size(); ++kk ) {

										matrix[ matrixIndex ] += weights[ ii ] * coefficients[ jj ] * coefficients[ kk ];
										++matrixIndex;
									}
								}
							}
						}
					}

					++index;
				}
			}
		}
	}

	// conjugate gradient
	SharedArray< double > solution( size );
	std::fill( &solution[ 0 ], &solution[ 0 ] + size, 0 );
	{	SharedArray< double > gradient( size );
		SharedArray< double > step(     size );

		std::fill( &gradient[ 0 ], &gradient[ 0 ] + size, 0 );
		std::fill( &step[     0 ], &step[     0 ] + size, 0 );

		double gradientNormSquared = 0;
		for ( unsigned int ii = 0; ii < size; ++ii ) {

			double const oldGradientNormSquared = gradientNormSquared;

			// calculate the gradient
			std::copy( &vector[ 0 ], &vector[ 0 ] + size, &gradient[ 0 ] );
			for ( unsigned int jj = 0; jj < size; ++jj ) {

				double accumulator = 0;
				double const* const row = &matrix[ 0 ] + jj * size;
				for ( unsigned int kk = 0; kk < size; ++kk )
					accumulator += row[ kk ] * solution[ kk ];
				gradient[ jj ] -= accumulator;
			}

			// project the gradient onto the constraints
			{	double average = 0;
				for ( unsigned int jj = 0; jj < size; ++jj )
					average += gradient[ jj ];
				average /= size;
				for ( unsigned int jj = 0; jj < size; ++jj )
					gradient[ jj ] -= average;
			}

			// find the gradient's squared norm
			gradientNormSquared = 0;
			for ( unsigned int jj = 0; jj < size; ++jj )
				gradientNormSquared += Square( gradient[ jj ] );

			// update the step
			double beta = 0;
			if ( std::fabs( oldGradientNormSquared ) > std::numeric_limits< double >::epsilon() )
				beta = gradientNormSquared / oldGradientNormSquared;
			for ( unsigned int jj = 0; jj < size; ++jj )
				step[ jj ] = beta * step[ jj ] + gradient[ jj ];

			// find the step size
			double alphaDenominator = 0;
			for ( unsigned int jj = 0; jj < size; ++jj ) {

				double accumulator = 0;
				double const* const row = &matrix[ 0 ] + jj * size;
				for ( unsigned int kk = 0; kk < size; ++kk )
					accumulator += row[ kk ] * step[ kk ];
				alphaDenominator += step[ jj ] * accumulator;
			}
			if ( std::fabs( alphaDenominator ) < std::numeric_limits< double >::epsilon() )
				break;
			double const alpha = gradientNormSquared / alphaDenominator;

			// take a step
			for ( unsigned int jj = 0; jj < size; ++jj )
				solution[ jj ] += alpha * step[ jj ];
		}
	}

	unsigned int minimumIndex = 0;
	unsigned int maximumIndex = 0;
	{	double minimumSolution =  std::numeric_limits< double >::infinity();
		double maximumSolution = -std::numeric_limits< double >::infinity();
		for ( unsigned int ii = 0; ii < size; ++ii ) {

			if ( solution[ ii ] < minimumSolution ) {

				minimumIndex = ii;
				minimumSolution = solution[ ii ];
			}
			if ( solution[ ii ] > maximumSolution ) {

				maximumIndex = ii;
				maximumSolution = solution[ ii ];
			}
		}
	}

	// create the output image
	SharedArray< float > result( rows * columns * 3 );
	float maximum = -std::numeric_limits< float >::infinity();
	#pragma omp parallel
	{	float localMaximum = -std::numeric_limits< float >::infinity();

		#pragma omp for schedule( static )
		for ( int ii = 0; ii < static_cast< int >( rows * columns ); ++ii ) {

			double red   = 0;
			double green = 0;
			double blue  = 0;
			double denominator = 0;

			for ( unsigned int jj = 0; jj < images.size(); ++jj ) {

				uint16_t const redValue   = images[ jj ][ ii * 3 + 0 ];
				uint16_t const greenValue = images[ jj ][ ii * 3 + 1 ];
				uint16_t const blueValue  = images[ jj ][ ii * 3 + 2 ];

				double redWeight   = 0;
				double greenWeight = 0;
				double blueWeight  = 0;
				{	bool low  = false;
					bool high = false;
					if ( jj == maximumIndex )    // lightest image has full information on shadows
						low = true;
					if ( jj == minimumIndex )    // darkest image has full information on highlights
						high = true;
					redWeight   = MergeWeightFunction( redValue,   low, high );
					greenWeight = MergeWeightFunction( greenValue, low, high );
					blueWeight  = MergeWeightFunction( blueValue,  low, high );
				}
				double const weight = redWeight * greenWeight * blueWeight + std::numeric_limits< float >::epsilon();

				red   += ( ( redValue   == 0 ) ? -solution[ jj ] : ( std::log( redValue   ) - solution[ jj ] ) ) * weight;
				green += ( ( greenValue == 0 ) ? -solution[ jj ] : ( std::log( greenValue ) - solution[ jj ] ) ) * weight;
				blue  += ( ( blueValue  == 0 ) ? -solution[ jj ] : ( std::log( blueValue  ) - solution[ jj ] ) ) * weight;
				denominator += weight;
			}

			result[ ii * 3 + 0 ] = red   / denominator;
			result[ ii * 3 + 1 ] = green / denominator;
			result[ ii * 3 + 2 ] = blue  / denominator;

			localMaximum = std::max( localMaximum, static_cast< float >( result[ ii * 3 + 0 ] ) );
			localMaximum = std::max( localMaximum, static_cast< float >( result[ ii * 3 + 1 ] ) );
			localMaximum = std::max( localMaximum, static_cast< float >( result[ ii * 3 + 2 ] ) );
		}

		#pragma omp critical
		maximum = std::max( maximum, localMaximum );
	}

	// scale the output to [0,1]
	#pragma omp parallel for schedule( static )
	for ( int ii = 0; ii < static_cast< int >( rows * columns * 3 ); ++ii )
		result[ ii ] = std::exp( result[ ii ] - maximum );

	return result;
}




//============================================================================
//    MergeNonlinearImages function
//============================================================================


/*
	Debevec, Malik. "Recovering high dynamic range radiance maps from
	photographs".

	The optimization problem is solved subject to monotonicity constraints, to
	which Debevec and Maloc allude, without making it clear whether they impose
	them, or merely hope that they will be satisfied at the solution. Actually
	imposing them as constraints significantly improves the quality of the
	solution, at the cost of making the problem much more difficult to
	optimize. We use the projected gradient algorithm with conjugate gradient
	inner iterations.
*/
SharedArray< float > MergeNonlinearImages(
	std::vector< SharedArray< uint8_t > > const& images,
	unsigned int const rows,
	unsigned int const columns,
	double const pixels,
	double const stops,
	double const lambda,
	double const epsilon
)
{
	if ( pixels <= 0 )
		throw std::runtime_error( "pixels must be positive" );
	if ( stops <= 0 )
		throw std::runtime_error( "stops must be positive" );
	if ( lambda <= 0 )
		throw std::runtime_error( "lambda must be positive" );
	if ( epsilon <= 0 )
		throw std::runtime_error( "epsilon must be positive" );

	unsigned int const size = 256 * 3 + images.size();

	SharedArray< double > matrix( Square( size ) );
	std::fill( &matrix[ 0 ], &matrix[ 0 ] + Square( size ), 0 );

	// calculate the data term of the matrix
	unsigned int truePixels = 0;
	{	unsigned int const sampleRows    = std::max( 1u, std::min( rows,    static_cast< unsigned int >( std::ceil( std::sqrt( pixels * rows    / columns ) ) ) ) );
		unsigned int const sampleColumns = std::max( 1u, std::min( columns, static_cast< unsigned int >( std::ceil( std::sqrt( pixels * columns / rows    ) ) ) ) );

		SharedArray< unsigned int > indices(      images.size() * 2 );
		SharedArray< double       > coefficients( images.size() * 2 );
		SharedArray< double       > weights(      images.size()     );

		for ( unsigned int sampleRow = 0; sampleRow < sampleRows; ++sampleRow ) {

			unsigned int const row = ( rows * ( 2 * sampleRow + 1 ) ) / ( 2 * sampleRows );

			for ( unsigned int sampleColumn = 0; sampleColumn < sampleColumns; ++sampleColumn ) {

				unsigned int const column = ( columns * ( 2 * sampleColumn + 1 ) ) / ( 2 * sampleColumns );

				unsigned int index = ( row * columns + column ) * 3;
				for ( unsigned int channel = 0; channel < 3; ++channel ) {

					double denominator = 0;
					for ( unsigned int ii = 0; ii < images.size(); ++ii ) {

						uint8_t const value = images[ ii ][ index ];
						weights[ ii ] = MergeWeightFunction( value );

						denominator += weights[ ii ];
					}

					if ( denominator > 0 ) {

						for ( unsigned int ii = 0; ii < images.size(); ++ii ) {

							indices[ ii * 2 + 0 ] = 256 * channel + images[ ii ][ index ];
							indices[ ii * 2 + 1 ] = 256 * 3 + ii;
						}

						for ( unsigned int ii = 0; ii < images.size(); ++ii ) {

							if ( weights[ ii ] > 0 ) {

								std::fill( &coefficients[ 0 ], &coefficients[ 0 ] + images.size() * 2, 0 );

								coefficients[ ii * 2 + 0 ] += 1;
								coefficients[ ii * 2 + 1 ] -= 1;
								for ( unsigned int jj = 0; jj < images.size(); ++jj ) {

									if ( weights[ jj ] > 0 ) {

										double const weight = weights[ jj ] / denominator;
										coefficients[ jj * 2 + 0 ] -= weight;
										coefficients[ jj * 2 + 1 ] += weight;
									}
								}

								for ( unsigned int jj = 0; jj < images.size() * 2; ++jj )
									for ( unsigned int kk = 0; kk < images.size() * 2; ++kk )
										matrix[ indices[ jj ] * size + indices[ kk ] ] += weights[ ii ] * coefficients[ jj ] * coefficients[ kk ];
							}
						}
					}

					++index;
				}
			}
		}

		truePixels = sampleRows * sampleColumns;
	}

	// normalize the data term of the matrix
	{	double const scale = 254.0 / ( images.size() * truePixels );
		#pragma omp parallel for schedule( static )
		for ( int ii = 0; ii < static_cast< int >( Square( size ) ); ++ii )
			matrix[ ii ] *= scale;
	}

	// add in the smoothness term
	for ( unsigned int channel = 0; channel < 3; ++channel ) {

		for ( unsigned int ii = 1; ii < 255; ++ii ) {

			unsigned int const offset = ( 256 * channel + ii ) * ( size + 1 );

			matrix[ offset - size - 1 ] += lambda;
			matrix[ offset - size     ] -= 2 * lambda;
			matrix[ offset - size + 1 ] += lambda;
			matrix[ offset        - 1 ] -= 2 * lambda;
			matrix[ offset            ] += 4 * lambda;
			matrix[ offset        + 1 ] -= 2 * lambda;
			matrix[ offset + size - 1 ] += lambda;
			matrix[ offset + size     ] -= 2 * lambda;
			matrix[ offset + size + 1 ] += lambda;
		}
	}

	unsigned int minimumIndex = 0;
	unsigned int maximumIndex = 0;
	{	double minimumIntensity =  std::numeric_limits< double >::infinity();
		double maximumIntensity = -std::numeric_limits< double >::infinity();
		for ( unsigned int ii = 0; ii < images.size(); ++ii ) {

			double intensity = 0;
			#pragma omp parallel for schedule( static ) reduction( + : intensity )
			for ( int jj = 0; jj < static_cast< int >( rows * columns * 3 ); ++jj )
				intensity += images[ ii ][ jj ];

			if ( intensity < minimumIntensity ) {

				minimumIndex     = ii;
				minimumIntensity = intensity;
			}
			else if ( intensity > maximumIntensity ) {

				maximumIndex     = ii;
				maximumIntensity = intensity;
			}
		}
	}
	if ( minimumIndex == maximumIndex )
		throw std::runtime_error( "unable to identify lightest and darkest images" );

	// projected gradient with conjugate gradient inner iterations
	SharedArray< double > solution( size );
	std::fill( &solution[ 0 ], &solution[ 0 ] + size, 0 );
	solution[ 256 * 3 + maximumIndex ] = stops * std::log( 2 );
	{	SharedArray< bool > constraints( 256 * 3 );
		SharedArray< double > gradient( size );
		SharedArray< double > step(     size );

		std::fill( &gradient[ 0 ], &gradient[ 0 ] + size, 0 );
		std::fill( &step[     0 ], &step[     0 ] + size, 0 );

		for ( bool progress = true; progress; ) {

			std::fill( &constraints[ 0 ], &constraints[ 0 ] + 256 * 3, false );

			// calculate the gradient
			#pragma omp parallel for schedule( static )
			for ( int ii = 0; ii < static_cast< int >( size ); ++ii ) {

				double accumulator = 0;
				double const* const row = &matrix[ 0 ] + ii * size;
				for ( unsigned int jj = 0; jj < size; ++jj )
					accumulator += row[ jj ] * solution[ jj ];
				gradient[ ii ] = -accumulator;
			}

			// project the gradient onto the equality constraints
			for ( unsigned int channel = 0; channel < 3; ++channel ) {

				double average = 0;
				for ( unsigned int ii = 256 * channel; ii < 256 * ( channel + 1 ); ++ii )
					average += gradient[ ii ];
				average /= 256;
				for ( unsigned int ii = 256 * channel; ii < 256 * ( channel + 1 ); ++ii )
					gradient[ ii ] -= average;
			}
			gradient[ 256 * 3 + minimumIndex ] = 0;
			gradient[ 256 * 3 + maximumIndex ] = 0;
			// project the gradient onto the inequality constraints
			for ( unsigned int channel = 0; channel < 3; ++channel ) {

				for ( bool done = false; ! done; ) {

					done = true;

					double numerator = 0;
					unsigned int denominator = 0;
					for ( unsigned int ii = 256 * channel; ii < 256 * ( channel + 1 ) ; ++ii ) {

						numerator += gradient[ ii ];
						++denominator;

						double const average = numerator / denominator;
						if ( ( ( ii & 255 ) < 255 ) && ( constraints[ ii ] || ( ( solution[ ii + 1 ] <= solution[ ii ] ) && ( gradient[ ii + 1 ] <= average ) ) ) ) {

							if ( ! constraints[ ii ] ) {

								done = false;
								constraints[ ii ] = true;
							}
						}
						else {

							if ( denominator > 1 ) {

								for ( unsigned int jj = 0; jj < denominator; ++jj )
									gradient[ ii - jj ] = average;
							}

							numerator   = 0;
							denominator = 0;
						}
					}
				}
			}

			// find the gradient's squared norm
			double gradientNormSquared = 0;
			for ( unsigned int ii = 0; ii < size; ++ii )
				gradientNormSquared += Square( gradient[ ii ] );
			if ( std::sqrt( gradientNormSquared ) < epsilon )
				break;

			std::copy( &gradient[ 0 ], &gradient[ 0 ] + size, &step[ 0 ] );

			for ( ; ; ) {

				double const oldGradientNormSquared = gradientNormSquared;

				// find the step size
				double alphaDenominator = 0;
				#pragma omp parallel for schedule( static ) reduction( + : alphaDenominator )
				for ( int ii = 0; ii < static_cast< int >( size ); ++ii ) {

					double accumulator = 0;
					double const* const row = &matrix[ 0 ] + ii * size;
					for ( unsigned int jj = 0; jj < size; ++jj )
						accumulator += row[ jj ] * step[ jj ];
					alphaDenominator += step[ ii ] * accumulator;
				}
				if ( std::fabs( alphaDenominator ) < std::numeric_limits< double >::epsilon() )
					break;
				double alpha = gradientNormSquared / alphaDenominator;

				// check the constraints
				int constraintIndex = -1;
				for ( unsigned int channel = 0; channel < 3; ++channel ) {

					for ( unsigned int ii = 256 * channel; ii < 256 * channel + 255; ++ii ) {

						if ( ! constraints[ ii ] ) {

							double const deltaStep     = step[     ii + 1 ] - step[     ii ];
							double const deltaSolution = solution[ ii + 1 ] - solution[ ii ];
							if ( deltaStep < 0 ) {

								double threshold = -deltaSolution / deltaStep;
								if ( alpha > threshold ) {

									constraintIndex = ii;
									alpha = threshold;
								}
							}
						}
					}
				}

				// take a step
				for ( unsigned int ii = 0; ii < size; ++ii )
					solution[ ii ] += alpha * step[ ii ];
				progress = true;

				// if we hit a constraint, stop taking conjugate gradient steps
				if ( constraintIndex >= 0 ) {

					constraints[ constraintIndex ] = true;
					break;
				}

				// calculate the gradient
				#pragma omp parallel for schedule( static )
				for ( int ii = 0; ii < static_cast< int >( size ); ++ii ) {

					double accumulator = 0;
					double const* const row = &matrix[ 0 ] + ii * size;
					for ( unsigned int jj = 0; jj < size; ++jj )
						accumulator += row[ jj ] * solution[ jj ];
					gradient[ ii ] = -accumulator;
				}

				// project the gradient onto the equality constraints
				for ( unsigned int channel = 0; channel < 3; ++channel ) {

					double average = 0;
					for ( unsigned int ii = 256 * channel; ii < 256 * ( channel + 1 ); ++ii )
						average += gradient[ ii ];
					average /= 256;
					for ( unsigned int ii = 256 * channel; ii < 256 * ( channel + 1 ); ++ii )
						gradient[ ii ] -= average;
				}
				gradient[ 256 * 3 + minimumIndex ] = 0;
				gradient[ 256 * 3 + maximumIndex ] = 0;
				// project the gradient onto the inequality constraints
				for ( unsigned int channel = 0; channel < 3; ++channel ) {

					for ( bool done = false; ! done; ) {

						done = true;

						double numerator = 0;
						unsigned int denominator = 0;
						for ( unsigned int ii = 256 * channel; ii < 256 * ( channel + 1 ) ; ++ii ) {

							numerator += gradient[ ii ];
							++denominator;

							if ( ! constraints[ ii ] ) {

								if ( denominator > 1 ) {

									double const average = numerator / denominator;
									for ( unsigned int jj = 0; jj < denominator; ++jj )
										gradient[ ii - jj ] = average;
								}

								numerator   = 0;
								denominator = 0;
							}
						}
					}
				}

				// find the gradient's squared norm
				gradientNormSquared = 0;
				for ( unsigned int ii = 0; ii < size; ++ii )
					gradientNormSquared += Square( gradient[ ii ] );
				if ( std::sqrt( gradientNormSquared ) < epsilon )
					break;

				// update the step
				double beta = 0;
				if ( std::fabs( oldGradientNormSquared ) > std::numeric_limits< double >::epsilon() )
					beta = gradientNormSquared / oldGradientNormSquared;
				for ( unsigned int ii = 0; ii < size; ++ii )
					step[ ii ] = beta * step[ ii ] + gradient[ ii ];
			}

			// project the solution onto the equality constraints
			for ( unsigned int channel = 0; channel < 3; ++channel ) {

				double average = 0;
				for ( unsigned int ii = 256 * channel; ii < 256 * ( channel + 1 ); ++ii )
					average += solution[ ii ];
				average /= 256;
				for ( unsigned int ii = 256 * channel; ii < 256 * ( channel + 1 ); ++ii )
					solution[ ii ] -= average;
			}
			// project the solution onto the inequality constraints
			for ( unsigned int channel = 0; channel < 3; ++channel ) {

				for ( bool done = false; ! done; ) {

					done = true;

					double numerator = 0;
					unsigned int denominator = 0;
					for ( unsigned int ii = 256 * channel; ii < 256 * ( channel + 1 ) ; ++ii ) {

						numerator += solution[ ii ];
						++denominator;

						double const average = numerator / denominator;
						if ( ( ( ii & 255 ) < 255 ) && ( constraints[ ii ] || ( solution[ ii + 1 ] <= average ) ) ) {

							if ( ! constraints[ ii ] ) {

								done = false;
								constraints[ ii ] = true;
							}
						}
						else {

							if ( denominator > 1 ) {

								for ( unsigned int jj = 0; jj < denominator; ++jj )
									solution[ ii - jj ] = average;
							}

							numerator   = 0;
							denominator = 0;
						}
					}
				}
			}
		}
	}

	// scale the response curves to have uniform weighted intensity (before this point, they have uniform unweighted log-intensity)
	{	for ( unsigned int channel = 0; channel < 3; ++channel ) {

			double const maximum = solution[ 256 * channel + 255 ];
			for ( unsigned int ii = 256 * channel; ii < 256 * ( channel + 1 ); ++ii )
				solution[ ii ] -= maximum;
		}

		double scales[ 3 ];
		for ( unsigned int channel = 0; channel < 3; ++channel ) {

			double total = 0;
			for ( unsigned int ii = 0; ii < 256; ++ii )
				total += std::exp( solution[ 256 * channel + ii ] ) * MergeWeightFunction( static_cast< uint8_t >( ii ) );
			scales[ channel ] = std::log( total );
		}

		double maximumScale = 0;
		for ( unsigned int channel = 0; channel < 3; ++channel )
			maximumScale = std::max( maximumScale, scales[ channel ] );
		for ( unsigned int channel = 0; channel < 3; ++channel )
			scales[ channel ] -= maximumScale;

		for ( unsigned int channel = 0; channel < 3; ++channel )
			for ( unsigned int ii = 256 * channel; ii < 256 * ( channel + 1 ); ++ii )
				solution[ ii ] -= scales[ channel ];
	}

	// create the output image
	SharedArray< float > result( rows * columns * 3 );
	float maximum = -std::numeric_limits< float >::infinity();
	#pragma omp parallel
	{	float localMaximum = -std::numeric_limits< float >::infinity();

		#pragma omp for schedule( static )
		for ( int ii = 0; ii < static_cast< int >( rows * columns ); ++ii ) {

			for ( unsigned int channel = 0; channel < 3; ++channel ) {

				double numerator   = 0;
				double denominator = 0;
				for ( unsigned int jj = 0; jj < images.size(); ++jj ) {

					uint8_t const value = images[ jj ][ ii * 3 + channel ];

					double weight = 0;
					{	bool low  = false;
						bool high = false;
						if ( jj == maximumIndex )    // lightest image has full information on shadows
							low = true;
						if ( jj == minimumIndex )    // darkest image has full information on highlights
							high = true;
						weight = MergeWeightFunction( value, low, high );
					}
					weight += std::numeric_limits< float >::epsilon();

					numerator += weight * ( solution[ 256 * channel + value ] - solution[ 256 * 3 + jj ] );
					denominator += weight;
				}
				numerator /= denominator;

				result[ ii * 3 + channel ] = numerator;
				localMaximum = std::max( localMaximum, static_cast< float >( numerator ) );
			}
		}

		#pragma omp critical
		maximum = std::max( maximum, localMaximum );
	}

	// scale the output to [0,1]
	#pragma omp parallel for schedule( static )
	for ( int ii = 0; ii < static_cast< int >( rows * columns * 3 ); ++ii )
		result[ ii ] = std::exp( result[ ii ] - maximum );

	return result;
}
