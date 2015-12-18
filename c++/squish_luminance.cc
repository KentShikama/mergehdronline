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
	\file squish_luminance.cc
	\brief implementation of SquishLuminance function
*/




#include "headers.hh"




namespace {




//============================================================================
//    FindHalfSizeImage function
//============================================================================


SharedArray< float > const FindHalfSizeImage(
	SharedArray< float > const& luminance,
	unsigned int const rows,
	unsigned int const columns
)
{
	unsigned int const destinationRows    = ( rows    + 1 ) / 2;
	unsigned int const destinationColumns = ( columns + 1 ) / 2;

	//double const sigma = 1.0 / std::sqrt( 2.0 );
	double const sigma = 1.0;
	unsigned int const size = static_cast< unsigned int >( std::ceil( sigma * 3 ) );

	SharedArray< double > filter( size );
	for ( unsigned int ii = 0; ii < size; ++ii )
		filter[ ii ] = std::exp( -0.5 / Square( ii / sigma ) );

	SharedArray< float > temp( rows * destinationColumns );
	#pragma omp parallel for schedule( static )
	for ( int ii = 0; ii < static_cast< int >( rows ); ++ii ) {

		float const* const pSourceBegin = &luminance[ 0 ] + ii * columns;
		float* pDestination = &temp[ 0 ] + ii * destinationColumns;
		for ( unsigned int jj = 0; jj < destinationColumns; ++jj, ++pDestination ) {

			double numerator   = 0;
			double denominator = 0;

			{	float const* pSource = pSourceBegin + ( jj * 2 );
				float const* const pSourceEnd = pSourceBegin;
				for ( unsigned int kk = 0; ( kk < size ) && ( pSource >= pSourceEnd ); ++kk, --pSource ) {

					numerator += filter[ kk ] * *pSource;
					denominator += filter[ kk ];
				}
			}

			{	float const* pSource = pSourceBegin + ( jj * 2 + 1 );
				float const* const pSourceEnd = pSourceBegin + ( columns - 1 );
				for ( unsigned int kk = 1; ( kk < size ) && ( pSource <= pSourceEnd ); ++kk, ++pSource ) {

					numerator += filter[ kk ] * *pSource;
					denominator += filter[ kk ];
				}
			}

			if ( denominator > 0 )
				numerator /= denominator;
			*pDestination = numerator;
		}
		assert( pDestination == &temp[ 0 ] + ( ii + 1 ) * destinationColumns );
	}

	SharedArray< float > result( destinationRows * destinationColumns );
	#pragma omp parallel for schedule( static )
	for ( int ii = 0; ii < static_cast< int >( destinationColumns ); ++ii ) {

		float const* const pSourceBegin = &temp[ 0 ] + ii;
		float* pDestination = &result[ 0 ] + ii;
		for ( unsigned int jj = 0; jj < destinationRows; ++jj, pDestination += destinationColumns ) {

			double numerator   = 0;
			double denominator = 0;

			{	float const* pSource = pSourceBegin + ( jj * 2 ) * destinationColumns;
				float const* const pSourceEnd = pSourceBegin;
				for ( unsigned int kk = 0; ( kk < size ) && ( pSource >= pSourceEnd ); ++kk, pSource -= destinationColumns ) {

					numerator += filter[ kk ] * *pSource;
					denominator += filter[ kk ];
				}
			}

			{	float const* pSource = pSourceBegin + ( jj * 2 + 1 ) * destinationColumns;
				float const* const pSourceEnd = pSourceBegin + ( rows - 1 ) * destinationColumns;;
				for ( unsigned int kk = 1; ( kk < size ) && ( pSource <= pSourceEnd ); ++kk, pSource += destinationColumns ) {

					numerator += filter[ kk ] * *pSource;
					denominator += filter[ kk ];
				}
			}

			if ( denominator > 0 )
				numerator /= denominator;
			*pDestination = numerator;
		}
	}

	return result;
}




//============================================================================
//    FindGradientAttenuation function
//============================================================================


SharedArray< float > const FindGradientAttenuation(
	SharedArray< float > const& luminance,
	unsigned int const rows,
	unsigned int const columns,
	double const alpha,
	double const beta,
	double const delta,
	double const theta
)
{
	if ( ( alpha <= 0 ) || ( alpha >= 1 ) )
		throw std::runtime_error( "alpha must be in the range (0,1)" );
	if ( ( beta <= 0 ) || ( beta >= 1 ) )
		throw std::runtime_error( "beta must be in the range (0,1)" );
	if ( delta < 1 )
		throw std::runtime_error( "delta must be at least 1" );
	if ( ( theta < 0 ) || ( theta > 1 ) )
		throw std::runtime_error( "theta must be in the range [0,1]" );

	if ( ( rows < 32 ) || ( columns < 32 ) )
		throw std::runtime_error( "image must be at least 32x32" );

	std::vector< SharedArray< float > > pyramid;
	std::vector< std::pair< unsigned int, unsigned int > > pyramidSizes;

	// create a pyramid of gradient magnitudes
	{	SharedArray< float > source = luminance;
		unsigned int sourceRows    = rows;
		unsigned int sourceColumns = columns;
		double sourceDenominator = 1;

		for ( ; ; ) {

			SharedArray< float > destination( sourceRows * sourceColumns );

			{	float const* pSource      = &source[ 0 ];
				float*       pDestination = &destination[ 0 ];
				{	double const horizontalGradient = ( *( pSource +             1 ) - *pSource ) / sourceDenominator;
					double const verticalGradient   = ( *( pSource + sourceColumns ) - *pSource ) / sourceDenominator;
					*pDestination = std::sqrt( Square( horizontalGradient ) + Square( verticalGradient ) );
				}
				++pSource;
				++pDestination;
				for ( unsigned int jj = 1; jj + 1 < sourceColumns; ++jj, ++pSource, ++pDestination ) {

					double const horizontalGradient = ( *( pSource +             1 ) - *( pSource - 1 ) ) / ( 2 * sourceDenominator );
					double const verticalGradient   = ( *( pSource + sourceColumns ) -   *pSource       ) /       sourceDenominator;
					*pDestination = std::sqrt( Square( horizontalGradient ) + Square( verticalGradient ) );
				}
				{	double const horizontalGradient = (   *pSource                   - *( pSource - 1 ) ) / sourceDenominator;
					double const verticalGradient   = ( *( pSource + sourceColumns ) -   *pSource       ) / sourceDenominator;
					*pDestination = std::sqrt( Square( horizontalGradient ) + Square( verticalGradient ) );
				}
				++pSource;
				++pDestination;
			}
			#pragma omp parallel for schedule( static )
			for ( int ii = 1; ii < static_cast< int >( sourceRows ) - 1; ++ii ) {

				float const* pSource      = &source[ 0 ]      + ii * sourceColumns;
				float*       pDestination = &destination[ 0 ] + ii * sourceColumns;
				{	double const horizontalGradient = ( *( pSource +             1 ) -   *pSource                   ) /       sourceDenominator;
					double const verticalGradient   = ( *( pSource + sourceColumns ) - *( pSource - sourceColumns ) ) / ( 2 * sourceDenominator );
					*pDestination = std::sqrt( Square( horizontalGradient ) + Square( verticalGradient ) );
				}
				++pSource;
				++pDestination;
				for ( unsigned int jj = 1; jj + 1 < sourceColumns; ++jj, ++pSource, ++pDestination ) {

					double const horizontalGradient = ( *( pSource +             1 ) - *( pSource -             1 ) ) / ( 2 * sourceDenominator );
					double const verticalGradient   = ( *( pSource + sourceColumns ) - *( pSource - sourceColumns ) ) / ( 2 * sourceDenominator );
					*pDestination = std::sqrt( Square( horizontalGradient ) + Square( verticalGradient ) );
				}
				{	double const horizontalGradient = (   *pSource                   - *( pSource -             1 ) ) /       sourceDenominator;
					double const verticalGradient   = ( *( pSource + sourceColumns ) - *( pSource - sourceColumns ) ) / ( 2 * sourceDenominator );
					*pDestination = std::sqrt( Square( horizontalGradient ) + Square( verticalGradient ) );
				}
				++pSource;
				++pDestination;
			}
			{	float const* pSource      = &source[ 0 ]      + ( sourceRows - 1 ) * sourceColumns;
				float*       pDestination = &destination[ 0 ] + ( sourceRows - 1 ) * sourceColumns;
				{	double const horizontalGradient = ( *( pSource + 1 ) -   *pSource                   ) / sourceDenominator;
					double const verticalGradient   = (   *pSource       - *( pSource - sourceColumns ) ) / sourceDenominator;
					*pDestination = std::sqrt( Square( horizontalGradient ) + Square( verticalGradient ) );
				}
				++pSource;
				++pDestination;
				for ( unsigned int jj = 1; jj + 1 < sourceColumns; ++jj, ++pSource, ++pDestination ) {

					double const horizontalGradient = ( *( pSource + 1 ) - *( pSource -             1 ) ) / ( 2 * sourceDenominator );
					double const verticalGradient   = (   *pSource       - *( pSource - sourceColumns ) ) /       sourceDenominator;
					*pDestination = std::sqrt( Square( horizontalGradient ) + Square( verticalGradient ) );
				}
				{	double const horizontalGradient = ( *pSource - *( pSource -             1 ) ) / sourceDenominator;
					double const verticalGradient   = ( *pSource - *( pSource - sourceColumns ) ) / sourceDenominator;
					*pDestination = std::sqrt( Square( horizontalGradient ) + Square( verticalGradient ) );
				}
				++pSource;
				++pDestination;
			}

			pyramid.push_back( destination );
			pyramidSizes.push_back( std::pair< unsigned int, unsigned int >( sourceRows, sourceColumns ) );

			if ( ( sourceRows < 63 ) || ( sourceColumns < 63 ) )
				break;

			source = FindHalfSizeImage( source, sourceRows, sourceColumns );
			sourceRows    = ( sourceRows    + 1 ) / 2;
			sourceColumns = ( sourceColumns + 1 ) / 2;
			sourceDenominator *= 2;
		}
	}

	double overallAverageGradientMagnitude = 0;
	{	unsigned int denominator = 0;

		unsigned int const size = pyramid.size();
		assert( size == pyramidSizes.size() );
		for ( unsigned int ii = 0; ii < size; ++ii ) {

			SharedArray< float > const& source = pyramid[ ii ];
			unsigned int sourceRows    = pyramidSizes[ ii ].first;
			unsigned int sourceColumns = pyramidSizes[ ii ].second;

			#pragma omp parallel for schedule( static ) reduction( + : overallAverageGradientMagnitude )
			for ( int jj = 0; jj < static_cast< int >( sourceRows * sourceColumns ); ++jj )
				overallAverageGradientMagnitude += source[ jj ];

			denominator += sourceRows * sourceColumns;
		}
		overallAverageGradientMagnitude /= denominator;
	}

	SharedArray< float > result;
	unsigned int resultRows    = 0;
	unsigned int resultColumns = 0;

	// calculate the gradient attenuation factors
	while ( ! pyramid.empty() ) {

		SharedArray< float > const& source = pyramid.back();
		unsigned int sourceRows    = pyramidSizes.back().first;
		unsigned int sourceColumns = pyramidSizes.back().second;

		double sourceAverageGradientMagnitude = 0;
		#pragma omp parallel for schedule( static ) reduction( + : sourceAverageGradientMagnitude )
		for ( int ii = 0; ii < static_cast< int >( sourceRows * sourceColumns ); ++ii )
			sourceAverageGradientMagnitude += source[ ii ];
		sourceAverageGradientMagnitude /= sourceRows * sourceColumns;

		double const averageGradientMagnitude = theta * sourceAverageGradientMagnitude + ( 1 - theta ) * overallAverageGradientMagnitude;

		SharedArray< float > destination( sourceRows * sourceColumns );
		#pragma omp parallel for schedule( static )
		for ( int ii = 0; ii < static_cast< int >( sourceRows * sourceColumns ); ++ii ) {

			double attenuation = std::numeric_limits< double >::infinity();
			if ( source[ ii ] > 0 )
				attenuation = std::pow( alpha * averageGradientMagnitude / source[ ii ], beta );

			// **TODO **FIXME: is this the right thing to do when we have a gradient of zero?
			destination[ ii ] = std::min( attenuation, delta );
		}

		if ( result != NULL ) {

			#pragma omp parallel for schedule( static )
			for ( int ii = 0; ii < static_cast< int >( sourceRows ); ++ii ) {

				unsigned int iiIndex = ii / 2;
				if ( iiIndex + 1 >= resultRows )
					iiIndex = resultRows - 2;
				assert( iiIndex + 1 < resultRows );

				double const rowLambda = 0.5 * ( ( iiIndex + 1 ) * 2.0 - ii );

				float const* const pSource      = &result[ 0 ]      + iiIndex * resultColumns;
				float*             pDestination = &destination[ 0 ] + ii      * sourceColumns;
				for ( unsigned int jj = 0; jj < sourceColumns; ++jj, ++pDestination ) {

					unsigned int jjIndex = jj / 2;
					if ( jjIndex + 1 >= resultColumns )
						jjIndex = resultColumns - 2;
					assert( jjIndex + 1 < resultColumns );

					double const columnLambda = 0.5 * ( ( jjIndex + 1 ) * 2.0 - jj );

					double value = 0;
					value += *( pSource + jjIndex                     ) *       rowLambda   *       columnLambda;
					value += *( pSource + jjIndex                 + 1 ) *       rowLambda   * ( 1 - columnLambda );
					value += *( pSource + jjIndex + resultColumns     ) * ( 1 - rowLambda ) *       columnLambda;
					value += *( pSource + jjIndex + resultColumns + 1 ) * ( 1 - rowLambda ) * ( 1 - columnLambda );

					*pDestination *= value;
				}
			}
		}

		result = destination;
		resultRows    = sourceRows;
		resultColumns = sourceColumns;

		pyramid.pop_back();
		pyramidSizes.pop_back();
	}

	assert( resultRows    == rows    );
	assert( resultColumns == columns );

	return result;
}




//============================================================================
//    FindPoissonRHS function
//============================================================================


SharedArray< float > const FindPoissonRHS(
	SharedArray< float > const& luminance,
	SharedArray< float > const& attenuation,
	unsigned int const rows,
	unsigned int const columns
)
{
	SharedArray< float > horizontalGradients( rows * columns );
	SharedArray< float > verticalGradients(   rows * columns );

	#pragma omp parallel for schedule( static )
	for ( int ii = 0; ii < static_cast< int >( rows ) - 1; ++ii ) {

		float const* pSource      = &luminance[ 0 ]           + ii * columns;
		float const* pAttenuation = &attenuation[ 0 ]         + ii * columns;
		float*       pHorizontal  = &horizontalGradients[ 0 ] + ii * columns;
		float*       pVertical    = &verticalGradients[ 0 ]   + ii * columns;
		for ( unsigned int jj = 0; jj < columns - 1; ++jj, ++pSource, ++pAttenuation, ++pHorizontal, ++pVertical ) {

			*pHorizontal = *pAttenuation * ( *( pSource +       1 ) - *pSource );
			*pVertical   = *pAttenuation * ( *( pSource + columns ) - *pSource );
		}
		{	*pHorizontal = 0;
			*pVertical   = *pAttenuation * ( *( pSource + columns ) - *pSource );
		}
	}
	{	float const* pSource      = &luminance[ 0 ]           + ( rows - 1 ) * columns;
		float const* pAttenuation = &attenuation[ 0 ]         + ( rows - 1 ) * columns;
		float*       pHorizontal  = &horizontalGradients[ 0 ] + ( rows - 1 ) * columns;
		float*       pVertical    = &verticalGradients[ 0 ]   + ( rows - 1 ) * columns;
		for ( unsigned int jj = 0; jj < columns - 1; ++jj, ++pSource, ++pAttenuation, ++pHorizontal, ++pVertical ) {

			*pHorizontal = *pAttenuation * ( *( pSource + 1 ) - *pSource );
			*pVertical   = 0;
		}
		{	*pHorizontal = 0;
			*pVertical   = 0;
		}
	}

	SharedArray< float > result( rows * columns );

	{	float const* pHorizontal  = &horizontalGradients[ 0 ];
		float const* pVertical    = &verticalGradients[ 0 ];
		float*       pDestination = &result[ 0 ];
		*pDestination = *pHorizontal + *pVertical;
		++pHorizontal;
		++pVertical;
		++pDestination;
		for ( unsigned int jj = 1; jj < columns; ++jj, ++pHorizontal, ++pVertical, ++pDestination )
			*pDestination = ( *pHorizontal - *( pHorizontal - 1 ) ) + *pVertical;
	}
	#pragma omp parallel for schedule( static )
	for ( int ii = 1; ii < static_cast< int >( rows ); ++ii ) {

		float const* pHorizontal  = &horizontalGradients[ 0 ] + ii * columns;
		float const* pVertical    = &verticalGradients[ 0 ]   + ii * columns;
		float*       pDestination = &result[ 0 ]              + ii * columns;
		*pDestination = *pHorizontal + ( *pVertical - *( pVertical - columns ) );
		++pHorizontal;
		++pVertical;
		++pDestination;
		for ( unsigned int jj = 1; jj < columns; ++jj, ++pHorizontal, ++pVertical, ++pDestination )
			*pDestination = ( *pHorizontal - *( pHorizontal - 1 ) ) + ( *pVertical - *( pVertical - columns ) );
	}

	return result;
}




}    // anonymous namespace




//============================================================================
//    SquishLuminance function
//============================================================================


/*
	Fattal, Lischinski, Werman. "Gradient Domain High Dynamic Range
	Compression".
*/
void SquishLuminance(
	SharedArray< float > const& luminance,    // input/output parameter
	unsigned int const rows,
	unsigned int const columns,
	double const epsilon,    // termination threshold for SolvePoisson
	double const alpha,
	double const beta,
	double const delta,
	double const theta
)
{
	SolvePoisson(
		luminance,
		FindPoissonRHS(
			luminance,
			FindGradientAttenuation(
				luminance,
				rows,
				columns,
				alpha,
				beta,
				delta,
				theta
			),
			rows,
			columns
		),
		rows,
		columns,
		epsilon
	);
}
