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
	\file fix_noise.cc
	\brief implementation of FixNoise function
*/




#include "headers.hh"




//============================================================================
//    FixNoise function
//============================================================================


SharedArray< float > const FixNoise(
	SharedArray< float > const& image,
	unsigned int const rows,
	unsigned int const columns,
	double const pixels,
	int const radius,
	int const dimension,
	double const alpha,
	double const beta,
	int const window,
	double const epsilon
)
{
	if ( pixels <= 0 )
		throw std::runtime_error( "pixels must be positive" );
	if ( radius < 0 )
		throw std::runtime_error( "radius must be nonnegative" );
	if ( dimension <= 0 )
		throw std::runtime_error( "dimension must be positive" );
	if ( ( alpha <= 0 ) || ( alpha >= 1 ) )
		throw std::runtime_error( "alpha must be in the range (0,1)" );
	if ( ( beta <= 0 ) || ( beta >= 1 ) )
		throw std::runtime_error( "beta must be in the range (0,1)" );
	if ( window <= 0 )
		throw std::runtime_error( "window must be positive" );
	if ( ( window > static_cast< int >( rows ) ) || ( window > static_cast< int >( columns ) ) )
		throw std::runtime_error( "window cannot be larger than either image dimension" );

	if ( 3 * Square( 2 * radius + 1 ) < dimension )
		throw std::runtime_error( "3*(2*radius+1)^2 must be no smaller than dimension" );

	// treat the green channel as twice as important (since digital cameras generally have twice as many green sensors)
	double const channelScales[ 3 ] = { 1, 2, 1 };

	SharedArray< double > filter( 2 * radius + 1 );
	{	double total = 0;
		for ( unsigned int ii = 0; ii < static_cast< unsigned int >( 2 * radius + 1 ); ++ii ) {

			// Hamming window
			filter[ ii ] = 0.54 - 0.46 * std::cos( M_PI * ii / radius );

			// Tukey-ish window
			filter[ ii ] = std::min( filter[ ii ], 0.5 );

			total += filter[ ii ];
		}
		for ( unsigned int ii = 0; ii < static_cast< unsigned int >( 2 * radius + 1 ); ++ii )
			filter[ ii ] /= total;
	}

	SharedArray< unsigned int > reflectedRows( rows + 2 * radius );
	for ( unsigned int ii = 0; ii < rows + 2 * radius; ++ii ) {

		int jj = static_cast< int >( ii ) - static_cast< int >( radius );
		if ( jj < 0 )
			jj = -jj;
		else if ( jj >= static_cast< int >( rows ) )
			jj = 2 * static_cast< int >( rows - 1 ) - jj;
		assert( ( jj >= 0 ) && ( jj < static_cast< int >( rows ) ) );

		reflectedRows[ ii ] = jj;
	}

	SharedArray< unsigned int > reflectedColumns( columns + 2 * radius );
	for ( unsigned int ii = 0; ii < columns + 2 * radius; ++ii ) {

		int jj = static_cast< int >( ii ) - static_cast< int >( radius );
		if ( jj < 0 )
			jj = -jj;
		else if ( jj >= static_cast< int >( columns ) )
			jj = 2 * static_cast< int >( columns - 1 ) - jj;
		assert( ( jj >= 0 ) && ( jj < static_cast< int >( columns ) ) );

		reflectedColumns[ ii ] = jj;
	}

	unsigned int const size = 3 * Square( 2 * radius + 1 );

	SharedArray< double > mean( size );
	SharedArray< double > covariance( Square( size ) );
	std::fill( &mean[ 0 ], &mean[ 0 ] + size, 0 );
	std::fill( &covariance[ 0 ], &covariance[ 0 ] + Square( size ), 0 );

	unsigned int truePixels = 0;
	{	unsigned int const sampleRows    = std::max( 1u, std::min( rows,    static_cast< unsigned int >( std::ceil( std::sqrt( pixels * rows    / columns ) ) ) ) );
		unsigned int const sampleColumns = std::max( 1u, std::min( columns, static_cast< unsigned int >( std::ceil( std::sqrt( pixels * columns / rows    ) ) ) ) );

		#pragma omp parallel
		{	SharedArray< float > vector( size );

			SharedArray< double > innerMean( size );
			SharedArray< double > innerCovariance( size * ( size + 1 ) / 2 );
			std::fill( &innerMean[ 0 ], &innerMean[ 0 ] + size, 0 );
			std::fill( &innerCovariance[ 0 ], &innerCovariance[ 0 ] + size * ( size + 1 ) / 2, 0 );

			#pragma omp for schedule( static )
			for ( int sampleRow = 0; sampleRow < static_cast< int >( sampleRows ); ++sampleRow ) {

				unsigned int const row = ( rows * ( 2 * sampleRow + 1 ) ) / ( 2 * sampleRows );

				for ( int sampleColumn = 0; sampleColumn < static_cast< int >( sampleColumns ); ++sampleColumn ) {

					unsigned int const column = ( columns * ( 2 * sampleColumn + 1 ) ) / ( 2 * sampleColumns );

					unsigned int index = 0;
					for ( unsigned int ii = 0; ii < static_cast< unsigned int >( 2 * radius + 1 ); ++ii ) {

						float const* pRow = &image[ reflectedRows[ row + ii ] * columns * 3 ];
						for ( unsigned int jj = 0; jj < static_cast< unsigned int >( 2 * radius + 1 ); ++jj ) {

							float const* pColumn = &pRow[ reflectedColumns[ column + jj ] * 3 ];
							for ( unsigned int channel = 0; channel < 3; ++channel, ++index )
								vector[ index ] = channelScales[ channel ] * pColumn[ channel ];
						}
					}

					index = 0;
					for ( unsigned int ii = 0; ii < size; ++ii ) {

						innerMean[ ii ] += vector[ ii ];
						for ( unsigned int jj = 0; jj <= ii; ++jj, ++index )
							innerCovariance[ index ] += vector[ ii ] * vector[ jj ];
					}
				}
			}

			#pragma omp critical
			{	unsigned int index = 0;
				for ( unsigned int ii = 0; ii < size; ++ii ) {

					mean[ ii ] += innerMean[ ii ];

					for ( unsigned int jj = 0; jj < ii; ++jj, ++index ) {

						double const value = innerCovariance[ index ];
						covariance[ ii * size + jj ] += value;
						covariance[ jj * size + ii ] += value;
					}

					covariance[ ii * size + ii ] += innerCovariance[ index ];
					++index;
				}
			}
		}

		truePixels = sampleRows * sampleColumns;
	}

	for ( unsigned int ii = 0; ii < size; ++ii ) {

		mean[ ii ] /= truePixels;
		for ( unsigned int jj = 0; jj < size; ++jj )
			covariance[ ii * size + jj ] /= truePixels;
	}
	for ( unsigned int ii = 0; ii < size; ++ii )
		for ( unsigned int jj = 0; jj < size; ++jj )
			covariance[ ii * size + jj ] -= mean[ ii ] * mean[ jj ];

	SharedArray< double > vectors( dimension * size );

	// subspace iteration method
	{	SharedArray< double > matrix( dimension * size );
		SharedArray< double > product( dimension * dimension );

		for ( unsigned int ii = 0; ii < dimension * size; ++ii )
			vectors[ ii ] = rand() / static_cast< double >( RAND_MAX );

		for ( ; ; ) {

			// Graham-Schmidt orthonormalization
			for ( unsigned int ii = 0; ii < static_cast< unsigned int >( dimension ); ++ii ) {

				double* pTargetVector = &vectors[ ii * size ];
				double* pSourceVector = &vectors[ 0 ];

				for ( unsigned int jj = 0; jj < ii; ++jj ) {

					double dot = 0;
					for ( unsigned int kk = 0; kk < size; ++kk )
						dot += pTargetVector[ kk ] * pSourceVector[ kk ];
					for ( unsigned int kk = 0; kk < size; ++kk )
						pTargetVector[ kk ] -= dot * pSourceVector[ kk ];

					pSourceVector += size;
				}

				double norm = 0;
				for ( unsigned int kk = 0; kk < size; ++kk )
					norm += Square( pTargetVector[ kk ] );
				norm = std::sqrt( norm );
				for ( unsigned int kk = 0; kk < size; ++kk )
					pTargetVector[ kk ] /= norm;
			}

			// matrix = vectors * covariance
			#pragma omp parallel for schedule( static )
			for ( int ii = 0; ii < dimension; ++ii ) {

				double const* pVector = &vectors[ ii * size ];
				double const* pCovariance = &covariance[ 0 ];
				double* pMatrix = &matrix[ ii * size ];
				for ( unsigned int jj = 0; jj < size; ++jj ) {

					double accumulator = 0;
					for ( unsigned int kk = 0; kk < size; ++kk )
						accumulator += pVector[ kk ] * pCovariance[ kk ];
					*pMatrix = accumulator;

					pCovariance += size;
					++pMatrix;
				}
			}

			// product = matrix * vectors'
			#pragma omp parallel for schedule( static )
			for ( int ii = 0; ii < dimension; ++ii ) {

				double const* pMatrix = &matrix[ ii * size ];
				double const* pVectors = &vectors[ 0 ];
				double* pProduct = &product[ ii * dimension ];
				for ( int jj = 0; jj < dimension; ++jj ) {

					double accumulator = 0;
					for ( unsigned int kk = 0; kk < size; ++kk )
						accumulator += pMatrix[ kk ] * pVectors[ kk ];
					*pProduct = accumulator;

					pVectors += size;
					++pProduct;
				}
			}

			// find the squared Frobenius norm (which upper bounds the squared spectral norm) of matrix - product * vectors
			double squaredNorm = 0;
			#pragma omp parallel for schedule( static ) reduction( + : squaredNorm )
			for ( int ii = 0; ii < dimension; ++ii ) {

				double const* pProduct = &product[ ii * dimension ];
				double const* pMatrix = &matrix[ ii * size ];
				for ( unsigned int jj = 0; jj < size; ++jj ) {

					double accumulator = 0;
					for ( int kk = 0; kk < dimension; ++kk )
						accumulator += pProduct[ kk ] * vectors[ kk * size + jj ];
					squaredNorm += Square( accumulator - *pMatrix );

					++pMatrix;
				}
			}
			if ( std::sqrt( squaredNorm ) <= epsilon )
				break;

			std::copy( &matrix[ 0 ], &matrix[ 0 ] + dimension * size, &vectors[ 0 ] );
		}
	}

	SharedArray< float > descriptors( rows * columns * dimension );
	#pragma omp parallel
	{	SharedArray< float > vector( size );

		#pragma omp for schedule( static )
		for ( int row = 0; row < static_cast< int >( rows ); ++row ) {

			float* pDescriptor = &descriptors[ row * columns * dimension ];
			for ( int column = 0; column < static_cast< int >( columns ); ++column ) {

				unsigned int index = 0;
				for ( unsigned int ii = 0; ii < static_cast< unsigned int >( 2 * radius + 1 ); ++ii ) {

					float const* pRow = &image[ reflectedRows[ row + ii ] * columns * 3 ];
					for ( unsigned int jj = 0; jj < static_cast< unsigned int >( 2 * radius + 1 ); ++jj ) {

						float const* pColumn = &pRow[ reflectedColumns[ column + jj ] * 3 ];
						for ( unsigned int channel = 0; channel < 3; ++channel, ++index )
							vector[ index ] = channelScales[ channel ] * pColumn[ channel ];
					}
				}

				for ( unsigned int ii = 0; ii < static_cast< unsigned int >( dimension ); ++ii ) {

					double accumulator = 0;
					for ( unsigned int jj = 0; jj < size; ++jj )
						accumulator += vectors[ ii * size + jj ] * vector[ jj ];
					*pDescriptor = accumulator;
					++pDescriptor;
				}
			}
		}
	}

	// find squared differences between many adjacent pixels
	std::vector< double > squaredDifferences;
	{	unsigned int const sampleRows    = std::max( 1u, std::min( rows,    static_cast< unsigned int >( std::ceil( std::sqrt( pixels * rows    / columns ) ) ) ) );
		unsigned int const sampleColumns = std::max( 1u, std::min( columns, static_cast< unsigned int >( std::ceil( std::sqrt( pixels * columns / rows    ) ) ) ) );

		{	SharedArray< float > vector( size );

			for ( unsigned int sampleRow = 0; sampleRow < sampleRows; ++sampleRow ) {

				unsigned int const row = ( rows * ( 2 * sampleRow + 1 ) ) / ( 2 * sampleRows );

				for ( unsigned int sampleColumn = 0; sampleColumn < sampleColumns; ++sampleColumn ) {

					unsigned int const column = ( columns * ( 2 * sampleColumn + 1 ) ) / ( 2 * sampleColumns );

					if ( ( row + 1 < rows ) && ( column + 1 < columns ) ) {

						float const* pSourceDescriptor = &descriptors[ ( row * columns + column ) * dimension ];

						{	float const* pTargetDescriptor = &descriptors[ ( ( row + 1 ) * columns + column ) * dimension ];
							double squaredDifference = 0;
							for ( unsigned int ii = 0; ii < static_cast< unsigned int >( dimension ); ++ii )
								squaredDifference += Square( pSourceDescriptor[ ii ] - pTargetDescriptor[ ii ] );
							squaredDifferences.push_back( squaredDifference );
						}

						{	float const* pTargetDescriptor = &descriptors[ ( row * columns + ( column + 1 ) ) * dimension ];
							double squaredDifference = 0;
							for ( unsigned int ii = 0; ii < static_cast< unsigned int >( dimension ); ++ii )
								squaredDifference += Square( pSourceDescriptor[ ii ] - pTargetDescriptor[ ii ] );
							squaredDifferences.push_back( squaredDifference );
						}
					}
				}
			}
		}

		truePixels = sampleRows * sampleColumns;
	}
	// choose sigma such that requested noise percentile has requested weight
	std::sort( squaredDifferences.begin(), squaredDifferences.end() );
	double const negativeGamma = std::log( beta ) / squaredDifferences[ std::floor( squaredDifferences.size() * alpha ) ];

	SharedArray< float > result( rows * columns * 3 );
	// ignore pixels which are weighted less than window^-2
	double const cutoff = -2 * std::log( window ) / negativeGamma;

	#pragma omp parallel for schedule( static )
	for ( int row = 0; row < static_cast< int >( rows ); ++row ) {

		float* pResult = &result[ row * columns * 3 ];
		float const* pSourceDescriptor = &descriptors[ row * columns * dimension ];
		for ( int column = 0; column < static_cast< int >( columns ); ++column ) {

			double numerators[ 3 ] = { 0, 0, 0 };
			double denominator = 0;

			int lowerRow = row - window / 2;
			int upperRow = lowerRow + window;
			if ( lowerRow < 0 ) {

				upperRow -= lowerRow;
				lowerRow = 0;
			}
			else if ( upperRow > static_cast< int >( rows ) ) {

				lowerRow -= upperRow - rows;
				upperRow = rows;
			}
			assert( ( lowerRow >= 0 ) && ( upperRow <= static_cast< int >( rows ) ) );

			int lowerColumn = column - window / 2;
			int upperColumn = lowerColumn + window;
			if ( lowerColumn < 0 ) {

				upperColumn -= lowerColumn;
				lowerColumn = 0;
			}
			else if ( upperColumn > static_cast< int >( columns ) ) {

				lowerColumn -= upperColumn - columns;
				upperColumn = columns;
			}
			assert( ( lowerColumn >= 0 ) && ( upperColumn <= static_cast< int >( columns ) ) );

			for ( unsigned int ii = lowerRow; ii < static_cast< unsigned int >( upperRow ); ++ii ) {

				float const* pImage = &image[ ( ii * columns + lowerColumn ) * 3 ];
				float const* pTargetDescriptor = &descriptors[ ( ii * columns + lowerColumn ) * dimension ];
				for ( unsigned int jj = lowerColumn; jj < static_cast< unsigned int >( upperColumn ); ++jj ) {

					double squaredDifference = 0;
					for ( unsigned int kk = 0; ( kk < static_cast< unsigned int >( dimension ) ) && ( squaredDifference < cutoff ); ++kk )
						squaredDifference += Square( pSourceDescriptor[ kk ] - pTargetDescriptor[ kk ] );

					if ( squaredDifference < cutoff ) {

						double const weight = FastExp( negativeGamma * squaredDifference );
						for ( unsigned int channel = 0; channel < 3; ++channel )
							numerators[ channel ] += weight * pImage[ channel ];
						denominator += weight;
					}

					pImage += 3;
					pTargetDescriptor += dimension;
				}
			}

			if ( denominator > 0 ) {

				for ( unsigned int channel = 0; channel < 3; ++channel )
					numerators[ channel ] /= denominator;
			}
			for ( unsigned int channel = 0; channel < 3; ++channel ) {

				*pResult += numerators[ channel ];
				++pResult;
			}

			pSourceDescriptor += dimension;
		}
	}

	return result;
}
