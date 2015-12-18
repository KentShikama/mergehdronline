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
	\file fix_chromatic_aberration.cc
	\brief implementation of FixChromaticAberration function
*/




#include "headers.hh"




namespace {




//============================================================================
//    FindEdges function
//============================================================================


SharedArray< float > const FindEdges(
	SharedArray< float > const& image,
	unsigned int const rows,
	unsigned int const columns,
	double const sigma1 = 1,
	double const sigma2 = 16,
	double const lambda = 0.1
)
{
	if ( sigma1 <= 0 )
		throw std::runtime_error( "sigma1 must be positive" );
	if ( lambda <= 0 )
		throw std::runtime_error( "lambda must be positive" );

	SharedArray< float > smoothed( rows * columns * 3 );

	// find smoothed log-image
	{	unsigned int const size = static_cast< unsigned int >( std::ceil( sigma1 * 3 ) );

		SharedArray< double > filter( size );
		for ( unsigned int ii = 0; ii < size; ++ii )
			filter[ ii ] = std::exp( -0.5 * Square( ii / sigma1 ) );

		SharedArray< float > half( rows * columns );

		for ( unsigned int channel = 0; channel < 3; ++channel ) {

			#pragma omp parallel for schedule( static )
			for ( int ii = 0; ii < static_cast< int >( rows ); ++ii ) {

				for ( unsigned int jj = 0; jj < columns; ++jj ) {

					double numerator   = 0;
					double denominator = 0;

					{	float const* pSource = &image[ 0 ] + ( ii * columns + jj ) * 3 + channel;
						float const* const pSourceEnd = &image[ 0 ] + ( ii * columns ) * 3 + channel;
						for ( unsigned int kk = 0; ( kk < size ) && ( pSource >= pSourceEnd ); ++kk, pSource -= 3 ) {

							numerator += filter[ kk ] * std::log( *pSource );
							denominator += filter[ kk ];
						}
					}

					{	float const* pSource = &image[ 0 ] + ( ii * columns + ( jj + 1 ) ) * 3 + channel;
						float const* const pSourceEnd = &image[ 0 ] + ( ii * columns + ( columns - 1 ) ) * 3 + channel;
						for ( unsigned int kk = 1; ( kk < size ) && ( pSource <= pSourceEnd ); ++kk, pSource += 3 ) {

							numerator += filter[ kk ] * std::log( *pSource );
							denominator += filter[ kk ];
						}
					}

					half[ ii * columns + jj ] = numerator / denominator;
				}
			}

			#pragma omp parallel for schedule( static )
			for ( int ii = 0; ii < static_cast< int >( rows ); ++ii ) {

				for ( unsigned int jj = 0; jj < columns; ++jj ) {

					double numerator   = 0;
					double denominator = 0;

					{	float const* pSource = &half[ 0 ] + ( ii * columns + jj );
						float const* const pSourceEnd = &half[ 0 ] + jj;
						for ( unsigned int kk = 0; ( kk < size ) && ( pSource >= pSourceEnd ); ++kk, pSource -= columns ) {

							numerator += filter[ kk ] * *pSource;
							denominator += filter[ kk ];
						}
					}

					{	float const* pSource = &half[ 0 ] + ( ( ii + 1 ) * columns + jj );
						float const* const pSourceEnd = &half[ 0 ] + ( ( rows - 1 ) * columns + jj );
						for ( unsigned int kk = 1; ( kk < size ) && ( pSource <= pSourceEnd ); ++kk, pSource += columns ) {

							numerator += filter[ kk ] * *pSource;
							denominator += filter[ kk ];
						}
					}

					smoothed[ ( ii * columns + jj ) * 3 + channel ] = numerator / denominator;
				}
			}
		}
	}

	SharedArray< float > result( rows * columns * 3 );

	// find Canny edge strengths
	#pragma omp parallel for schedule( static )
	for ( int ii = 0; ii < static_cast< int >( rows ); ++ii ) {

		int rowIndex2 = ii;
		int rowIndex1 = rowIndex2 - 1;
		int rowIndex3 = rowIndex2 + 1;
		if ( rowIndex1 < 0 )
			rowIndex1 = rowIndex2;
		if ( rowIndex3 >= static_cast< int >( rows ) )
			rowIndex3 = rowIndex2;
		rowIndex1 *= columns * 3;
		rowIndex2 *= columns * 3;
		rowIndex3 *= columns * 3;

		for ( unsigned int jj = 0; jj < columns; ++jj ) {

			int columnIndex2 = jj;
			int columnIndex1 = columnIndex2 - 1;
			int columnIndex3 = columnIndex2 + 1;
			if ( columnIndex1 < 0 )
				columnIndex1 = columnIndex2;
			if ( columnIndex3 >= static_cast< int >( columns ) )
				columnIndex3 = columnIndex2;
			columnIndex1 *= 3;
			columnIndex2 *= 3;
			columnIndex3 *= 3;

			for ( unsigned int channel = 0; channel < 3; ++channel ) {

				double const dx = -0.5 * smoothed[ rowIndex2 + columnIndex1 + channel ] + 0.5 * smoothed[ rowIndex2 + columnIndex3 + channel ];
				double const dy = -0.5 * smoothed[ rowIndex1 + columnIndex2 + channel ] + 0.5 * smoothed[ rowIndex3 + columnIndex2 + channel ];

				double const dxdx = smoothed[ rowIndex2 + columnIndex1 + channel ] - 2 * smoothed[ rowIndex2 + columnIndex2 + channel ] + smoothed[ rowIndex2 + columnIndex3 + channel ];
				double const dydy = smoothed[ rowIndex1 + columnIndex2 + channel ] - 2 * smoothed[ rowIndex2 + columnIndex2 + channel ] + smoothed[ rowIndex3 + columnIndex2 + channel ];

				double const dxdy = -0.25 * smoothed[ rowIndex1 + columnIndex1 + channel ] + 0.25 * smoothed[ rowIndex1 + columnIndex3 + channel ] + 0.25 * smoothed[ rowIndex3 + columnIndex1 + channel ] - 0.25 * smoothed[ rowIndex3 + columnIndex3 + channel ];

				double const squaredNorm = Square( dx ) + Square( dy );

				double numerator = Square( squaredNorm );
				double const denominator = std::abs( Square( dx ) * dxdx + 2 * dx * dy * dxdy + Square( dy ) * dydy ) + lambda * squaredNorm;
				if ( numerator > std::numeric_limits< double >::epsilon() )
					numerator /= denominator;

				result[ rowIndex2 + columnIndex2 + channel ] = numerator;
			}
		}
	}

	return result;
}




//============================================================================
//    FindWarpedImage function
//============================================================================


SharedArray< float > const FindWarpedImage(
	SharedArray< float > const& image,
	unsigned int const rows,
	unsigned int const columns,
	SharedArray< double > const& redCoefficients,
	SharedArray< double > const& greenCoefficients,
	SharedArray< double > const& blueCoefficients,
	unsigned int const degree
)
{
	SharedArray< double > const* coefficients[] = { &redCoefficients, &greenCoefficients, &blueCoefficients };

	double const scale = 0.5 * std::sqrt( Square( static_cast< double >( rows ) ) + Square( static_cast< double >( columns ) ) );

	SharedArray< float > result( rows * columns * 3 );
	#pragma omp parallel for schedule( static )
	for ( int ii = 0; ii < static_cast< int >( rows ); ++ii ) {

		for ( unsigned int jj = 0; jj < columns; ++jj ) {

			double const yy = ( ii + 0.5 ) - 0.5 * rows;
			double const xx = ( jj + 0.5 ) - 0.5 * columns;
			double const radius = std::sqrt( Square( xx ) + Square( yy ) ) / scale;

			for ( unsigned int channel = 0; channel < 3; ++channel ) {

				double scale = 0;
				{	double term = 1;
					for ( unsigned int kk = 0; kk < degree + 1; ++kk, term *= radius )
						scale += ( *coefficients[ channel ] )[ kk ] * term;
				}

				double const row    = yy * scale + 0.5 * rows;
				double const column = xx * scale + 0.5 * columns;

				int rowIndex = std::floor( row - 0.5 );
				if ( rowIndex < 0 )
					rowIndex = 0;
				else if ( rowIndex + 1 >= static_cast< int >( rows ) )
					rowIndex = rows - 2;

				int columnIndex = std::floor( column - 0.5 );
				if ( columnIndex < 0 )
					columnIndex = 0;
				else if ( columnIndex + 1 >= static_cast< int >( columns ) )
					columnIndex = columns - 2;

				double const rowLambda    = rowIndex    + 1.5 - row;
				double const columnLambda = columnIndex + 1.5 - column;

				unsigned int const index = rowIndex * columns * 3 + columnIndex * 3;

				double value = 0;
				value += image[ index                   + channel ] *       rowLambda   *       columnLambda;
				value += image[ index               + 3 + channel ] *       rowLambda   * ( 1 - columnLambda );
				value += image[ index + columns * 3     + channel ] * ( 1 - rowLambda ) *       columnLambda;
				value += image[ index + columns * 3 + 3 + channel ] * ( 1 - rowLambda ) * ( 1 - columnLambda );

				result[ ( ii * columns + jj ) * 3 + channel ] = value;
			}
		}
	}

	return result;
}




}    // anomymous namespace




//============================================================================
//    FixChromaticAberration function
//============================================================================


SharedArray< float > const FixChromaticAberration(
	SharedArray< float > const& image,
	unsigned int const rows,
	unsigned int const columns,
	int const degree,
	double const sigma1,
	double const sigma2,
	double const lambda,
	double const alpha,
	double const gamma,
	double const rho,
	double const sigma,
	double const delta,
	double const epsilon
)
{
	if ( degree < 0 )
		throw std::runtime_error( "degree must be positive" );

	SharedArray< float > edges = FindEdges( image, rows, columns, sigma1, sigma2, lambda );

	SharedArray< double > redCoefficients(  degree + 1 );
	SharedArray< double > blueCoefficients( degree + 1 );

	SharedArray< double > greenCoefficients( degree + 1 );
	std::fill( &greenCoefficients[ 0 ], &greenCoefficients[ 0 ] + degree + 1, 0 );
	greenCoefficients[ 0 ] = 1;

	// Nelder-Mead simplex algorithm
	{	unsigned int const size = ( degree + 1 ) * 2;

		SharedArray< double > coefficients( ( size + 1 ) * size );
		std::fill( &coefficients[ 0 ], &coefficients[ 0 ] + ( size + 1 ) * size, 0 );
		for ( unsigned int ii = 0; ii < ( size + 1 ); ++ii ) {

			coefficients[ ii * size                  ] = 1;
			coefficients[ ii * size + ( degree + 1 ) ] = 1;
		}
		for ( unsigned int ii = 0; ii < ( size + 1 ); ++ii ) {

			double* const pCoefficients = &coefficients[ 0 ] + ii * size;
			for ( unsigned int kk = 0; kk < ii; ++kk )
				pCoefficients[ kk ] -= delta;
			if ( ii < size )
				pCoefficients[ ii ] += delta;
		}

		SharedArray< double > scores( size + 1 );
		for ( unsigned int ii = 0; ii < ( size + 1 ); ++ii ) {

			double* const pCoefficients = &coefficients[ 0 ] + ii * size;
			std::copy( pCoefficients,                  pCoefficients + ( degree + 1 ), &redCoefficients[  0 ] );
			std::copy( pCoefficients + ( degree + 1 ), pCoefficients + size,           &blueCoefficients[ 0 ] );

			SharedArray< float > warped = FindWarpedImage(
				edges,
				rows,
				columns,
				redCoefficients,
				greenCoefficients,
				blueCoefficients,
				degree
			);

			double score = 0;
			#pragma omp parallel for schedule( static ) reduction( + : score )
			for ( int jj = 0; jj < static_cast< int >( rows * columns ); ++jj )
				score += Square( warped[ jj * 3 + 0 ] - warped[ jj * 3 + 1 ] ) + Square( warped[ jj * 3 + 2 ] - warped[ jj * 3 + 1 ] );

			scores[ ii ] = score;
		}

		SharedArray< double > centroid(   size );
		SharedArray< double > reflected(  size );
		SharedArray< double > expanded(   size );
		SharedArray< double > contracted( size );

		for ( ; ; ) {

			unsigned int minimumIndex = 0;
			unsigned int maximumIndex = 0;
			for ( unsigned int ii = 0; ii < ( size + 1 ); ++ii ) {

				if ( scores[ ii ] < scores[ minimumIndex ] )
					minimumIndex = ii;
				if ( scores[ ii ] > scores[ maximumIndex ] )
					maximumIndex = ii;
			}

			unsigned int submaximumIndex = minimumIndex;
			for ( unsigned int ii = 0; ii < ( size + 1 ); ++ii ) {

				if ( ( ii != maximumIndex ) && ( scores[ ii ] > scores[ submaximumIndex ] ) )
					submaximumIndex = ii;
			}

			double maximumError = 0;
			for ( unsigned int ii = 0; ii < ( size + 1 ); ++ii ) {

				if ( ii != minimumIndex ) {

					double error = 0;
					{	double const* pCoefficients1 = &coefficients[ 0 ] + ii           * size;
						double const* pCoefficients2 = &coefficients[ 0 ] + minimumIndex * size;
						for ( unsigned int jj = 0; jj < size; ++jj )
							error += Square( pCoefficients1[ jj ] - pCoefficients2[ jj ] );
					}

					if ( error > maximumError )
						maximumError = error;
				}
			}
			maximumError = std::sqrt( maximumError );
			if ( maximumError < epsilon )
				break;

			std::fill( &centroid[ 0 ], &centroid[ 0 ] + size, 0 );
			for ( unsigned int ii = 0; ii < ( size + 1 ); ++ii ) {

				if ( ii != maximumIndex ) {

					double const* const pCoefficients = &coefficients[ 0 ] + ii * size;
					for ( unsigned int jj = 0; jj < size; ++jj )
						centroid[ jj ] += pCoefficients[ jj ];
				}
			}
			for ( unsigned int ii = 0; ii < size; ++ii )
				centroid[ ii ] /= size;

			double reflectedScore = 0;
			{	double const* const pCoefficients = &coefficients[ 0 ] + maximumIndex * size;
				for ( unsigned int ii = 0; ii < size; ++ii )
					reflected[ ii ] = centroid[ ii ] + alpha * ( centroid[ ii ] - pCoefficients[ ii ] );

				std::copy( &reflected[ 0 ],                  &reflected[ 0 ] + ( degree + 1 ), &redCoefficients[  0 ] );
				std::copy( &reflected[ 0 ] + ( degree + 1 ), &reflected[ 0 ] + size,           &blueCoefficients[ 0 ] );

				SharedArray< float > warped = FindWarpedImage(
					edges,
					rows,
					columns,
					redCoefficients,
					greenCoefficients,
					blueCoefficients,
					degree
				);

				#pragma omp parallel for schedule( static ) reduction( + : reflectedScore )
				for ( int jj = 0; jj < static_cast< int >( rows * columns ); ++jj )
					reflectedScore += Square( warped[ jj * 3 + 0 ] - warped[ jj * 3 + 1 ] ) + Square( warped[ jj * 3 + 2 ] - warped[ jj * 3 + 1 ] );
			}

			if ( reflectedScore < scores[ minimumIndex ] ) {    // expansion

				double expandedScore = 0;
				{	double const* const pCoefficients = &coefficients[ 0 ] + maximumIndex * size;
					for ( unsigned int ii = 0; ii < size; ++ii )
						expanded[ ii ] = centroid[ ii ] + gamma * ( centroid[ ii ] - pCoefficients[ ii ] );

					std::copy( &expanded[ 0 ],                  &expanded[ 0 ] + ( degree + 1 ), &redCoefficients[  0 ] );
					std::copy( &expanded[ 0 ] + ( degree + 1 ), &expanded[ 0 ] + size,           &blueCoefficients[ 0 ] );

					SharedArray< float > warped = FindWarpedImage(
						edges,
						rows,
						columns,
						redCoefficients,
						greenCoefficients,
						blueCoefficients,
						degree
					);

					#pragma omp parallel for schedule( static ) reduction( + : expandedScore )
					for ( int jj = 0; jj < static_cast< int >( rows * columns ); ++jj )
						expandedScore += Square( warped[ jj * 3 + 0 ] - warped[ jj * 3 + 1 ] ) + Square( warped[ jj * 3 + 2 ] - warped[ jj * 3 + 1 ] );
				}

				if ( expandedScore <= reflectedScore ) {

					std::copy( &expanded[ 0 ], &expanded[ 0 ] + size, &coefficients[ 0 ] + maximumIndex * size );
					scores[ maximumIndex ] = expandedScore;
				}
				else {

					std::copy( &reflected[ 0 ], &reflected[ 0 ] + size, &coefficients[ 0 ] + maximumIndex * size );
					scores[ maximumIndex ] = reflectedScore;
				}
			}
			else if ( reflectedScore < scores[ submaximumIndex ] ) {    // reflection

				std::copy( &reflected[ 0 ], &reflected[ 0 ] + size, &coefficients[ 0 ] + maximumIndex * size );
				scores[ maximumIndex ] = reflectedScore;
			}
			else {    // contraction

				double contractedScore = 0;
				{	double const* const pCoefficients = &coefficients[ 0 ] + maximumIndex * size;
					for ( unsigned int ii = 0; ii < size; ++ii )
						contracted[ ii ] = pCoefficients[ ii ] + rho * ( centroid[ ii ] - pCoefficients[ ii ] );

					std::copy( &contracted[ 0 ],                  &contracted[ 0 ] + ( degree + 1 ), &redCoefficients[  0 ] );
					std::copy( &contracted[ 0 ] + ( degree + 1 ), &contracted[ 0 ] + size,           &blueCoefficients[ 0 ] );

					SharedArray< float > warped = FindWarpedImage(
						edges,
						rows,
						columns,
						redCoefficients,
						greenCoefficients,
						blueCoefficients,
						degree
					);

					#pragma omp parallel for schedule( static ) reduction( + : contractedScore )
					for ( int jj = 0; jj < static_cast< int >( rows * columns ); ++jj )
						contractedScore += Square( warped[ jj * 3 + 0 ] - warped[ jj * 3 + 1 ] ) + Square( warped[ jj * 3 + 2 ] - warped[ jj * 3 + 1 ] );
				}

				if ( contractedScore < scores[ maximumIndex ] ) {

					std::copy( &contracted[ 0 ], &contracted[ 0 ] + size, &coefficients[ 0 ] + maximumIndex * size );
					scores[ maximumIndex ] = contractedScore;
				}
				else {    // reduction

					for ( unsigned int ii = 0; ii < ( size + 1 ); ++ii ) {

						if ( ii != minimumIndex ) {

							double*       pCoefficients1 = &coefficients[ 0 ] + ii           * size;
							double const* pCoefficients2 = &coefficients[ 0 ] + minimumIndex * size;
							for ( unsigned int jj = 0; jj < size; ++jj )
								pCoefficients1[ jj ] = pCoefficients2[ jj ] + sigma * ( pCoefficients1[ jj ] - pCoefficients2[ jj ] );

							std::copy( pCoefficients1,                  pCoefficients1 + ( degree + 1 ), &redCoefficients[  0 ] );
							std::copy( pCoefficients1 + ( degree + 1 ), pCoefficients1 + size,           &blueCoefficients[ 0 ] );

							SharedArray< float > warped = FindWarpedImage(
								edges,
								rows,
								columns,
								redCoefficients,
								greenCoefficients,
								blueCoefficients,
								degree
							);

							double score = 0;
							#pragma omp parallel for schedule( static ) reduction( + : score )
							for ( int jj = 0; jj < static_cast< int >( rows * columns ); ++jj )
								score += Square( warped[ jj * 3 + 0 ] - warped[ jj * 3 + 1 ] ) + Square( warped[ jj * 3 + 2 ] - warped[ jj * 3 + 1 ] );

							scores[ ii ] = score;
						}
					}
				}
			}
		}

		{	unsigned int minimumIndex = 0;
			for ( unsigned int ii = 0; ii < ( size + 1 ); ++ii ) {

				if ( scores[ ii ] < scores[ minimumIndex ] )
					minimumIndex = ii;
			}

			double* const pCoefficients = &coefficients[ 0 ] + minimumIndex * size;
			std::copy( pCoefficients,                  pCoefficients + ( degree + 1 ), &redCoefficients[  0 ] );
			std::copy( pCoefficients + ( degree + 1 ), pCoefficients + size,           &blueCoefficients[ 0 ] );
		}
	}

	return FindWarpedImage(
		image,
		rows,
		columns,
		redCoefficients,
		greenCoefficients,
		blueCoefficients,
		degree
	);
}
