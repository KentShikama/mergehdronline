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
	\file solve_poisson.cc
	\brief implementation of SolvePoisson function
*/




#include "headers.hh"




namespace {




//============================================================================
//    AddMatrixVectorProduct function
//============================================================================


// adds scalar * L * rightVector to result, where L is the Laplacian operator
void AddMatrixVectorProduct(
	SharedArray< float > const& result,
	SharedArray< float > const& rightVector,
	unsigned int const rows,
	unsigned int const columns,
	float const scalar = 1
)
{
	{	float const* pRightVector = &rightVector[ 0 ];
		float*       pResult      = &result[ 0 ];
		*pResult -= scalar * ( 2.0 * *pRightVector - *( pRightVector + 1 ) - *( pRightVector + columns ) );
		++pRightVector;
		++pResult;
		for ( unsigned int jj = 1; jj < columns - 1; ++jj, ++pRightVector, ++pResult )
			*pResult -= scalar * ( 3.0 * *pRightVector - *( pRightVector - 1 ) - *( pRightVector + 1 ) - *( pRightVector + columns ) );
		*pResult -= scalar * ( 2.0 * *pRightVector - *( pRightVector - 1 ) - *( pRightVector + columns ) );
	}
	#pragma omp parallel for schedule( static )
	for ( int ii = 1; ii < static_cast< int >( rows ) - 1; ++ii ) {

		float const* pRightVector = &rightVector[ 0 ] + ii * columns;
		float*       pResult      = &result[ 0 ]      + ii * columns;
		*pResult -= scalar * ( 3.0 * *pRightVector - *( pRightVector + 1 ) - *( pRightVector - columns ) - *( pRightVector + columns ) );
		++pRightVector;
		++pResult;
		for ( unsigned int jj = 1; jj < columns - 1; ++jj, ++pRightVector, ++pResult )
			*pResult -= scalar * ( 4.0 * *pRightVector - *( pRightVector - 1 ) - *( pRightVector + 1 ) - *( pRightVector - columns ) - *( pRightVector + columns ) );
		*pResult -= scalar * ( 3.0 * *pRightVector - *( pRightVector - 1 ) - *( pRightVector - columns ) - *( pRightVector + columns ) );
	}
	{	float const* pRightVector = &rightVector[ 0 ] + ( rows - 1 ) * columns;
		float*       pResult      = &result[ 0 ]      + ( rows - 1 ) * columns;
		*pResult -= scalar * ( 2.0 * *pRightVector - *( pRightVector + 1 ) - *( pRightVector - columns ) );
		++pRightVector;
		++pResult;
		for ( unsigned int jj = 1; jj < columns - 1; ++jj, ++pRightVector, ++pResult )
			*pResult -= scalar * ( 3.0 * *pRightVector - *( pRightVector - 1 ) - *( pRightVector + 1 ) - *( pRightVector - columns ) );
		*pResult -= scalar * ( 2.0 * *pRightVector - *( pRightVector - 1 ) - *( pRightVector - columns ) );
	}
}




//============================================================================
//    VectorMatrixVectorProduct function
//============================================================================


// calculates leftVector * L * rightVector, where L is the Laplacian operator
double const VectorMatrixVectorProduct(
	SharedArray< float > const& leftVector,
	SharedArray< float > const& rightVector,
	unsigned int const rows,
	unsigned int const columns
)
{
	double result     = 0;
	double edgeResult = 0;

	{	float const* pLeftVector  = &leftVector[ 0 ];
		float const* pRightVector = &rightVector[ 0 ];
		edgeResult -= *pLeftVector * ( 2.0 * *pRightVector - *( pRightVector + 1 ) - *( pRightVector + columns ) );
		++pLeftVector;
		++pRightVector;
		for ( unsigned int jj = 1; jj < columns - 1; ++jj, ++pLeftVector, ++pRightVector )
			edgeResult -= *pLeftVector * ( 3.0 * *pRightVector - *( pRightVector - 1 ) - *( pRightVector + 1 ) - *( pRightVector + columns ) );
		edgeResult -= *pLeftVector * ( 2.0 * *pRightVector - *( pRightVector - 1 ) - *( pRightVector + columns ) );
	}
	#pragma omp parallel for schedule( static ) reduction( + : result )
	for ( int ii = 1; ii < static_cast< int >( rows ) - 1; ++ii ) {

		float const* pLeftVector  = &leftVector[ 0 ]  + ii * columns;
		float const* pRightVector = &rightVector[ 0 ] + ii * columns;
		result -= *pLeftVector * ( 3.0 * *pRightVector - *( pRightVector + 1 ) - *( pRightVector - columns ) - *( pRightVector + columns ) );
		++pLeftVector;
		++pRightVector;
		for ( unsigned int jj = 1; jj < columns - 1; ++jj, ++pLeftVector, ++pRightVector )
			result -= *pLeftVector * ( 4.0 * *pRightVector - *( pRightVector - 1 ) - *( pRightVector + 1 ) - *( pRightVector - columns ) - *( pRightVector + columns ) );
		result -= *pLeftVector * ( 3.0 * *pRightVector - *( pRightVector - 1 ) - *( pRightVector - columns ) - *( pRightVector + columns ) );
	}
	{	float const* pLeftVector  = &leftVector[ 0 ]  + ( rows - 1 ) * columns;
		float const* pRightVector = &rightVector[ 0 ] + ( rows - 1 ) * columns;
		edgeResult -= *pLeftVector * ( 2.0 * *pRightVector - *( pRightVector + 1 ) - *( pRightVector - columns ) );
		++pLeftVector;
		++pRightVector;
		for ( unsigned int jj = 1; jj < columns - 1; ++jj, ++pLeftVector, ++pRightVector )
			edgeResult -= *pLeftVector * ( 3.0 * *pRightVector - *( pRightVector - 1 ) - *( pRightVector + 1 ) - *( pRightVector - columns ) );
		edgeResult -= *pLeftVector * ( 2.0 * *pRightVector - *( pRightVector - 1 ) - *( pRightVector - columns ) );
	}

	return( result + edgeResult );
}




//============================================================================
//    ConjugateGradient function
//============================================================================


// performs the specified number of conjugate gradient iterations on the problem L * solution = rhs, where L is the Laplacian operator
void ConjugateGradient(
	SharedArray< float > const& solution,
	SharedArray< float > const& rhs,
	unsigned int const rows,
	unsigned int const columns,
	unsigned int const iterations
)
{
	SharedArray< float > gradient( rows * columns );
	SharedArray< float > step(     rows * columns );

	std::fill( &gradient[ 0 ], &gradient[ 0 ] + rows * columns, 0 );
	std::fill( &step[     0 ], &step[     0 ] + rows * columns, 0 );

	double gradientNormSquared = 0;
	for ( unsigned int ii = 0; ii < iterations; ++ii ) {

		double const oldGradientNormSquared = gradientNormSquared;

		// calculate the gradient
		std::copy( &rhs[ 0 ], &rhs[ 0 ] + rows * columns, &gradient[ 0 ] );
		AddMatrixVectorProduct( gradient, solution, rows, columns, -1 );

		// center the gradient, and find its squared norm
		gradientNormSquared = 0;
		{	double gradientAverage = 0;
			#pragma omp parallel for schedule( static ) reduction( + : gradientAverage )
			for ( int jj = 0; jj < static_cast< int >( rows * columns ); ++jj )
				gradientAverage += gradient[ jj ];
			gradientAverage /= rows * columns;
			#pragma omp parallel for schedule( static ) reduction( + : gradientNormSquared )
			for ( int jj = 0; jj < static_cast< int >( rows * columns ); ++jj ) {

				gradient[ jj ] -= gradientAverage;
				gradientNormSquared += Square( gradient[ jj ] );
			}
		}

		// update the step
		double beta = 0;
		if ( std::fabs( oldGradientNormSquared ) > std::numeric_limits< double >::epsilon() )
			beta = gradientNormSquared / oldGradientNormSquared;
		#pragma omp parallel for schedule( static )
		for ( int jj = 0; jj < static_cast< int >( rows * columns ); ++jj )
			step[ jj ] = beta * step[ jj ] + gradient[ jj ];

		// find the step size
		double const alphaDenominator = VectorMatrixVectorProduct( step, step, rows, columns );
		if ( std::fabs( alphaDenominator ) < std::numeric_limits< double >::epsilon() )
			break;
		double const alpha = gradientNormSquared / alphaDenominator;

		// take a step
		#pragma omp parallel for schedule( static )
		for ( int jj = 0; jj < static_cast< int >( rows * columns ); ++jj )
			solution[ jj ] += alpha * step[ jj ];
	}
}




//============================================================================
//    Restrict function
//============================================================================


void Restrict(
	SharedArray< float > const& restriction,
	SharedArray< float > const& rhs,
	unsigned int const rows,
	unsigned int const columns
)
{
	// scale RHS by a factor of 4, since we don't scale the Laplacian operator

	unsigned int const restrictionRows    = ( rows    + 1 ) / 2;
	unsigned int const restrictionColumns = ( columns + 1 ) / 2;

	{	unsigned int const restrictionRowsEnd    = rows    / 2;
		unsigned int const restrictionColumnsEnd = columns / 2;

		{	float const* pRHS         = &rhs[ 0 ];
			float*       pRestriction = &restriction[ 0 ];
			*pRestriction = (
				*pRHS +
				0.5 * (
					*( pRHS + 1 ) +
					*( pRHS + columns )
				) +
				0.25 * (
					*( pRHS + 1 + columns )
				)
			);
			++pRestriction;
			pRHS += 2;
			for ( unsigned int jj = 1; jj < restrictionColumnsEnd; ++jj, ++pRestriction, pRHS += 2 ) {

				*pRestriction = (
					*pRHS +
					0.5 * (
						*( pRHS - 1 ) +
						*( pRHS + 1 ) +
						*( pRHS + columns )
					) +
					0.25 * (
						*( pRHS - 1 + columns ) +
						*( pRHS + 1 + columns )
					)
				);
			}
			if ( restrictionColumnsEnd < restrictionColumns ) {

				*pRestriction = (
					*pRHS +
					0.5 * (
						*( pRHS - 1 ) +
						*( pRHS + columns )
					) +
					0.25 * (
						*( pRHS - 1 + columns )
					)
				);
			}
		}
		#pragma omp parallel for schedule( static )
		for ( int ii = 1; ii < static_cast< int >( restrictionRowsEnd ); ++ii ) {

			float const* pRHS         = &rhs[ 0 ]         + ( ii * 2 ) * columns;
			float*       pRestriction = &restriction[ 0 ] +   ii       * restrictionColumns;
			*pRestriction = (
				*pRHS +
				0.5 * (
					*( pRHS + 1 ) +
					*( pRHS - columns ) +
					*( pRHS + columns )
				) +
				0.25 * (
					*( pRHS + 1 - columns ) +
					*( pRHS + 1 + columns )
				)
			);
			++pRestriction;
			pRHS += 2;
			for ( unsigned int jj = 1; jj < restrictionColumnsEnd; ++jj, ++pRestriction, pRHS += 2 ) {

				*pRestriction = (
					*pRHS +
					0.5 * (
						*( pRHS - 1 ) +
						*( pRHS + 1 ) +
						*( pRHS - columns ) +
						*( pRHS + columns )
					) +
					0.25 * (
						*( pRHS - 1 - columns ) +
						*( pRHS + 1 - columns ) +
						*( pRHS - 1 + columns ) +
						*( pRHS + 1 + columns )
					)
				);
			}
			if ( restrictionColumnsEnd < restrictionColumns ) {

				*pRestriction = (
					*pRHS +
					0.5 * (
						*( pRHS - 1 ) +
						*( pRHS - columns ) +
						*( pRHS + columns )
					) +
					0.25 * (
						*( pRHS - 1 - columns ) +
						*( pRHS - 1 + columns )
					)
				);
			}
		}
		if ( restrictionRowsEnd < restrictionRows ) {

			float const* pRHS         = &rhs[ 0 ]         + ( ( restrictionRows - 1 ) * 2 ) * columns;
			float*       pRestriction = &restriction[ 0 ] +   ( restrictionRows - 1 )       * restrictionColumns;
			*pRestriction = (
				*pRHS +
				0.5 * (
					*( pRHS + 1 ) +
					*( pRHS - columns )
				) +
				0.25 * (
					*( pRHS + 1 - columns )
				)
			);
			++pRestriction;
			pRHS += 2;
			for ( unsigned int jj = 1; jj < restrictionColumnsEnd; ++jj, ++pRestriction, pRHS += 2 ) {

				*pRestriction = (
					*pRHS +
					0.5 * (
						*( pRHS - 1 ) +
						*( pRHS + 1 ) +
						*( pRHS - columns )
					) +
					0.25 * (
						*( pRHS - 1 - columns ) +
						*( pRHS + 1 - columns )
					)
				);
			}
			if ( restrictionColumnsEnd < restrictionColumns ) {

				*pRestriction = (
					*pRHS +
					0.5 * (
						*( pRHS - 1 ) +
						*( pRHS - columns )
					) +
					0.25 * (
						*( pRHS - 1 - columns )
					)
				);
			}
		}
	}
}




//============================================================================
//    Prolong function
//============================================================================


void Prolong(
	SharedArray< float > const& solution,
	SharedArray< float > const& restriction,
	unsigned int const rows,
	unsigned int const columns
)
{
	unsigned int const restrictionRows    = ( rows    + 1 ) / 2;
	unsigned int const restrictionColumns = ( columns + 1 ) / 2;

	#pragma omp parallel for schedule( static )
	for ( int ii = 0; ii < static_cast< int >( rows ); ++ii ) {

		unsigned int iiIndex = ii / 2;
		if ( iiIndex + 1 >= restrictionRows )
			iiIndex = restrictionRows - 2;
		assert( iiIndex + 1 < restrictionRows );

		double const rowLambda = 0.5 * ( ( iiIndex + 1 ) * 2.0 - ii );

		float const* const pSource      = &restriction[ 0 ] + iiIndex * restrictionColumns;
		float*             pDestination = &solution[ 0 ]    + ii      * columns;
		for ( unsigned int jj = 0; jj < columns; ++jj, ++pDestination ) {

			unsigned int jjIndex = jj / 2;
			if ( jjIndex + 1 >= restrictionColumns )
				jjIndex = restrictionColumns - 2;
			assert( jjIndex + 1 < restrictionColumns );

			double const columnLambda = 0.5 * ( ( jjIndex + 1 ) * 2.0 - jj );

			double value = 0;
			value += *( pSource + jjIndex                          ) *       rowLambda   *       columnLambda;
			value += *( pSource + jjIndex                      + 1 ) *       rowLambda   * ( 1 - columnLambda );
			value += *( pSource + jjIndex + restrictionColumns     ) * ( 1 - rowLambda ) *       columnLambda;
			value += *( pSource + jjIndex + restrictionColumns + 1 ) * ( 1 - rowLambda ) * ( 1 - columnLambda );

			*pDestination += value;
		}
	}
}




}    // anomymous namespace




//============================================================================
//    SolvePoisson function
//============================================================================


/*
	Optimizes via the full multigrid algorithm (tweaked to make it adaptive),
	using conjugate gradient smoothing iterations (which work a lot better than
	Gauss-Seidel).

	Press, Teukolsky, Vetterling and Flannery. "Numerical recipies in C++,
	Second Edition". Pages 873-884.
*/
void SolvePoisson(
	SharedArray< float > const& solution,
	SharedArray< float > const& rhs,
	unsigned int const rows,
	unsigned int const columns,
	double const epsilon,    // termination criterion is this proportion of the L2 norm of rhs
	unsigned int const innerIterations
)
{
	if ( epsilon <= 0 )
		throw std::runtime_error( "epsilon must be positive" );

	// center the solution
	{	double average = 0;
		#pragma omp parallel for schedule( static ) reduction( + : average )
		for ( int ii = 0; ii < static_cast< int >( rows * columns ); ++ii )
			average += solution[ ii ];
		average /= rows * columns;
		#pragma omp parallel for schedule( static )
		for ( int ii = 0; ii < static_cast< int >( rows * columns ); ++ii )
			solution[ ii ] -= average;
	}

	// calculate the norm of the RHS
	double rhsNorm = 0;
	#pragma omp parallel for schedule( static ) reduction( + : rhsNorm )
	for ( int ii = 0; ii < static_cast< int >( rows * columns ); ++ii )
		rhsNorm += Square( rhs[ ii ] );
	rhsNorm = std::sqrt( rhsNorm / ( rows * columns ) );

	std::vector< SharedArray< float > > pyramidSolutions;
	std::vector< SharedArray< float > > pyramidRHSs;
	std::vector< std::pair< unsigned int, unsigned int > > pyramidSizes;

	pyramidSolutions.push_back( solution );
	pyramidRHSs.push_back( rhs );
	pyramidSizes.push_back( std::pair< unsigned int, unsigned int >( rows, columns ) );

	for ( ; ; ) {

		unsigned int const sourceRows    = pyramidSizes.back().first;
		unsigned int const sourceColumns = pyramidSizes.back().second;

		unsigned int const destinationRows    = ( sourceRows    + 1 ) / 2;
		unsigned int const destinationColumns = ( sourceColumns + 1 ) / 2;
		if ( ( destinationRows < 2 ) || ( destinationColumns < 2 ) )
			break;

		SharedArray< float > destinationSolution( destinationRows * destinationColumns );
		SharedArray< float > destinationRHS(      destinationRows * destinationColumns );

		pyramidSolutions.push_back( destinationSolution );
		pyramidRHSs.push_back( destinationRHS );
		pyramidSizes.push_back( std::pair< unsigned int, unsigned int >( destinationRows, destinationColumns ) );
	}

	unsigned int const size = pyramidSolutions.size();
	assert( size == pyramidRHSs.size()  );
	assert( size == pyramidSizes.size() );

	SharedArray< float > residual( rows * columns );

	unsigned int depth = 0;

	for ( ; ; ) {

		SharedArray< float > const& pyramidSolution = pyramidSolutions[ depth ];
		SharedArray< float > const& pyramidRHS      = pyramidRHSs[      depth ];
		unsigned int const pyramidRows    = pyramidSizes[ depth ].first;
		unsigned int const pyramidColumns = pyramidSizes[ depth ].second;

		ConjugateGradient( pyramidSolution, pyramidRHS, pyramidRows, pyramidColumns, innerIterations );

		// find the residual
		std::copy( &pyramidRHS[ 0 ], &pyramidRHS[ 0 ] + pyramidRows * pyramidColumns, &residual[ 0 ] );
		AddMatrixVectorProduct( residual, pyramidSolution, pyramidRows, pyramidColumns, -1 );

		// center the residual, and find its norm
		double residualNorm = 0;
		{	double average = 0;
			#pragma omp parallel for schedule( static )reduction( + : average )
			for ( int ii = 0; ii < static_cast< int >( pyramidRows * pyramidColumns ); ++ii )
				average += residual[ ii ];
			average /= pyramidRows * pyramidColumns;
			#pragma omp parallel for schedule( static ) reduction( + : residualNorm )
			for ( int ii = 0; ii < static_cast< int >( pyramidRows * pyramidColumns ); ++ii ) {

				residual[ ii ] -= average;
				residualNorm += Square( residual[ ii ] );
			}
			residualNorm = std::sqrt( residualNorm / ( pyramidRows * pyramidColumns ) );
		}

		if ( ( residualNorm <= epsilon * rhsNorm ) || ( residualNorm < std::numeric_limits< double >::epsilon() ) ) {

			if ( depth > 0 ) {

				SharedArray< float > const& nextPyramidSolution = pyramidSolutions[ depth - 1 ];
				unsigned int const nextPyramidRows    = pyramidSizes[ depth - 1 ].first;
				unsigned int const nextPyramidColumns = pyramidSizes[ depth - 1 ].second;

				Prolong( nextPyramidSolution, pyramidSolution, nextPyramidRows, nextPyramidColumns );

				// center the prolonged solution
				{	double average = 0;
					#pragma omp parallel for schedule( static ) reduction( + : average )
					for ( int ii = 0; ii < static_cast< int >( nextPyramidRows * nextPyramidColumns ); ++ii )
						average += nextPyramidSolution[ ii ];
					average /= nextPyramidRows * nextPyramidColumns;
					#pragma omp parallel for schedule( static )
					for ( int ii = 0; ii < static_cast< int >( nextPyramidRows * nextPyramidColumns ); ++ii )
						nextPyramidSolution[ ii ] -= average;
				}

				--depth;
			}
			else
				break;
		}
		else {

			if ( depth + 1 < size ) {

				SharedArray< float > const& nextPyramidSolution = pyramidSolutions[ depth + 1 ];
				SharedArray< float > const& nextPyramidRHS      = pyramidRHSs[      depth + 1 ];
				unsigned int const nextPyramidRows    = pyramidSizes[ depth + 1 ].first;
				unsigned int const nextPyramidColumns = pyramidSizes[ depth + 1 ].second;

				std::fill( &nextPyramidSolution[ 0 ], &nextPyramidSolution[ 0 ] + nextPyramidRows * nextPyramidColumns, 0 );
				Restrict( nextPyramidRHS, residual, pyramidRows, pyramidColumns );

				++depth;
			}
		}
	}
}
