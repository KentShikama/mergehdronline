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
	\file helpers.hh
	\brief Helper macros, functions and classes
*/




#ifndef __HELPERS_HH__
#define __HELPERS_HH__




#include <limits>
#include <cmath>

#include <stdint.h>




//============================================================================
//    ARRAYLENGTH macro
//============================================================================


#define ARRAYLENGTH( array )  \
	( sizeof( array ) / sizeof( array[ 0 ] ) )




//============================================================================
//    LIKELY and UNLIKELY macros
//============================================================================


#if defined( __GNUC__ ) && ( __GNUC__ >= 3 )

#define LIKELY( boolean ) __builtin_expect( ( boolean ), 1 )
#define UNLIKELY( boolean ) __builtin_expect( ( boolean ), 0 )

#else    /* defined( __GNUC__ ) && ( __GNUC__ >= 3 ) */

#define LIKELY( boolean ) ( boolean )
#define UNLIKELY( boolean ) ( boolean )

#endif    /* defined( __GNUC__ ) && ( __GNUC__ >= 3 ) */




#ifdef __cplusplus




//============================================================================
//    Power helper template
//============================================================================


template< unsigned int t_Number, unsigned int t_Power >
struct Power {

	enum { RESULT = t_Number * Power< t_Number, t_Power - 1 >::RESULT };
};


template< unsigned int t_Number >
struct Power< t_Number, 0 > {

	enum { RESULT = 1 };
};




//============================================================================
//    Signum helper functions
//============================================================================


template< typename t_Type >
inline t_Type const Signum( t_Type const& value ) {

	t_Type result = 0;
	if ( value < 0 )
		result = -1;
	else if ( value > 0 )
		result = 1;
	return result;
}


template< typename t_Type >
inline t_Type const Signum( t_Type const& value, t_Type const& scale ) {

	t_Type result = 0;
	if ( value < 0 )
		result = -scale;
	else if ( value > 0 )
		result = scale;
	return result;
}




//============================================================================
//    Square helper function
//============================================================================


template< typename t_Type >
inline t_Type const Square( t_Type const& value ) {

	return( value * value );
}




//============================================================================
//    Cube helper function
//============================================================================


template< typename t_Type >
inline t_Type const Cube( t_Type const& value ) {

	return( value * value * value );
}




//============================================================================
//    Pow helper function
//============================================================================


template< typename t_Type >
inline t_Type const Pow( t_Type const& value, int const power ) {

	if ( power < 0 )
		return( 1 / Pow( value, -power ) );

	t_Type result( 1 );
	if ( power > 0 ) {

		t_Type accumulator( value );
		unsigned int bits = power;
		for ( unsigned int bit = 1; bits != 0; bit += bit, accumulator *= accumulator ) {

			if ( ( bits & bit ) != 0 ) {

				result *= accumulator;
				bits &= ~bit;
			}
		}
	}
	return result;
}




//============================================================================
//    FastExp function
//============================================================================


inline double FastExp( double exponent ) {

	static double const scale = ( 1ull << 52 ) / logl( 2.0 );
	static int64_t const shift = roundl( ( 1ull << 52 ) * ( 0x03ff - log2l( 0.5 + ( 1.0 / ( expl( 1.0 ) * logl( 2.0 ) ) ) ) ) );

	// 709 = floor( 2^10 * log( 2 ) )
	if ( exponent < -709 )
		return 0;
	else if ( exponent > 709 )
		return std::numeric_limits< double >::infinity();
	else {

		union { uint64_t ii; double ff; } uu;
		uu.ii = exponent * scale + 0.5;
		uu.ii += shift;
		return uu.ff;
	}
}




#endif    /* __cplusplus */




#endif    /* __HELPERS_HH__ */
