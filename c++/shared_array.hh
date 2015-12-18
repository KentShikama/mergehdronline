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
	\file shared_array.hh
	\brief definition of SharedArray class
*/




#ifndef __SHARED_ARRAY_HH__
#define __SHARED_ARRAY_HH__




#include <assert.h>
#include <stdlib.h>




//============================================================================
//    SharedArray class
//============================================================================


template< typename t_Type >
struct SharedArray {

	inline SharedArray();
	inline explicit SharedArray( unsigned int const size );
	inline SharedArray( SharedArray const& other );

	inline ~SharedArray();


	inline SharedArray const& operator=( SharedArray const& other );


	inline t_Type& operator[]( unsigned int const index ) const;


	inline bool const operator==( t_Type const* const other ) const;
	inline bool const operator==( SharedArray const& other ) const;

	inline bool const operator!=( t_Type const* const other ) const;
	inline bool const operator!=( SharedArray const& other ) const;


private:

	unsigned int mutable* m_pRefcount;
	t_Type* m_array;
};




//============================================================================
//    SharedArray inline methods
//============================================================================


template< typename t_Type >
SharedArray< t_Type >::SharedArray() : m_pRefcount( NULL ), m_array( NULL ) {
}


template< typename t_Type >
SharedArray< t_Type >::SharedArray( unsigned int const size ) : m_pRefcount( new unsigned int ), m_array( new t_Type[ size ] ) {

	*m_pRefcount = 1;
}


template< typename t_Type >
SharedArray< t_Type >::SharedArray( SharedArray const& other ) : m_pRefcount( other.m_pRefcount ), m_array( other.m_array ) {

	assert( ( m_pRefcount == NULL ) == ( m_array == NULL ) );
	if ( m_pRefcount != NULL )
		++( *m_pRefcount );
}


template< typename t_Type >
SharedArray< t_Type >::~SharedArray() {

	assert( ( m_pRefcount == NULL ) == ( m_array == NULL ) );
	if ( m_pRefcount != NULL ) {

		assert( *m_pRefcount > 0 );
		if ( --( *m_pRefcount ) == 0 ) {

			delete m_pRefcount;
			delete[] m_array;
		}
	}
}


template< typename t_Type >
SharedArray< t_Type > const& SharedArray< t_Type >::operator=( SharedArray const& other ) {

	assert( ( m_pRefcount == other.m_pRefcount ) == ( m_array == other.m_array ) );
	if ( m_pRefcount != other.m_pRefcount ) {

		assert( ( m_pRefcount == NULL ) == ( m_array == NULL ) );
		if ( m_pRefcount != NULL ) {

			assert( *m_pRefcount > 0 );
			if ( --( *m_pRefcount ) == 0 ) {

				delete m_pRefcount;
				delete[] m_array;
			}
		}

		m_pRefcount = other.m_pRefcount;
		m_array     = other.m_array;

		assert( ( m_pRefcount == NULL ) == ( m_array == NULL ) );
		if ( m_pRefcount != NULL )
			++( *m_pRefcount );
	}

	return *this;
}


template< typename t_Type >
t_Type& SharedArray< t_Type >::operator[]( unsigned int const index ) const {

	assert( ( m_pRefcount == NULL ) == ( m_array == NULL ) );
	assert( m_array != NULL );
	return m_array[ index ];
}


template< typename t_Type >
bool const SharedArray< t_Type >::operator==( t_Type const* const other ) const {

	return( m_array == other );
}


template< typename t_Type >
bool const SharedArray< t_Type >::operator==( SharedArray const& other ) const {

	assert( ( m_pRefcount == other.m_pRefcount ) == ( m_array == other.m_array ) );
	return( m_array == other.m_array );
}


template< typename t_Type >
bool const SharedArray< t_Type >::operator!=( t_Type const* const other ) const {

	return( m_array != other );
}


template< typename t_Type >
bool const SharedArray< t_Type >::operator!=( SharedArray const& other ) const {

	assert( ( m_pRefcount == other.m_pRefcount ) == ( m_array == other.m_array ) );
	return( m_array != other.m_array );
}




#endif    /* __SHARED_ARRAY_HH__ */
