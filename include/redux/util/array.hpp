#ifndef REDUX_UTIL_ARRAY_HPP
#define REDUX_UTIL_ARRAY_HPP

#include "redux/util/convert.hpp"

#include <stdexcept>
#include <memory>
#include <cstring>
#include <vector>
#include <iostream>

namespace redux {

    namespace util {

        /*!  @ingroup util
         *  @{
         */


        class IterationOutOfBounds : public std::out_of_range {
        public:
            IterationOutOfBounds( void ) : out_of_range( "Array::BadIteration" ) {}
            IterationOutOfBounds( const std::string &message ) : out_of_range( message ) {}

            virtual ~IterationOutOfBounds( void ) throw() {}
        };


        /*! @brief A wrapper class for a data-block of arbitrary dimensions/type
         *  @details The data is stored in a shared_ptr and accessed via iterators or some (variadic template) methods
         *  @todo Implement iterators for sparse data (i.e. when it's a sub-array accessing a subset of the datablock.
         */
        template <class T = double>
        class Array {

        public:

            typedef size_t size_type;
            typedef T value_type;
            typedef T& reference;
            typedef T* pointer;
            typedef const T* const_pointer;
            typedef const T& const_reference;
            typedef std::bidirectional_iterator_tag iterator_category;
            typedef ptrdiff_t difference_type;

            class const_iterator {
            public:
                const_iterator( const Array<T>& a, size_t begin, size_t end ) : m_Begin( begin ), m_Current( begin ),
                    m_End( end ), m_ConstArrayPtr( &a ), m_ConstData( a.data.get() ) {
                    if( m_Current > m_End ) throw std::out_of_range( "Array::const_iterator:  begin>end" );
                }
                // prefix
                const_iterator operator++() { increment(); return *this; }
                const_iterator operator--() { decrement(); return *this; }
                // postfix
                const_iterator operator++( int ) { const_iterator i = *this; increment(); return i; }
                const_iterator operator--( int ) { const_iterator i = *this; decrement(); return i; }

                const_reference operator*() { return *( m_ConstData + m_Current ); }
                const_pointer operator->() { return m_ConstData + m_Current; }
                bool operator==( const const_iterator& rhs ) const { return m_Current == rhs.m_Current; }
                bool operator!=( const const_iterator& rhs ) const { return m_Current != rhs.m_Current; }

            protected:
                void decrement( void ) {
                    if( m_Current <= m_Begin ) {
                        throw std::out_of_range( "Iterating past begin()." );
                    }
                    if( m_Current > m_End ) {
                        m_Current = m_End;
                    }
                    else {
                        --m_Current;
                    }
                }
                void increment( void ) {
                    if( m_Current >= m_End ) {
                        throw std::out_of_range( "Iterating past end()." );
                    }
                    ++m_Current;

                }

                size_t m_Begin, m_Current, m_End;
                const Array<T>* m_ConstArrayPtr;
                const T* m_ConstData;

            };

            class iterator : public const_iterator {
            public:
                iterator( Array<T>& a, size_t begin, size_t end ) : const_iterator( a, begin, end ),
                    m_ArrayPtr( &a ), data( a.data.get() ) {}

                // prefix
                iterator operator++() { this->increment(); return *this; }
                iterator operator--() { this->decrement(); return *this; }
                // postfix
                iterator operator++( int ) { iterator i = *this; this->increment(); return i; }
                iterator operator--( int ) { iterator i = *this; this->decrement(); return i; }

                reference operator*() { return *( data + this->m_Current ); }
                pointer operator->() { return ( data + this->m_Current ); }

            private:
                Array<T>* m_ArrayPtr;
                T* data;
            };

            template <typename ...S>
            Array( S ...sizes ) {
                reset( sizes... );
            }

            template <typename ...S>
            Array( T* ptr, S ...sizes ) {
                setSizes( sizes... );
                calcStrides();
                if( dataSize ) {
                    datablock.reset( ptr, []( T * p ) {} );
                }
                denseData = true;
            }

            /*! @brief Copy constructor and sub-array constructor, depending on the number of indices passed.
             */
            template <typename ...S>
            Array( Array<T>& rhs, S ...s ) {

                nDims = rhs.dimSizes.size();
                dataSize = rhs.dataSize;
                datablock = rhs.datablock;
                dimSizes = rhs.dimSizes;
                dimStrides = rhs.dimStrides;
                if( nDims && ( sizeof...( s ) ) == 2 * nDims ) {
                    //std::vector<int64_t> tmp = {static_cast<int64_t>( s )...};
                    int64_t tmp[] = {s...};
                    dimFirst.resize( nDims );
                    dimLast.resize( nDims );
                    nElements = 1;
                    for( size_t i = 0; i < nDims; ++i ) {
                        if( !tmp[nDims + i] ) {
                            throw std::logic_error( "Attempting to create a sub-array with a dimension of 0 size." );
                        }
                        dimFirst[i] = redux::util::bound_cast<size_t>( rhs.dimFirst[i] + tmp[i], 0, rhs.dimSizes[i] );
                        dimLast[i] = redux::util::bound_cast<size_t>( dimFirst[i] + tmp[nDims + i] - 1, 0, rhs.dimSizes[i] );
                        if( dimFirst[i] > dimLast[i] ) std::swap( dimFirst[i], dimLast[i] );
                        nElements *= ( dimLast[i] - dimFirst[i] + 1 );
                    }

                }
                else {  // ordinary copy constructor
                    dimFirst = rhs.dimFirst;
                    dimLast = rhs.dimLast;
                    nElements = rhs.m_nElements;
                }
                denseData = !( rhs.denseData && ( nElements == rhs.m_nElements ) );
            }

            void reset( const std::vector<size_t>& sizes ) {
                setSizes( sizes );
                calcStrides();
                create();
            }

            template <typename ...S>
            void reset( S ...sizes ) { reset( {static_cast<size_t>( sizes )...} ); }

            size_t size( size_t i = 0 ) {
                if( i < dimSizes.size() ) return dimSizes[i];
                else return 0;
            }

            T* cloneData( void ) {
                T* data = nullptr;
                if( nElements ) {
                    data = new T[nElements];
                    T* tmpPtr = data;
                    for( auto & it : *this ) {
                        *tmpPtr++ = it;
                    }
                }
                return data;
            }

            Array<T> copy( void ) {
                std::vector<size_t> tmpDims = dimLast;
                for( size_t i = 0; i < tmpDims.size(); ++i ) {
                    tmpDims[i] -= ( dimFirst[i] - 1 );
                }
                Array<T> tmp( tmpDims );
                const_iterator it = begin();
                iterator tmpit = tmp.begin();
                for( size_t i = 0; i < nElements; ++i ) {
                    *tmpit++ = *it++;
                }
                return tmp;
            }

            template <typename U>
            void rawCopy( void* ) {
                // TODO implement
            }

            template <typename ...S>
            T* ptr( S ...s ) {
                return ( datablock.get() + getOffset( {static_cast<size_t>( s )...} ) );
            }

            template <typename ...S>
            const T* ptr( S ...s ) const {
                return ( datablock.get() + getOffset( {static_cast<size_t>( s )...} ) );
            }

            T* ptr( const std::vector<size_t>& indices ) {
                return ( datablock.get() + getOffset( indices ) );
            }

            const T* ptr( const std::vector<size_t>& indices ) const {
                return ( datablock.get() + getOffset( indices ) );
            }

            template <typename ...S>
            T& operator()( S ...s ) {
                return *ptr( s... );
            }

            template <typename ...S>
            const T& operator()( S ...s ) const {
                return *ptr( s... );
            }

            template <typename ...S>
            void setBlock( T* data, size_t count, S ...s ) const {
                memcpy( ptr( s... ), data, count * sizeof( T ) );
            }

            void set( const std::vector<T>& values ) {
                if( values.size() == nElements ) {
                    T* ptr = datablock.get();
                    for( auto & i : values ) {
                        *ptr++ = i;
                    }
                }
                else {
                    throw std::logic_error( "Number of total elements must match when setting a list of values." );
                }
            }

            template <typename ...S>
            void set( S ...values ) {   set( {static_cast<T>( values )...} ); }

            void set( Array<T>& rhs ) {
                if( ( this != &rhs ) && sameSizes( rhs ) ) {
                    iterator it = begin();
                    for( auto & rhsit : rhs ) {
                        *it++ = rhsit;
                    }
                }
            }

            bool operator==( const Array<T>& rhs ) const {
                if( ! sameSizes( rhs ) ) {
                    return false;
                }
                if( ptr() == rhs.ptr() ) {      // shared data, check if it is the same sub-array.
                    if( dimFirst == rhs.dimFirst ) {
                        return true;
                    }
                }
                // compare data.
                const_iterator it = begin();
                const_iterator theend = end();
                const_iterator it2 = rhs.begin();
                while( it != theend ) {
                    if( *it++ != *it2++ ) return false;
                }

                return true;
            }
            bool operator!=( const Array<T>& rhs ) const { return !( *this == rhs ); }

            iterator begin( void ) {
                size_t first = getOffset( {0} );
                if( denseData ) {
                    return iterator( *this, first, first + nElements );
                }
                return iterator( *this, first, first + nElements );
            }
            const_iterator begin( void ) const {
                size_t first = getOffset( {0} );
                if( denseData ) {
                    return const_iterator( *this, first, first + nElements );
                }
                return const_iterator( *this, first, first + nElements );

            }
            iterator end( void ) {
                size_t first = getOffset( {0} );
                return iterator( *this, first + nElements, first + nElements );
            }
            const_iterator end( void ) const {
                size_t first = getOffset( {0} );
                return const_iterator( *this, first + nElements, first + nElements );
            }

            T* get( void ) { return datablock.get(); };
            const T* get( void ) const { return datablock.get(); };

        private:
            template <typename U>
            bool sameSizes( const Array<U>& rhs ) const {
                if( nDims != rhs.m_nDims ) {
                    return false;
                }
                for( size_t i = 0; i < nDims; ++i ) {
                    if( ( dimLast[i] - dimFirst[i] ) != ( rhs.dimLast[i] - rhs.dimFirst[i] ) ) {
                        return false;
                    }
                }
                return true;
            }

            void setSizes( const std::vector<size_t>& sizes ) {
                dimSizes = sizes;
                nDims = dimSizes.size();
                dimFirst.resize( nDims, 0 );
                dimLast = dimSizes;
                for( auto & it : dimLast ) it--;
            }

            template <typename ...S>
            void setSizes( S ...sizes ) { setSizes( {static_cast<size_t>( sizes )...} );  }

            void create( void ) {
                if( dataSize ) {
                    datablock.reset( new T[dataSize], []( T * p ) { delete[] p; } );
                }
                else {
                    datablock.reset();
                }
                denseData = true;
            }

            void calcStrides( void ) {
                dimStrides.resize( dimSizes.size(), 1 );
                nElements = 0;
                dataSize = 0;
                if( dimSizes.size() > 0 ) {
                    nElements = 1;
                    for( int i = dimSizes.size() - 1; i > 0; --i ) {
                        dimStrides[i - 1] = dimSizes[i] * dimStrides[i];
                        nElements *= ( dimLast[i] - dimFirst[i] + 1 );
                    }
                    nElements *= ( dimLast[0] - dimFirst[0] + 1 );
                    dataSize = dimSizes[0] * dimStrides[0];
                }
            }

            size_t getOffset( const std::vector<size_t>& indices ) const {
                size_t offset = 0;
                if( indices.size() > dimSizes.size() ) {
                    throw std::logic_error( "More indices than dimensions." );
                }
                else if( indices.size() ) {
                    size_t dimDiff = nDims - indices.size();
                    for( size_t i = 0; i < indices.size(); ++i ) {
                        if( indices[i] > dimLast[dimDiff + i] ) {
                            throw std::out_of_range( indices[i], dimLast[dimDiff + i] );
                        }
                        offset += indices[i] * dimStrides[dimDiff + i];
                    }
                }
                return offset;
            }

            std::vector<size_t> dimSizes;
            std::vector<size_t> dimStrides;

            std::vector<size_t> dimFirst;
            std::vector<size_t> dimLast;

            std::shared_ptr<T> datablock;
            bool denseData;
            size_t nDims;
            size_t nElements;
            size_t dataSize;

            friend class const_iterator;
            friend class iterator;

        };


        /*! @} */


    }

}

#endif // REDUX_UTIL_ARRAY_HPP

