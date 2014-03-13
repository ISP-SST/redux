#ifndef REDUX_UTIL_ARRAY_HPP
#define REDUX_UTIL_ARRAY_HPP

#include "redux/util/convert.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"

#include <stdexcept>
#include <memory>
#include <cstring>
#include <vector>
#include <iostream>

namespace redux {

    namespace util {

        /*! @defgroup util Util
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
                typedef void(const_iterator::*StepFunction)(void);
            public:
                const_iterator( const Array<T>& a, size_t begin, size_t end ) : current_( begin ), end_( end ), m_ConstArrayPtr( &a ), cdata( a.datablock.get() ) {

                    if( m_ConstArrayPtr->nElements_ == m_ConstArrayPtr->dataSize ) {    // dense array
                        incrementor = &redux::util::Array<T>::const_iterator::fastIncrement;
                        decrementor = &redux::util::Array<T>::const_iterator::fastDecrement;
                    } else {                                                            // segmented array
                        if (begin < end) indices = m_ConstArrayPtr->dimFirst;
                        else indices = m_ConstArrayPtr->dimLast;
                        incrementor = &redux::util::Array<T>::const_iterator::sparseIncrement;
                        decrementor = &redux::util::Array<T>::const_iterator::sparseDecrement;
                    }
                }
                
                // prefix
                const_iterator operator++() { (this->*incrementor)(); return *this; }
                const_iterator operator--() { (this->*decrementor)(); return *this; }
                // postfix
                const_iterator operator++( int ) { const_iterator i = *this; (this->*incrementor)(); return i; }
                const_iterator operator--( int ) { const_iterator i = *this; (this->*decrementor)(); return i; }

                const_reference operator*() { return *( cdata + current_ ); }
                const_pointer operator->()  { return cdata + current_; }
                bool operator==( const const_iterator& rhs ) const { return current_ == rhs.current_; }
                bool operator!=( const const_iterator& rhs ) const { return current_ != rhs.current_; }
                bool operator>=( const const_iterator& rhs ) const { return current_ >= rhs.current_; }
                bool operator<=( const const_iterator& rhs ) const { return current_ <= rhs.current_; }
                bool operator> ( const const_iterator& rhs ) const { return current_ > rhs.current_; }
                bool operator< ( const const_iterator& rhs ) const { return current_ < rhs.current_; }

            protected:
                void fastDecrement( void ) {
                    if( current_ == 0 ) {
                        throw std::out_of_range( "Iterating past begin()." );
                    }
                    if( current_ > end_ ) {
                        current_ = end_;
                    }
                    else {
                        --current_;
                    }
                }
                void fastIncrement( void ) {
                    if( current_ >= end_ ) {
                        throw std::out_of_range( "Iterating past end()." );
                    }
                    current_++;
                }

                void sparseDecrement( void ) {
                    if( current_ == 0 ) {
                        throw std::out_of_range( "Iterating past begin()." );
                    } else if ( current_ >= end_ ) {
                        current_ = end_ -1;
                        indices = m_ConstArrayPtr->dimLast;
                        return;
                    }
                    for( int d = m_ConstArrayPtr->nDims_ - 1; d >= 0; --d ) {
                        if( indices[d] == m_ConstArrayPtr->dimFirst[d] ) {
                            if( d == 0 ) {
                                current_ = 0;
                                return;
                            }
                            indices[d] = m_ConstArrayPtr->dimLast[d];
                            current_ += ( m_ConstArrayPtr->dimSizes[d] * m_ConstArrayPtr->dimStrides[d]);
                        } else {
                            indices[d]--;
                            while (d<m_ConstArrayPtr->nDims_) {
                                current_ -= m_ConstArrayPtr->dimStrides[d++];
                            }
                            break;
                        }
                    }
                }
                
                void sparseIncrement( void ) {

                    if( current_ >= end_ ) {
                        throw std::out_of_range( "Iterating past end()." );
                    }
                    for( int d = m_ConstArrayPtr->nDims_ - 1; d >= 0; --d ) {
                        if( indices[d] == m_ConstArrayPtr->dimLast[d] ) {
                            if( d == 0 ) {
                                current_ = end_;
                                return;
                            }
                            indices[d] = m_ConstArrayPtr->dimFirst[d];
                            current_ -= ( m_ConstArrayPtr->dimSizes[d] * m_ConstArrayPtr->dimStrides[d]);
                        } else {
                            indices[d]++;
                            while (d<m_ConstArrayPtr->nDims_) {
                                current_ += m_ConstArrayPtr->dimStrides[d++];
                            }
                            break;
                        }
                    }
                }

                size_t current_, end_;
                const Array<T>* m_ConstArrayPtr;
                const T* cdata;
                StepFunction incrementor, decrementor;
                std::vector<size_t> indices;
            };

            class iterator : public const_iterator {
            public:
                iterator( Array<T>& a, size_t begin, size_t end ) : const_iterator( a, begin, end ), data( a.datablock.get() ) {}

                // prefix
                iterator operator++() { const_iterator::operator++(); return *this; }
                iterator operator--() { const_iterator::operator--(); return *this; }
                // postfix
                iterator operator++( int ) { iterator i = *this; const_iterator::operator++(); return i; }
                iterator operator--( int ) { iterator i = *this; const_iterator::operator--(); return i; }

                reference operator*() { return *( data + this->current_ ); }
                pointer operator->() { return ( data + this->current_ ); }

            private:
                T* data;
            };

            template <typename ...S>
            Array( S ...sizes ) {
                resize( sizes... );
            }

            template <typename ...S>
            Array( T* ptr, S ...sizes ) {
                setSizes( sizes... );
                calcStrides();
                if( dataSize ) {
                    datablock.reset( ptr, []( T * p ) {} );
                }
            }

            /*! @brief Copy constructor and sub-array constructor, depending on the number of indices passed.
             */
            template <typename ...S>
            Array( Array<T>& rhs, S ...s ) {

                nDims_ = rhs.dimSizes.size();
                dataSize = rhs.dataSize;
                datablock = rhs.datablock;
                dimStrides = rhs.dimStrides;
                nElements_ = 0;
                if( ( sizeof...( s ) ) > 0 ) {
                    std::vector<int64_t> tmp = {static_cast<int64_t>( s )...};
                    dimFirst.resize( nDims_ );
                    dimLast.resize( nDims_ );
                    dimSizes.resize( nDims_ );
                    nElements_ = 1;
                    for( size_t i = 0; i < nDims_; ++i ) {
                        dimFirst[i] = redux::util::bound_cast<size_t>( rhs.dimFirst[i] + tmp[2 * i], 0, rhs.dimSizes[i] );
                        dimLast[i] = redux::util::bound_cast<size_t>( rhs.dimFirst[i] + tmp[2 * i + 1], 0, rhs.dimSizes[i] );
                        if( dimFirst[i] > dimLast[i] ) std::swap( dimFirst[i], dimLast[i] );
                        dimSizes[i] = ( dimLast[i] - dimFirst[i] + 1 );
                        if( !dimSizes[i] ) {
                            throw std::logic_error( "Attempting to create a sub-array with a dimension of 0 size. " + printArray( tmp, "dims" ) );
                        }
                        nElements_ *= dimSizes[i];
                    }

                }
                else {  // ordinary copy constructor
                    dimSizes = rhs.dimSizes;
                    dimFirst = rhs.dimFirst;
                    dimLast = rhs.dimLast;
                    nElements_ = rhs.nElements_;
                }
            }

            size_t size( void ) const {
                size_t sz = sizeof( size_t );                                                          // blockSize
                if( dimSizes.size() ) {
                    sz += 1 + dimSizes.size() * sizeof( size_t ) + nElements_ * sizeof( T );            // nDims + dimensions + dataSize
                }
                return sz;
            }

            char* pack( char* ptr ) const {
                using redux::util::pack;
                ptr = pack( ptr, size() );
                if( uint8_t ndims = dimSizes.size() ) {
                    ptr = pack( ptr, ndims );
                    ptr = pack( ptr, dimSizes );
                    ptr = pack( ptr, datablock.get(), nElements_ );
                }
                return ptr;

            }

            const char* unpack( const char* ptr, bool swap_endian ) {
                using redux::util::unpack;
                size_t sz;
                ptr = unpack( ptr, sz, swap_endian );
                if( sz > 8 ) {
                    uint8_t ndims;
                    ptr = unpack( ptr, ndims );
                    std::vector<size_t> tmp( ndims );
                    ptr = unpack( ptr, tmp, swap_endian );
                    resize( tmp );
                    ptr = unpack( ptr, datablock.get(), nElements_, swap_endian );
                }
                return ptr;
            }


            void resize( const std::vector<size_t>& sizes ) {
                setSizes( sizes );
                calcStrides();
                create();
            }

            template <typename ...S>
            void resize( S ...sizes ) { resize( {static_cast<size_t>( sizes )...} ); }

            const std::vector<size_t>& dimensions( void ) const { return dimSizes; }
            size_t nDimensions( void ) const { return dimSizes.size(); }
            size_t nElements( void ) const { return nElements_; }

            size_t dimSize( size_t i = 0 ) const {
                if( i < dimSizes.size() ) return dimSizes[i];
                else return 0;
            }

            T* cloneData( void ) {
                T* data = nullptr;
                if( nElements_ ) {
                    data = new T[nElements_];
                    T* tmpPtr = data;
                    for( auto & it : *this ) {
                        *tmpPtr++ = it;
                    }
                }
                return data;
            }

            template <typename U>
            void copyFrom( void* data ) {
                U* ptr = reinterpret_cast<U*>( data );
                for( auto & it : *this ) {
                    it = static_cast<T>( *ptr++ );
                }
            }

            template <typename U=T>
            Array<U> copy( void ) const {
                std::vector<size_t> newDimSizes;
                for (auto& it: dimSizes) {
                    if(it>1) {
                        newDimSizes.push_back(it);
                    }
                }
                Array<U> tmp( newDimSizes );
                const_iterator cit = begin();
                for( auto& it : tmp ) {
                    it = static_cast<U>( *cit );
                    ++cit;
                }
                return tmp;
            }

            template <typename U>
            void copy( Array<U> arr ) const {
                std::vector<size_t> newDimSizes;
                for (auto& it: dimSizes) {
                    if(it>1) {
                        newDimSizes.push_back(it);
                    }
                }
                arr.resize( newDimSizes );
                const_iterator cit = begin();
                for( auto& it : arr ) {
                    it = static_cast<U>( *cit );
                    ++cit;
                }

            }

            template <typename ...S>
            T* ptr( S ...s ) {
                return ( datablock.get() + getOffset( {static_cast<size_t>( s )...}, dimFirst ) );
            }

            template <typename ...S>
            const T* ptr( S ...s ) const {
                return ( datablock.get() + getOffset( {static_cast<size_t>( s )...}, dimFirst ) );
            }

            T* ptr( const std::vector<size_t>& indices ) {
                return ( datablock.get() + getOffset( indices, dimFirst) );
            }

            const T* ptr( const std::vector<size_t>& indices ) const {
                return ( datablock.get() + getOffset( indices, dimFirst ) );
            }

            operator bool() const { return static_cast<bool>(datablock); }

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
                if( values.size() == nElements_ ) {
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

            // operators with scalars
            const Array<T>& operator= ( const T& rhs ) { for( auto & it : *this ) it  = rhs; return *this; };
            const Array<T>& operator+=( const T& rhs ) { for( auto & it : *this ) it += rhs; return *this; };
            const Array<T>& operator-=( const T& rhs ) { for( auto & it : *this ) it -= rhs; return *this; };
            const Array<T>& operator*=( const T& rhs ) { for( auto & it : *this ) it *= rhs; return *this; };
            const Array<T>& operator/=( const T& rhs ) { for( auto & it : *this ) it /= rhs; return *this; };
            virtual void zero( void ) { if( nElements_ ) memset( datablock.get(), 0, nElements_ * sizeof( T ) ); }

            template <typename U>
            const Array<T>& operator=( const Array<U>& rhs ) {
                if( !sameSizes( rhs ) ) this->resize( rhs.dimSizes );
                typename Array<U>::const_iterator rhsit = rhs.begin();
                for( auto & it : *this ) it = redux::util::bound_cast<T>( *rhsit++ );
                return *this;
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
                return iterator( *this, getOffset(dimFirst), getOffset( dimLast ) + 1 );
            }
            const_iterator begin( void ) const {
                return const_iterator( *this, getOffset(dimFirst), getOffset( dimLast ) + 1 );
            }
            iterator end( void ) {
                size_t lastp1 = getOffset( dimLast ) + 1;
                return iterator( *this, lastp1, lastp1 );
            }
            const_iterator end( void ) const {
                size_t lastp1 = getOffset( dimLast ) + 1;
                return const_iterator( *this, lastp1, lastp1 );
            }

            T* get( void ) { return datablock.get(); };
            const T* get( void ) const { return datablock.get(); };

            bool dense(void) const { return (nElements_ == dataSize); }
            
            template <typename U>
            bool sameSize( const Array<U>& rhs ) const {
                return ( nElements_ == rhs.nElements_ );
            }

            template <typename U>
            bool sameSizes( const Array<U>& rhs ) const {
                if( nDims_ != rhs.nDims_ ) {
                    return false;
                }
                for( size_t i = 0; i < nDims_; ++i ) {
                    if( ( dimSizes[i] ) != ( rhs.dimSizes[i] ) ) {
                        return false;
                    }
                }
                return true;
            }
            
        private:
            void setSizes( const std::vector<size_t>& sizes ) {
                dimSizes = sizes;
                nDims_ = dimSizes.size();
                dimFirst.resize( nDims_, 0 );
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
            }

            void calcStrides( void ) {
                dimStrides.resize( dimSizes.size(), 1 );
                nElements_ = 0;
                dataSize = 0;
                if( dimSizes.size() > 0 ) {
                    nElements_ = 1;
                    for( int i = dimSizes.size() - 1; i > 0; --i ) {
                        dimStrides[i - 1] = dimSizes[i] * dimStrides[i];
                        nElements_ *= ( dimLast[i] - dimFirst[i] + 1 );
                    }
                    nElements_ *= ( dimLast[0] - dimFirst[0] + 1 );
                    dataSize = dimSizes[0] * dimStrides[0];
                }
            }

            size_t getOffset( const std::vector<size_t>& indices, const std::vector<size_t>& offsets ) const {
                size_t offset = 0;
                if( indices.size() > dimSizes.size() ) {
                    throw std::logic_error( "More indices than dimensions." );
                }
                else {
                    size_t dimDiff = nDims_ - indices.size();
                    for( size_t i = 0; i < indices.size(); ++i ) {
                        size_t tmp = indices[i] + offsets[i];
                        if( tmp > dimLast[i] ) {
                            throw std::out_of_range( "Index out of range: " + std::to_string(tmp)  );
                        }
                        offset += tmp * dimStrides[dimDiff + i];
                    }
                }
                return offset;
            }
            
            size_t getOffset( const std::vector<size_t>& indices ) const {
                size_t offset = 0;
                if( indices.size() > dimSizes.size() ) {
                    throw std::logic_error( "More indices than dimensions." );
                }
                else {
                    size_t dimDiff = nDims_ - indices.size();
                    for( size_t i = 0; i < indices.size(); ++i ) {
                        if(  indices[i] > dimLast[i] ) {
                            throw std::out_of_range( "Index out of range: " + std::to_string( indices[i])  );
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
            size_t nDims_;
            size_t nElements_;
            size_t dataSize;

            template<typename U> friend class Array;
            friend class const_iterator;
            friend class iterator;

        };


        /*! @} */


    }

}

#endif // REDUX_UTIL_ARRAY_HPP

