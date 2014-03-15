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
                typedef void( const_iterator::*StepFunction )( void );
            public:
                const_iterator( const Array<T>& a, int64_t pos, int64_t begin, int64_t end ) : pos_( pos ), begin_( begin ), end_( end ), m_ConstArrayPtr( &a ), cdata( a.datablock.get() ) {
                    if( a.dense_ ) {
                        incrementor = &redux::util::Array<T>::const_iterator::fastIncrement;
                        decrementor = &redux::util::Array<T>::const_iterator::fastDecrement;
                    } else {
                        incrementor = &redux::util::Array<T>::const_iterator::sparseIncrement;
                        decrementor = &redux::util::Array<T>::const_iterator::sparseDecrement;
                        padding.resize(a.nDims_,a.dataSize); //a.dataSize-a.dimSizes[0]); // padding for the fastest dimension
                        for(size_t i=1; i<a.nDims_; ++i) {
                            padding[i] = a.dimStrides[i-1]-a.dimSizes[i]*a.dimStrides[i];
                        }
                       
                    }
                }

                // prefix
                const_iterator operator++() { ( this->*incrementor )(); return *this; }
                const_iterator operator--() { ( this->*decrementor )(); return *this; }
                // postfix
                const_iterator operator++( int ) { const_iterator i = *this; ( this->*incrementor )(); return i; }
                const_iterator operator--( int ) { const_iterator i = *this; ( this->*decrementor )(); return i; }

                const_reference operator*() { return *( cdata + pos_ ); }
                const_pointer operator->()  { return cdata + pos_; }

                int64_t pos( void ) { return pos_; };

                const_iterator& step( const std::vector<int64_t>& ind ) {
                    size_t nIndices = ind.size();
                    if( nIndices == 0 ) {
                        ++pos_;
                    }
                    else  if( nIndices > m_ConstArrayPtr->nDims_ ) {
                        throw std::out_of_range( "Array::iterator::step() called with more indices than dimensions." );
                    }
                    else {
                        int64_t delta = 0;
                        for( size_t i = 0; i < nIndices; ++i ) {
                            delta += ind[nIndices - i - 1] * m_ConstArrayPtr->dimStrides[m_ConstArrayPtr->nDims_ - i - 1];
                        }
                        pos_ += delta;
                    }
                    return *this;
                }
                template <typename ...S> const_iterator& step( S ...s ) { return step( {s...} ); }

                bool operator==( const const_iterator& rhs ) const { return pos_ == rhs.pos_; }
                bool operator!=( const const_iterator& rhs ) const { return pos_ != rhs.pos_; }
                bool operator>=( const const_iterator& rhs ) const { return pos_ >= rhs.pos_; }
                bool operator<=( const const_iterator& rhs ) const { return pos_ <= rhs.pos_; }
                bool operator> ( const const_iterator& rhs ) const { return pos_ > rhs.pos_; }
                bool operator< ( const const_iterator& rhs ) const { return pos_ < rhs.pos_; }

            protected:
                void fastDecrement( void ) { --pos_; }
                void fastIncrement( void ) { ++pos_; }

                void sparseDecrement( void ) {
                    --pos_;
                    for( int d = m_ConstArrayPtr->nDims_ - 1; d > 0; --d ) {
                        int64_t remainder = (pos_ % m_ConstArrayPtr->dimStrides[d-1])/m_ConstArrayPtr->dimStrides[d];
                        // if we went out-of-bounds for dimension d, subtract "padding"
                        if( remainder > m_ConstArrayPtr->dimLast[d] ||
                            remainder < m_ConstArrayPtr->dimFirst[d] ) {
                            pos_ -= padding[d];
                        }
                    }
                }

                void sparseIncrement( void ) {
                    ++pos_;
                    for( int d = m_ConstArrayPtr->nDims_ - 1; d > 0; --d ) {
                        int64_t remainder = (pos_ % m_ConstArrayPtr->dimStrides[d-1])/m_ConstArrayPtr->dimStrides[d];
                        // if we went out-of-bounds for dimension d, add "padding"
                        if( remainder > m_ConstArrayPtr->dimLast[d] ||
                            remainder < m_ConstArrayPtr->dimFirst[d] ) {
                            pos_ += padding[d];
                        }
                    }
                }

                int64_t pos_, begin_, end_;
                const Array<T>* m_ConstArrayPtr;
                const T* cdata;
                StepFunction incrementor, decrementor;
                std::vector<int64_t> padding;
            };

            class iterator : public const_iterator {
            public:
                iterator( Array<T>& a, int64_t pos, int64_t begin, int64_t end ) : const_iterator( a, pos, begin, end ), data( a.datablock.get() ) { }

                // prefix
                iterator operator++() { const_iterator::operator++(); return *this; }
                iterator operator--() { const_iterator::operator--(); return *this; }
                // postfix
                iterator operator++( int ) { iterator i = *this; const_iterator::operator++(); return i; }
                iterator operator--( int ) { iterator i = *this; const_iterator::operator--(); return i; }

                reference operator*() { return *( data + this->pos_ ); }
                pointer operator->() { return ( data + this->pos_ ); }

            private:
                T* data;
            };

            template <typename ...S> Array( S ...sizes ) : begin_(0) { resize( sizes... ); }

            template <typename ...S>
            Array( T* ptr, S ...sizes ) : begin_(0) {
                setSizes( sizes... );
                calcStrides();
                if( dataSize ) {
                    datablock.reset( ptr, []( T * p ) {} );
                }
            }

            /*! @brief Copy constructor and sub-array constructor, depending on the number of indices passed.
             */
            template <typename ...S>
            Array( Array<T>& rhs, S ...s ) : dense_(true), begin_(0) {

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
                        if( i > 0 && (dimStrides[i-1] > dimSizes[i]) ) dense_ = false;
                        nElements_ *= dimSizes[i];
                    }
                    begin_ = getOffset(dimFirst);
                    auto it = const_iterator(*this, getOffset(dimLast), begin_, getOffset(dimLast)+nElements_);
                    end_ = (++it).pos();    // 1 step after last element;
                }
                else {  // ordinary copy constructor
                    dimSizes = rhs.dimSizes;
                    dimFirst = rhs.dimFirst;
                    dimLast = rhs.dimLast;
                    nElements_ = rhs.nElements_;
                    begin_ = rhs.nElements_;
                    end_ = rhs.end_;
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
            template <typename ...S> void resize( S ...sizes ) { resize( {static_cast<size_t>( sizes )...} ); }

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

            template <typename U = T>
            Array<U> copy( void ) const {
                std::vector<size_t> newDimSizes;
                for( auto & it : dimSizes ) {
                    if( it > 1 ) {
                        newDimSizes.push_back( it );
                    }
                }
                Array<U> tmp( newDimSizes );
                const_iterator cit = begin();
                for( auto & it : tmp ) {
                    it = static_cast<U>( *cit );
                    ++cit;
                }
                return tmp;
            }

            template <typename U>
            void copy( Array<U> arr ) const {
                std::vector<size_t> newDimSizes;
                for( auto & it : dimSizes ) {
                    if( it > 1 ) {
                        newDimSizes.push_back( it );
                    }
                }
                arr.resize( newDimSizes );
                const_iterator cit = begin();
                for( auto & it : arr ) {
                    it = static_cast<U>( *cit );
                    ++cit;
                }

            }

            template <typename ...S>
            T* ptr( S ...s ) {
                int64_t offset = getOffset<int64_t>( {static_cast<int64_t>( s )...}, dimFirst );
                if( offset < 0 || offset > dataSize ) {
                    throw std::out_of_range( "Offset out of range: " + std::to_string( offset ) );
                }
                return ( datablock.get() + offset );
            }
            template <typename ...S> const T* ptr( S ...s ) const { return const_cast<Array<T>*>( this )->ptr( s... ); }

            T* ptr( const std::vector<int64_t>& indices ) {
                int64_t offset = getOffset( indices, dimFirst );
                if( offset < 0 || offset > dataSize ) {
                    throw std::out_of_range( "Offset out of range: " + std::to_string( offset ) );
                }
                return ( datablock.get() + offset );
            }
            const T* ptr( const std::vector<int64_t>& indices ) const { return const_cast<Array<T>*>( this )->ptr( indices ); }

            operator bool() const { return static_cast<bool>( datablock ); }

            template <typename ...S> T& operator()( S ...indices ) { return *ptr( indices... ); }
            template <typename ...S> const T& operator()( S ...indices ) const { return *ptr( indices... ); }

            template <typename ...S>
            T& at( S ...indices ) {
                std::vector<int64_t> tmp({indices...});
                size_t sz = tmp.size();
                if( sz > nDims_ ) {
                    throw std::out_of_range( "Array::at() called with more indices than dimensions: " + std::to_string( sz ) );
                }
                else if( sz == 0 ) {
                    return *datablock.get();
                }
                for( size_t i = 0; i < sz; ++i ) {
                    if( tmp[i] >= dimSizes[nDims_ - i - 1] ) {
                        throw std::out_of_range( "Index out of range. " + printArray(tmp,"indices") );
                    }
                }
                int64_t offset = getOffset( tmp, dimFirst );
                return *( datablock.get() + offset );
            }
            template <typename ...S> const T& at( S ...indices ) const { return const_cast<Array<T>*>( this )->at( indices... ); }

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
                const_iterator rhsit = rhs.begin();
                for( const auto & it : *this ) {
                    if( it != *rhsit++ ) return false;
                }

                return true;
            }
            bool operator!=( const Array<T>& rhs ) const { return !( *this == rhs ); }

            template <typename ...S>
            iterator pos( S ...indices ) {
                return iterator( *this, getOffset<int64_t>( {static_cast<int64_t>( indices )...}, dimFirst ), begin_, end_ );
            }
            template <typename ...S>
            const_iterator pos( S ...indices ) const { return const_cast<Array<T>*>( this )->pos( indices... ); }

            iterator begin( void ) {
                return iterator( *this, begin_, begin_, end_ );
            }
            const_iterator begin( void ) const { return const_iterator( *this, begin_, begin_, end_ ); }

            iterator end( void ) { return iterator( *this, end_ , begin_, end_ ); }
            const_iterator end( void ) const { return const_iterator( *this, end_ , begin_, end_ ); }

            T* get( void ) { return datablock.get(); };
            const T* get( void ) const { return datablock.get(); };

            bool dense( void ) const { return dense_; }

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
                dense_ = true;
                dataSize = 1;
                for( auto & it : dimLast ) {
                    dataSize *= it;
                    it--;
                }
                begin_ = 0;
                end_ = dataSize;
            }
            template <typename ...S> void setSizes( S ...sizes ) { setSizes( {static_cast<size_t>( sizes )...} );  }

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
                if( dimSizes.size() > 0 ) {
                    nElements_ = 1;
                    for( int i = dimSizes.size() - 1; i > 0; --i ) {
                        dimStrides[i - 1] = dimSizes[i] * dimStrides[i];
                        nElements_ *= ( dimLast[i] - dimFirst[i] + 1 );
                    }
                    nElements_ *= ( dimLast[0] - dimFirst[0] + 1 );
                }
            }

            template <typename U>
            int64_t getOffset( const std::vector<U>& indices, const std::vector<size_t>& offsets ) const {
                int64_t offset = 0;
                if( indices.size() > offsets.size() ) {
                    return getOffset( indices );
                }
                else if( indices.size() > dimSizes.size() ) {
                    throw std::logic_error( "Array::getOffset() called with more indices than dimensions: " + printArray( indices, "indices" ) );
                }
                else {
                    size_t dimDiff = nDims_ - indices.size();
                    for( size_t i = 0; i < indices.size(); ++i ) {
                        int64_t tmp = indices[i] + offsets[i];
                        offset += tmp * dimStrides[dimDiff + i];
                    }
                }
                return offset;
            }

            template <typename U>
            int64_t getOffset( const std::vector<U>& indices ) const {
                int64_t offset = 0;
                if( indices.size() > dimSizes.size() ) {
                    throw std::logic_error( "Array::getOffset() called with more indices than dimensions: " + printArray( indices, "indices" ) );
                }
                else {
                    size_t dimDiff = nDims_ - indices.size();
                    for( size_t i = 0; i < indices.size(); ++i ) {
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
            bool dense_;
            int64_t begin_, end_;

            template<typename U> friend class Array;
            friend class const_iterator;
            friend class iterator;

        };


        /*! @} */


    }

}

#endif // REDUX_UTIL_ARRAY_HPP

