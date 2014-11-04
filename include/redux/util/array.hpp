#ifndef REDUX_UTIL_ARRAY_HPP
#define REDUX_UTIL_ARRAY_HPP

#include "redux/util/arrayutil.hpp"
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


        /*! @brief A wrapper class for a data-block of arbitrary dimensions/type
         *  @details The data is stored in a shared_ptr and accessed via iterators or some (variadic template) methods
         *  @note This class is intended for convenient access to sub-arrays, it is not optimized for speed.
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
                const_iterator& operator++() { ( this->*incrementor )(); return *this; }
                const_iterator& operator--() { ( this->*decrementor )(); return *this; }
                // postfix
                const_iterator operator++( int ) { const_iterator i = *this; ( this->*incrementor )(); return i; }
                const_iterator operator--( int ) { const_iterator i = *this; ( this->*decrementor )(); return i; }

                const_reference operator*() const { return *( cdata + pos_ ); }
                const_pointer operator->() const { return cdata + pos_; }
                operator const_pointer() const { return ( cdata + pos_ ); }

                const_iterator& operator+=( const const_iterator& rhs ) { pos_ += rhs.pos_; trimUp(); return *this; }
                const_iterator& operator-=( const const_iterator& rhs ) { pos_ -= rhs.pos_; trimUp(); return *this; }
                
                const_iterator operator+( difference_type step ) const { const_iterator i = *this; i.pos_ += step; i.trimUp(); return i; }
                const_iterator operator-( difference_type step ) const { const_iterator i = *this; i.pos_ -= step; i.trimUp(); return i; }
                
                difference_type operator+( const const_iterator& rhs ) const { return (this->pos_ + rhs.pos()); }
                difference_type operator-( const const_iterator& rhs ) const { return (this->pos_ - rhs.pos()); }
                
                
                int64_t pos( void ) const { return pos_; }
                template <typename U>
                const Array<T> neighbourhood(const std::vector<U>& maxDistance) {
                    std::vector<int64_t> first = m_ConstArrayPtr->indexFromPos(pos_);
                    if( first.size() != maxDistance.size() ) {
                        throw std::logic_error("Array::iterator::neighbourhood() Dimension mismatch:  " + printArray(maxDistance,"maxDistance"));
                    }
                    std::vector<int64_t> last = first;
                    for (size_t i=0; i<first.size(); ++i) {
                        first[i] -= maxDistance[i];
                        last[i]  += maxDistance[i];
                    }
                    return Array<T>(*const_cast<Array<T>*>(m_ConstArrayPtr),first,last);
                }
                template <typename ...S> const Array<T> neighbourhood( S ...s ) { return neighbourhood<int64_t>({static_cast<int64_t>(s)...}); }

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
                    trimDown();
                }

                void sparseIncrement( void ) {
                    ++pos_;
                    trimUp();
                }

                void trimUp(void) {
                    for( int d = m_ConstArrayPtr->nDims_ - 1; d > 0; --d ) {
                        int64_t remainder = (pos_ % m_ConstArrayPtr->dimStrides[d-1])/m_ConstArrayPtr->dimStrides[d];
                        // if we went out-of-bounds for dimension d, add "padding"
                        if( remainder > m_ConstArrayPtr->dimLast[d] ||
                            remainder < m_ConstArrayPtr->dimFirst[d] ) {
                            pos_ += padding[d];
                        }
                    }

                }
                void trimDown(void) {
                    for( int d = m_ConstArrayPtr->nDims_ - 1; d > 0; --d ) {
                        int64_t remainder = (pos_ % m_ConstArrayPtr->dimStrides[d-1])/m_ConstArrayPtr->dimStrides[d];
                        // if we went out-of-bounds for dimension d, subtract "padding"
                        if( remainder > m_ConstArrayPtr->dimLast[d] ||
                            remainder < m_ConstArrayPtr->dimFirst[d] ) {
                            pos_ -= padding[d];
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
                iterator& operator++() { const_iterator::operator++(); return *this; }
                iterator& operator--() { const_iterator::operator--(); return *this; }
                // postfix
                iterator operator++( int ) { iterator i = *this; const_iterator::operator++(); return i; }
                iterator operator--( int ) { iterator i = *this; const_iterator::operator--(); return i; }

                reference operator*() { return *( data + this->pos_ ); }
                pointer operator->() { return ( data + this->pos_ ); }
                operator const_pointer() { return ( data + this->pos_ ); }

                iterator& operator+=( const const_iterator& rhs ) { this->pos_ += rhs.pos(); this->trimUp(); return *this; }
                iterator& operator-=( const const_iterator& rhs ) { this->pos_ -= rhs.pos(); this->trimUp(); return *this; }
                
                iterator& operator+=( difference_type step ) { this->pos_ += step; this->trimUp(); return *this; }
                iterator& operator-=( difference_type step ) { this->pos_ -= step; this->trimUp(); return *this; }
                
                iterator operator+( difference_type step ) const { iterator i = *this; i += step; return i; }
                iterator operator-( difference_type step ) const { iterator i = *this; i -= step; return i; }

                
            private:
                T* data;
            };

            template <typename ...S>
            Array( T* ptr, S ...sizes ) : begin_(0) {
                setSizes( sizes... );
                calcStrides();
                if( dataSize ) {
                    datablock.reset( ptr, []( T * p ) {} );
                }
            }

            /*! @name Copy constructors
             *  @note The copy constructors will share the data-block, for a "deep copy" use
             *  @code
             *  Array<double> original,deepcopy;
             *  original.copy(deepcopy);
             *  @endcode
             */ 
            //@{
            //! @brief Default, plain copy.
            Array( const Array<T>& rhs ) : dimSizes(rhs.dimSizes), dimStrides(rhs.dimStrides),
                                          dimFirst(rhs.dimFirst), dimLast(rhs.dimLast),
                                          datablock(rhs.datablock), nDims_(rhs.nDims_),
                                          nElements_(rhs.nElements_), dataSize(rhs.dataSize),
                                          dense_(rhs.dense_), begin_(rhs.begin_), end_(rhs.end_) {}
            /*! @brief The copy will be a sub-array ranging from "first" to "last"
             *  @note The dimensions of the vectors must agree with the dimensions of the array.
             */
            template <typename U>
            Array( const Array<T>& rhs, const std::vector<U>& first, const std::vector<U>& last ) : Array<T>(rhs) {
                setLimits(first,last);
            }

            /*! @brief The copy will be a sub-array.
             *  @note The dimensions of the vector must be 2x the dimensions of the array. Given as: "first1, last1, first2,...,.lastN"
             */
            template <typename U>
            Array( const Array<T>& rhs, const std::vector<U>& indices ) : Array<T>(rhs) {

                size_t nIndices = indices.size();
                if( nIndices&1 || (nIndices > 2*nDims_) )  {  // odd number of indices, or too many indices
                    throw std::logic_error("Array: Copy (sub-)constructor: Too many indices or odd number of indices: " + printArray(indices,"indices"));
                }

                if ( nIndices ) {
                    setLimits(indices);
                } else {
                    begin_ = getOffset( dimFirst );
                    auto it = const_iterator( *this, getOffset( dimLast ), begin_, getOffset( dimLast ) + nElements_ );
                    end_ = ( ++it ).pos();  // 1 step after last element;
                }
            }
            /*! @brief The copy will be a sub-array.
             *  @note The number of arguments must be 2x the dimensions of the array. Given as: "first1, last1, first2,...,.lastN"
             */
            template <typename ...S>
            Array( const Array<T>& rhs, S ...s ) : Array( rhs, std::vector<int64_t>({static_cast<int64_t>( s )...}) ) {}
            //@}
            
            /*! @name Copy/convert constructors
             *  @details Construct an array from an array of different type.
             */
            //@{
            template <typename U> Array( const Array<U>& rhs) { rhs.copy(*this); }
            template <typename U, typename V> Array( const Array<U>& rhs, const std::vector<V>& indices ) : Array<T>(rhs) { setLimits(indices); }
            template <typename U, typename ...S> Array( const Array<U>& rhs, S ...s ) : Array( rhs, std::vector<int64_t>({static_cast<int64_t>( s )...}) ) {}
            //@}
            
            /*! Construct from a generic list of dimension-sizes.
             * @note This has to be last in the header because of it's "catch all" properties.
             */
            template <typename ...S> Array( S ...sizes ) : begin_(0) { resize( sizes... ); }

            
            /*! @brief Get the @e packed size of this array
             *  @details This is only for getting the size needed to store the packed array.
             */
            size_t size( void ) const {
                size_t sz = sizeof( size_t );                                                          // blockSize
                if( dimSizes.size() ) {
                    sz += 1 + dimSizes.size() * sizeof( size_t ) + nElements_ * sizeof( T );            // nDims + dimensions + dataSize
                }
                return sz;
            }

            /*! @brief Pack the array into a dense string of characters, suitable for sending across the network or write to file.
             *  @details For sub-arrays, only the accessed region is packed, so the receiving end will unpack it into a dense array
             *  of the same dimensions as the sub-array.
             */
            uint64_t pack( char* ptr ) const {
                using redux::util::pack;
                uint64_t count = pack( ptr, size() );
                if( uint8_t ndims = dimSizes.size() ) {
                    count += pack( ptr+count, ndims );
                    count += pack( ptr+count, dimSizes );
                    count += pack( ptr+count, datablock.get(), nElements_ );
                }
                return count;

            }

            /*! @brief Unpack an array from a string of characters, swapping endianess if necessary.
             *  @details It is always the receiving/unpacking side which checks/fixes the endianess.
             */
            uint64_t unpack( const char* ptr, bool swap_endian ) {
                using redux::util::unpack;
                size_t sz;
                uint64_t count = unpack( ptr, sz, swap_endian );
                if( sz > 8 ) {
                    uint8_t ndims;
                    count += unpack( ptr+count, ndims );
                    std::vector<size_t> tmp( ndims );
                    count += unpack( ptr+count, tmp, swap_endian );
                    resize( tmp );
                    count += unpack( ptr+count, datablock.get(), nElements_, swap_endian );
                }
                return count;
            }


            /*! @name resize
             *  @brief Resize the array.
             *  @warning This will de-allocate the memory-block and allocate a new one, data is not kept/copied !!
             */
            //@{
            void resize( const std::vector<size_t>& sizes ) {
                setSizes( sizes );
                calcStrides();
                create();
            }
            template <typename ...S> void resize( S ...sizes ) { resize( {static_cast<size_t>( sizes )...} ); }
            //@}
            
            /*! @name permuteDimensions
             *  @brief Permute the dimensions specified
             *  @note This will only shuffle the dimensional information around, the datablock will not be touched.
             */
            //@{
            void permuteDimensions( const std::vector<size_t>& dims ) {
                if( dims.size() < 2 || dims.size() > nDims_ || dims[dims.size()-1] >= nDims_ ) return;
                size_t tmp = dims[dims.size()-1];
                size_t tmpSize = dimSizes[tmp];
                int64_t tmpFirst = dimFirst[tmp];
                int64_t tmpLast = dimLast[tmp];
                for(size_t i=0; i<dims.size(); ++i) {
                    if( dims[i] >= nDims_ ) throw std::out_of_range("permuteDimensions: Supplied dimension is out of range: " + printArray(dims,"dims"));
                    std::swap(dimSizes[dims[i]],tmpSize);
                    std::swap(dimFirst[dims[i]],tmpFirst);
                    std::swap(dimLast[dims[i]],tmpLast);
                }
                calcStrides();
                if( ! dense() ) {
                    begin_ = getOffset( dimFirst );
                    auto it = const_iterator( *this, getOffset( dimLast ), begin_, getOffset( dimLast ) + nElements_ );
                    end_ = ( ++it ).pos();  // 1 step after last element;
                }
            }
            
            template <typename ...S> void permuteDimensions( S ...dims ) { permuteDimensions( {static_cast<size_t>( dims )...} ); }
            //@}
            
            const std::vector<int64_t>& first( void ) const { return dimFirst; }
            const std::vector<int64_t>& last( void ) const { return dimLast; }
            const std::vector<size_t>& dimensions(void) const { return dimSizes; }
            const std::vector<size_t> dimensions( bool skipTrivialDims ) const {
                if(skipTrivialDims) {
                    std::vector<size_t> newDimSizes;
                    for( auto & it : dimSizes ) {
                        if( it > 1 || !skipTrivialDims) {
                            newDimSizes.push_back( it );
                        }
                    }
                    return std::move(newDimSizes);
                }
                return dimSizes;
            }
            size_t nDimensions( void ) const { return dimSizes.size(); }
            size_t nElements( void ) const { return nElements_; }

            size_t dimSize( size_t i = 0 ) const {
                if( i < dimSizes.size() ) {
                    return dimSizes[i];
                }
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
            Array<U> copy( bool skipTrivialDims=true ) const {
                std::vector<size_t> newDimSizes;
                for( auto & it : dimSizes ) {
                    if( it > 1 || !skipTrivialDims) {
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
            void copy( Array<U>& arr ) const {
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
            
            Array<T>& trim( bool skipTrivialDims=true ) {
                Array<T> tmp = copy(skipTrivialDims);
                *this = tmp;
                return *this;
            }

            template <typename ...S>
            T* ptr( S ...s ) {
                uint64_t offset = getOffset<uint64_t>( {static_cast<uint64_t>( s )...}, dimFirst );
                if( offset > dataSize ) {
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

            bool valid(void) const { return static_cast<bool>( datablock ); }

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
                        throw std::out_of_range( "Array::at() Index out of range: " + printArray(tmp,"indices") );
                    }
                }
                int64_t offset = getOffset( tmp, dimFirst );
                return *( datablock.get() + offset );
            }
            template <typename ...S> const T& at( S ...indices ) const { return const_cast<Array<T>*>( this )->at( indices... ); }

            std::vector<int64_t> indexFromPos(int64_t pos) const {
                std::vector<int64_t> tmp(nDims_);
                for (size_t i=0; i<nDims_; ++i) {
                    tmp[i] = pos / static_cast<int64_t>(dimStrides[i]);
                    pos -= tmp[i]*dimStrides[i];
                }
                return tmp;
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
            const Array<T>& operator+=( const Array<U>& rhs ) {
                if( this->sameSize( rhs ) ) {
                    typename Array<U>::const_iterator rhsit = rhs.begin();
                    for( auto & it : *this ) it += *rhsit++;
                }
                else {
                    throw std::invalid_argument( "image dimensions does not match." );
                }
                return *this;
            }

            template <typename U>
            const Array<T>& operator-=( const Array<U>& rhs ) {
                if( this->sameSize( rhs ) ) {
                    typename Array<U>::const_iterator rhsit = rhs.begin();
                    for( auto & it : *this ) it -= *rhsit++;
                }
                else {
                    throw std::invalid_argument( "image dimensions does not match." );
                }
                return *this;
            }

            template <typename U>
            const Array<T>& operator*=( const Array<U>& rhs ) {
                if( this->sameSize( rhs ) ) {
                    typename Array<U>::const_iterator rhsit = rhs.begin();
                    for( auto & it : *this ) it *= *rhsit++;
                }
                else {
                    throw std::invalid_argument( "image dimensions does not match." );
                }
                return *this;
            }

            template <typename U>
            const Array<T>& operator/=( const Array<U>& rhs ) {
                if( this->sameSize( rhs ) ) {
                    typename Array<U>::const_iterator rhsit = rhs.begin();
                    for( auto & it : *this ) it /= *rhsit++;
                }
                else {
                    throw std::invalid_argument( "image dimensions does not match." );
                }
                return *this;
            }

            
            template <typename U, typename V>
            const Array<T>& set( const Array<U>& rhs, V weight ) {
                if( this->sameSize( rhs ) ) {
                    typename Array<U>::const_iterator rhsit = rhs.begin();
                    for( auto & it : *this ) it = (*rhsit++ * weight);
                }
                else {
                    throw std::invalid_argument( "image dimensions does not match." );
                }
                return *this;
            }

            
            template <typename U, typename V>
            const Array<T>& add( const Array<U>& rhs, V weight ) {
                if( this->sameSize( rhs ) ) {
                    typename Array<U>::const_iterator rhsit = rhs.begin();
                    for( auto & it : *this ) it += (*rhsit++ * weight);
                }
                else {
                    throw std::invalid_argument( "image dimensions does not match." );
                }
                return *this;
            }

            
            template <typename U, typename V>
            const Array<T>& mult( const Array<U>& rhs, V weight ) {
                if( this->sameSize( rhs ) ) {
                    typename Array<U>::const_iterator rhsit = rhs.begin();
                    for( auto & it : *this ) it *= (*rhsit++ * weight);
                }
                else {
                    throw std::invalid_argument( "image dimensions does not match." );
                }
                return *this;
            }

            
            
            template <typename U>
            const Array<T>& operator=( const Array<U>& rhs ) {
                if( !sameSizes( rhs ) ) this->resize( rhs.dimSizes );
                typename Array<U>::const_iterator rhsit = rhs.begin();
                for( auto & it : *this ) it = *rhsit++;
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

            iterator mid( void ) {
                return iterator( *this, (end_-begin_)/2, begin_, end_ );
            }
            const_iterator mid( void ) const { return const_iterator( *this, (end_-begin_)/2, begin_, end_ ); }

            iterator end( void ) { return iterator( *this, end_ , begin_, end_ ); }
            const_iterator end( void ) const { return const_iterator( *this, end_ , begin_, end_ ); }


            /*!
             * @name Get a raw pointer to the datablock
             * @note Raw access only works for dense memory-blocks, so don't use it for sub-arrays !!
             * @todo Throw if non-dense
             *
             */
            //@{
            T* get( void ) { return datablock.get(); };
            const T* get( void ) const { return datablock.get(); };
            template <typename ...S>
            auto get(S ...s) -> std::shared_ptr<typename redux::util::detail::Dummy<T,S...>::dataType> {
                return redux::util::reshapeArray(datablock.get(),s...);
            }
            //@}
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
            
            template <typename U>
            void setLimits( const std::vector<U>& first, const std::vector<U>& last) {
                if( first.size() != nDims_ || last.size() != nDims_ ) {
                    throw std::logic_error("Array::setLimits: Dimensions does not match:  " + printArray(first,"first") + printArray(last,"  last"));
                }
                nElements_ = 1;
                for( size_t i = 0; i < nDims_; ++i ) {
                    int64_t tmpLast = redux::util::bound_cast<int64_t>( dimFirst[i] + last[i], 0, dimLast[i] );
                    int64_t tmpFirst = redux::util::bound_cast<int64_t>( dimFirst[i] + first[i], 0, dimLast[i] );
                    dimLast[i] = tmpLast;
                    dimFirst[i] = tmpFirst;
                    if( dimFirst[i] > dimLast[i] ) std::swap( dimFirst[i], dimLast[i] );
                    dimSizes[i] = ( dimLast[i] - dimFirst[i] + 1 );
                    if( i > 0 && ( dimStrides[i - 1] > dimSizes[i] ) ) dense_ = false;
                    nElements_ *= dimSizes[i];
                }
                begin_ = getOffset( dimFirst );
                auto it = const_iterator( *this, getOffset( dimLast ), begin_, getOffset( dimLast ) + nElements_ );
                end_ = ( ++it ).pos();  // 1 step after last element;
            }
            template <typename U>
            void setLimits( const std::vector<U>& limits ) {
                if( limits.size() != 2*nDims_ )  {  // odd number of indices, or too many indices
                    throw std::logic_error("Array::setLimits: Dimensions does not match:  " + printArray(limits,"limits"));
                }
                std::vector<U> first(nDims_);
                std::vector<U> last(nDims_);
                for( size_t i = 0; i < nDims_; ++i ) {
                    last[i] = limits[2 * i + 1];
                    first[i] = limits[2 * i];
                }
                setLimits(first,last);
            }
            template <typename ...S> void setLimits( S ...s ) { setLimits<int64_t>({static_cast<int64_t>(s)...}); }
            
            void resetLimits( void ) {
                nElements_ = 1;
                dense_ = true;
                for( size_t i = 1; i < nDims_; ++i ) {
                    dimSizes[i] = dimStrides[i-1];
                    dimFirst[i] = 0;
                    dimLast[i] = dimSizes[i]-1;
                    if( ( dimStrides[i - 1] > dimSizes[i] ) ) dense_ = false;
                    nElements_ *= dimSizes[i];
                }
                dimSizes[0] = dataSize/nElements_;
                dimFirst[0] = 0;
                dimLast[0] = dimSizes[0]-1;
                nElements_ *= dimSizes[0];
                begin_ = getOffset( dimFirst );
                auto it = const_iterator( *this, getOffset( dimLast ), begin_, getOffset( dimLast ) + nElements_ );
                end_ = ( ++it ).pos();  // 1 step after last element;
            }
            
        private:
            void setSizes( const std::vector<size_t>& sizes ) {
                begin_ = end_ = dataSize = 0;
                dimSizes = sizes;
                nDims_ = dimSizes.size();
                dimFirst.resize( nDims_, 0 );
                dimLast.assign(dimSizes.begin(), dimSizes.end());
                dense_ = true;
                if( nDims_ > 0 ) {
                    dataSize = 1;
                    for( auto & it : dimLast ) {
                        dataSize *= it;
                        it--;
                    }
                    end_ = dataSize;
                }
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
            int64_t getOffset( const std::vector<U>& indices, const std::vector<int64_t>& offsets ) const {
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

            std::vector<int64_t> dimFirst;
            std::vector<int64_t> dimLast;

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

