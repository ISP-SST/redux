#ifndef REDUX_UTIL_ARRAY_HPP
#define REDUX_UTIL_ARRAY_HPP

#include "redux/types.hpp"
#include "redux/math/functions.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"

#include <algorithm>
#include <stdexcept>
#include <memory>
#include <cstddef>
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
                        padding.resize(a.nDims_,a.dataSize); //a.dataSize-a.currentSizes[0]); // padding for the fastest dimension
                        for(size_t i=1; i<a.nDims_; ++i) {
                            padding[i] = a.dimStrides[i-1]-a.currentSizes[i]*a.dimStrides[i];
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
                void fastDecrement( void ) {
                    if(pos_ <= begin_) throw std::out_of_range("Array::iterator:  Iterating past begin().");
                    --pos_;
                }
                void fastIncrement( void ) { 
                    if(pos_ >= end_) {
                        std::cout << "pos=" << pos_ << "  end=" << end_ << std::endl;
                        throw std::out_of_range("Array::iterator:  Iterating past end().");
                    }
                    ++pos_;
                }

                void sparseDecrement( void ) {
                    if(pos_ <= begin_) throw std::out_of_range("Array::iterator:  Iterating past begin().");
                    --pos_;
                    trimDown();
                }

                void sparseIncrement( void ) {
                    if(pos_ >= end_) throw std::out_of_range("Array::iterator:  Iterating past end().");
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
                    pos_ = std::min(pos_,end_);
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
                    pos_ = std::max(pos_,begin_);
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
            Array( T* data, S ...sizes ) : begin_(0), end_(0) {
                wrap( data, sizes... );
            }

            /*! @name Copy constructors
             *  @note The copy constructors will share the data-block, for a "deep copy" use
             *  @code
             *  Array<double> original,deepcopy;
             *  original.copy(deepcopy);
             *  @endcode
             */ 
            //@{
            //! @brief Default, plain move/copy.
            Array( Array<T>&& rhs ) : dimSizes(std::move(rhs.dimSizes)), dimStrides(std::move(rhs.dimStrides)), nDims_(std::move(rhs.nDims_)),
                                      dataSize(std::move(rhs.dataSize)), dimFirst(std::move(rhs.dimFirst)), dimLast(std::move(rhs.dimLast)),
                                      currentSizes(std::move(rhs.currentSizes)), nElements_(std::move(rhs.nElements_)), dense_(std::move(rhs.dense_)),
                                      begin_(std::move(rhs.begin_)), end_(std::move(rhs.end_)), datablock(std::move(rhs.datablock)) { }
                                          
            Array( const Array<T>& rhs ) : dimSizes(rhs.dimSizes), dimStrides(rhs.dimStrides), nDims_(rhs.nDims_), dataSize(rhs.dataSize),
                                          dimFirst(rhs.dimFirst), dimLast(rhs.dimLast), currentSizes(rhs.currentSizes), nElements_(rhs.nElements_),
                                          dense_(rhs.dense_), begin_(rhs.begin_), end_(rhs.end_), datablock(rhs.datablock) { }
                                          
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
                    throw std::logic_error("Array: Copy (sub-)constructor: Too many indices or odd number of indices.");
                }

                if ( nIndices ) {
                    setLimits(indices);
                } else {
                    begin_ = getOffset( dimFirst );
                    auto it = const_iterator( *this, getOffset( dimLast ), begin_, getOffset( dimLast ) + nElements_ );
                    end_ = it.pos()+1;  // 1 step after last element;
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
            template <typename ...S> explicit Array( S ...sizes ) : begin_(0) { resize( sizes... ); }

            virtual ~Array() {};
            
            long int use_count(void) const { return datablock.use_count(); };
            
            inline bool empty( void ) const { return (nElements_ == 0); };
            
            /*! @brief Get the @e packed size of this array
             *  @details This is only for getting the size needed to store the packed array.
             */
            uint64_t size( void ) const {
                uint64_t sz = sizeof( uint64_t );                                                          // blockSize
                if( nElements_ ) {
                    std::vector<size_t> tmp = dimensions();
                    sz += sizeof( uint64_t ) + tmp.size() * sizeof( size_t ) + nElements_ * sizeof( T );            // nDims + dimensions + dataSize
                }
                return sz;
            }

            /*! @brief Pack the array into a dense string of characters, suitable for sending across the network or write to file.
             *  @details For sub-arrays, only the accessed region is packed, so the receiving end will unpack it into a dense array
             *  of the same dimensions as the sub-array. Trivial dimensions (size=1) will be stripped in the packing.
             */
            uint64_t pack( char* dataPtr ) const {
                using redux::util::pack;
                uint64_t count = pack( dataPtr, size() );
                if( nElements_ ) {
                    count += pack( dataPtr+count, dimensions() );
                    if(dense_) {
                        count += pack( dataPtr+count, get()+begin_, nElements_ );
                    } else {
                        T* dptr = new T[nElements_];
                        std::transform(begin(), end(), dptr, [](const T& a) { return a; });
                        //std::copy(begin(), end(), dptr);
                        count += pack( dataPtr+count, dptr, nElements_ );
                        delete[] dptr;
                    }
                }
                return count;

            }

            /*! @brief Unpack an array from a string of characters, swapping endianess if necessary.
             *  @details It is always the receiving/unpacking side which checks/fixes the endianess.
             */
            uint64_t unpack( const char* dataPtr, bool swap_endian ) {
                using redux::util::unpack;
                uint64_t sz;
                uint64_t count = unpack( dataPtr, sz, swap_endian );
                if( sz > sizeof(uint64_t) ) {          // 8 means an empty array was transferred
                    std::vector<size_t> tmp;
                    count += unpack( dataPtr+count, tmp, swap_endian );
                    resize( tmp );
                    count += unpack( dataPtr+count, datablock.get(), nElements_, swap_endian );
                } else resize();
                return count;
            }


            /*! @name resize
             *  @brief Resize the array.
             *  @warning This will de-allocate the memory-block and allocate a new one, data is not kept/copied !!
             */
            //@{
            void resize( const std::vector<size_t>& sizes ) {
                setSizes( sizes );
                setStrides();
                countElements();
                create();
            }
            template <typename ...S> void resize( S ...sizes ) { resize( {static_cast<size_t>( sizes )...} ); }
            //@}

            void clear(void) { resize(); }
            
            /*! @name wrap
             *  @brief Set the array to wrap a new raw datablock with the specified sizes
             */
            //@{
            template <typename ...S>
            void wrap( T* data, S ...sizes ) {
                setSizes( sizes... );
                setStrides();
                countElements();
                if( dataSize ) {
                    datablock.reset( data, []( T * p ) {} );
                }
                dense_ = true;
            }
            //@}

                        
            /*! @name permuteDimensions
             *  @brief Permute the dimensions specified
             *  @note This will only shuffle the dimensional information around, the datablock will not be touched.
             */
            //@{
            void permuteDimensions( const std::vector<size_t>& dims ) {
                if( dims.size() < 2 || dims.size() > nDims_ || dims[dims.size()-1] >= nDims_ ) {
                    return;
                }
                size_t tmp = dims[dims.size()-1];
                size_t tmpCurrent = currentSizes[tmp];
                size_t tmpSize = dimSizes[tmp];
                int64_t tmpFirst = dimFirst[tmp];
                int64_t tmpLast = dimLast[tmp];
                for(size_t i=0; i<dims.size(); ++i) {
                    if( dims[i] >= nDims_ ) throw std::out_of_range("permuteDimensions: Supplied dimension is out of range.");
                    std::swap(dimSizes[dims[i]],tmpSize);
                    std::swap(currentSizes[dims[i]],tmpCurrent);
                    std::swap(dimFirst[dims[i]],tmpFirst);
                    std::swap(dimLast[dims[i]],tmpLast);
                }
                setStrides();
                begin_ = getOffset( dimFirst );
                auto it = const_iterator( *this, getOffset( dimLast ), begin_, getOffset( dimLast ) + nElements_ );
                end_ = it.pos()+1;      // 1 step after last element;
                dense_ = (int64_t(begin_+nElements_) == end_);
            }
            
            template <typename ...S> void permuteDimensions( S ...dims ) { permuteDimensions( {static_cast<size_t>( dims )...} ); }
            //@}
            
            const std::vector<int64_t>& first( void ) const { return dimFirst; }
            const std::vector<int64_t>& last( void ) const { return dimLast; }
            const std::vector<size_t>& dimensions(void) const { return currentSizes; }
            const std::vector<size_t>& strides(void) const { return dimStrides; }
            const std::vector<size_t> dimensions( bool skipTrivialDims ) const {
                if(skipTrivialDims) {
                    std::vector<size_t> newDimSizes;
                    for( auto & it : currentSizes ) {
                        if( it > 1 || !skipTrivialDims) {
                            newDimSizes.push_back( it );
                        }
                    }
                    return std::move(newDimSizes);
                }
                return currentSizes;
            }
            size_t nDimensions( void ) const { return currentSizes.size(); }
            size_t nElements( void ) const { return nElements_; }

            size_t dimSize( size_t i = 0 ) const {
                if( i < currentSizes.size() ) {
                    return currentSizes[i];
                }
                else return 0;
            }

            T* cloneData( void ) {
                T* dptr = nullptr;
                if( nElements_ ) {
                    dptr = new T[nElements_];
                    if( dense_ ) {
                        std::copy( get()+begin_, get()+end_, dptr );
                    } else {
                        std::transform( begin(), end(), dptr, dptr, [](const T& a, const T& b) { return a; } );
                    }
                }
                return dptr;
            }

            template <typename U>
            void copyFrom( const void* data ) {
                const U* dptr = reinterpret_cast<const U*>( data );
                if( dense_ ) {
                    std::transform( get()+begin_, get()+end_, dptr, get()+begin_, [](const T& a, const T& b) { return static_cast<T>(b); } );
                } else {
                    std::transform( begin(), end(), dptr, begin(), [](const T& a, const T& b) { return static_cast<T>(b); }  );
                }
            }

            template <typename U>
            void copyTo( void* data ) {
                U* dptr = reinterpret_cast<U*>( data );
                if( dense_ ) {
                    std::copy( get()+begin_, get()+end_, dptr );
                    //std::transform( get()+begin_, get()+end_, dptr, get()+begin_, [](const T& a, const T& b) { return static_cast<T>(b); } );
                } else {
                    std::transform( begin(), end(), dptr, dptr, [](const T& a, const T& b) { return static_cast<U>(a); }  );
                }
            }

            template <typename U = T>
            Array<U> copy( bool skipTrivialDims=true ) const {
                std::vector<size_t> newDimSizes = dimensions(skipTrivialDims);
                if(newDimSizes.empty()) return Array<U>();
                Array<U> tmp( newDimSizes );
                tmp.assign(*this);
                return std::move(tmp);
            }


            void copy( Array<T>& out ) const {
                if( this == &out ) return;                                 // check for self-assignment
                if( (get() == out.get()) || !sameSizes( out ) ) out.resize(dimensions(true));     // if shared or wrong size: re-allocate
                if( dense_ && out.dense_ ) {
                    std::copy( get()+begin_, get()+end_, out.get()+out.begin_ );
                } else {
                    std::transform(begin(), end(), out.begin(), out.begin(), [](const T& a, const T& b) { return a; });
                }
            }

            template <typename U>
            void copy( Array<U>& out ) const {
                if( !sameSizes( out ) ) out.resize(dimensions(true));
                if( dense_ && out.dense_ ) {
                    std::copy( get()+begin_, get()+end_, out.get()+out.begin_ );
                    //std::transform(out.ptr(), out.ptr()+nElements_, get()+begin_, out.ptr(), redux::math::assign<U,T>());
                } else {
                    std::transform(out.begin(), out.end(), begin(), out.begin(), redux::math::assign<U,T>());
                }
            }



            Array<T>& trim( bool skipTrivialDims=true ) {
                Array<T> tmp = copy(skipTrivialDims);
                *this = std::move(tmp);
                return *this;
            }

            template <typename U> 
            T* ptr( const std::vector<U>& indices ) {
                int64_t offset = begin_ + getOffset( indices, dimFirst );
                if( offset < 0 || offset > dataSize ) {
                    throw std::out_of_range( "Offset out of range: " + std::to_string( offset ) );
                }
                return ( datablock.get() + offset );
            }
            template <typename U> const T* ptr( const std::vector<U>& indices ) const { return const_cast<Array<T>*>( this )->ptr( indices ); }

            template <typename ...S>
            T* ptr( S ...s ) {
                uint64_t offset = begin_ + getOffset<uint64_t>( {static_cast<uint64_t>( s )...}, dimFirst );
                if( offset > dataSize ) {
                    throw std::out_of_range( "Offset out of range: " + std::to_string( offset ) );
                }
                return ( datablock.get() + offset );
            }
            template <typename ...S> const T* ptr( S ...s ) const { return const_cast<Array<T>*>( this )->ptr( s... ); }

            bool valid(void) const { return static_cast<bool>( datablock ); }

            template <typename U> T& operator()( const std::vector<U>& indices ) { return *ptr( indices ); }
            template <typename U> const T& operator()( const std::vector<U>& indices ) const { return *const_cast<Array<T>*>( this )->ptr( indices ); }
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
                    if( tmp[i] < 0 || abs(tmp[i]) >= currentSizes[nDims_ - i - 1] ) {
                        throw std::out_of_range( "Array::at() Index out of range." );
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

            void assign( const std::vector<T>& values ) {
                if( values.size() <= nElements_ ) {
                    if(dense_) {
                        std::copy( values.begin(), values.end(), get()+begin_);
                    } else {
                        //std::copy(rhs.begin(), rhs.end(), begin());
                        std::transform(values.begin(), values.end(), begin(), begin(), [](const T& a, const T& b) { return a; });
                    }
                }
                else {
                    throw std::logic_error( "Number of total elements must match when setting a list of values." );
                }
            }

            template <typename U, typename V>
            const Array<T>& assign( const Array<U>& rhs, V weight ) {
                if( sameSize( rhs ) ) {
                   if( this != &rhs ) {
                        if(dense_ && rhs.dense_) {
                            std::transform(rhs.get()+rhs.begin_, rhs.get()+rhs.end_, get()+begin_, get()+begin_, [weight](const T& a, const T& b) { return b*weight; });
                        } else {
                            std::transform(begin(), end(), rhs.begin(), begin(), [weight](const T& a, const T& b) { return b*weight; });
                        }
                    }
                }
                else {
                    throw std::invalid_argument( "Array::assign<v>:  dimensions does not match." );
                }
                return *this;
            }

            void assign( const Array<T>& rhs ) {
                if( sameSize( rhs ) ) {
                    if( this != &rhs ) {
                        if( dense_ && rhs.dense_ ) {
                            std::copy(rhs.get()+rhs.begin_, rhs.get()+rhs.end_, get()+begin_);
                        } else {
                            //std::copy(rhs.begin(), rhs.end(), begin());
                            std::transform(begin(), end(), rhs.begin(), begin(), [](const T& a, const T& b) { return b; });
                        }
                    }
                } else {
                    throw std::invalid_argument( "Array<T>::assign:  dimensions does not match: " + printArray(dimensions(),"dims")
                                                 + printArray(rhs.dimensions(),"  rhsdims") );
                }
            }
            
            template <typename U>
            void assign( const Array<U>& rhs ) {
                if( sameSize( rhs ) ) {
                    if( dense_ && rhs.dense_ ) {
                       //std::copy(rhs.ptr(),rhs.ptr()+nElements_,get()+begin_);
                       std::transform(get()+begin_, get()+end_, rhs.get()+rhs.begin_, get()+begin_, redux::math::assign<T,U>());
                    } else {
                        std::transform(begin(), end(), rhs.begin(), begin(), redux::math::assign<T,U>());
                    }
                } else {
                    throw std::invalid_argument( "Array<U>::assign:  dimensions does not match: " + printArray(dimensions(),"dims")
                                                 + printArray(rhs.dimensions(),"  rhsdims") );
                }
            }

            
            template <typename ...S>
            void assign( S ...values ) {   assign( std::vector<T>({static_cast<T>( values )...}) ); }


            template <typename Predicate>
            void fill( T val, Predicate predicate = std::bind2nd( std::less_equal<T>(), 0 ) ) {
                if( dense_ ) {
                    std::fill(get()+begin_, get()+end_, val );
                } else {
                    //std::copy(rhs.begin(), rhs.end(), begin());
                    std::transform(begin(), end(), begin(), [val](const T& a){ return val; } );
                }
            }
            
            template <typename FillFunction, typename Predicate>
            void fill( FillFunction* filler, Predicate predicate = std::bind2nd( std::less_equal<T>(), 0 ) ) {
                if( dense_ ) {
                    std::transform(get()+begin_, get()+end_, get()+begin_, filler );
                } else {
                    //std::copy(rhs.begin(), rhs.end(), begin());
                    std::transform(begin(), end(), begin(), filler);
                }
            }
            
            /*!
             *   Assignment operators. The default assignment will make a shallow copy (shared datablock)
             */
            Array<T>& operator=( const Array<T>& rhs ) = default;

            /*!
             *  Move assignment.
             */
            Array<T>& operator=( Array<T>&& rhs )  {
                dimSizes = std::move(rhs.dimSizes);
                dimStrides = std::move(rhs.dimStrides);
                nDims_ = std::move(rhs.nDims_);
                dataSize = std::move(rhs.dataSize);
                dimFirst = std::move(rhs.dimFirst);
                dimLast = std::move(rhs.dimLast);
                currentSizes = std::move(rhs.currentSizes);
                nElements_ = std::move(rhs.nElements_);
                dense_ = std::move(rhs.dense_);
                begin_ = std::move(rhs.begin_);
                end_ = std::move(rhs.end_);
                datablock = std::move(rhs.datablock);
                return *this;
            }
            
            /*!
             *  Scalar assignment
             */
            const Array<T>& operator=( T rhs ) {
                    if( dense_ ) {
                        std::fill( get()+begin_, get()+end_, rhs );
                    } else {
                        std::fill( begin(), end(), rhs );
                    }
                return *this;
            }
            
            /*!
             *  Array assignment.
             */
            template <typename U>
            const Array<T>& operator=( const Array<U>& rhs ) {
                if( !sameSize( rhs ) ) this->resize( rhs.currentSizes );                // TBD: should the array be silently resized, or warn/throw when sizes doesn't match?
                if( dense_ && rhs.dense_ ) {
                    //std::transform(rawPtr, rawPtr+nElements_, rhs.get(), rawPtr, redux::math::assign<T,U>());
                    std::copy(rhs.get()+rhs.begin_, rhs.get()+rhs.end_, get()+begin_);
               } else {
                    std::transform(begin(), end(), rhs.begin(), begin(), redux::math::assign<T,U>());           // using the (slow) iterators defined above.
                }
                return *this;
            }
            
            
            /*!
             *  Scalar operations
             */
            Array<T> operator+( const T& rhs ) {
                Array<T> tmp;
                this->copy(tmp);
                return std::move(tmp+=rhs);
            }
            const Array<T>& operator+=( T rhs ) {
                if(dense_) {
                    std::transform(get()+begin_, get()+end_, get()+begin_, std::bind2nd<>(std::plus<T>(),rhs));
                } else {
                    std::transform(begin(), end(), begin(), std::bind2nd<>(std::plus<T>(),rhs));
                }
                return *this;
            }

            Array<T> operator-( const T& rhs ) {
                Array<T> tmp;
                this->copy(tmp);
                return std::move(tmp-=rhs);
            }
            const Array<T>& operator-=( T rhs ) {
                if(dense_) {
                    std::transform(get()+begin_, get()+end_, get()+begin_, std::bind2nd<>(std::minus<T>(),rhs));
                } else {
                    std::transform(begin(), end(), begin(), std::bind2nd<>(std::minus<T>(),rhs));
                }
                return *this;
            }
            
            Array<T> operator*( const T& rhs ) {
                Array<T> tmp;
                this->copy(tmp);
                return std::move(tmp*=rhs);
            }
            const Array<T>& operator*=( const T& rhs ) {
                if(dense_) {
                    std::transform(get()+begin_, get()+end_, get()+begin_, std::bind2nd<>(std::multiplies<T>(),rhs));
                } else {
                    std::transform(begin(), end(), begin(), std::bind2nd<>(std::multiplies<T>(),rhs));
                }
                return *this;
            }
            
            Array<T> operator/( const T& rhs ) {
                Array<T> tmp;
                this->copy(tmp);
                return std::move(tmp/=rhs);
            }
            const Array<T>& operator/=( const T& rhs ) {
                if(dense_) {
                    std::transform(get()+begin_, get()+end_, get()+begin_, std::bind2nd<>(std::divides<T>(),rhs));
                } else {
                    std::transform(begin(), end(), begin(), std::bind2nd<>(std::divides<T>(),rhs));
                }
                return *this;
            }

            
            /*!
             *  Array operations
             */
            template <typename U>
            Array<T> operator+( const Array<U>& rhs ) const {
                Array<T> tmp;
                this->copy(tmp);
                return std::move(tmp+=rhs);
            }
            template <typename U>
            const Array<T>& operator+=( const Array<U>& rhs ) {
                if( this->sameSize( rhs ) ) {
                    if(dense_ && rhs.dense_) {
                        std::transform(get()+begin_, get()+end_, rhs.get()+rhs.begin_, get()+begin_, redux::math::add<T,U>());
                    } else {
                        std::transform(begin(), end(), rhs.begin(), begin(), redux::math::add<T,U>());
                    }
                }
                else {
                    throw std::invalid_argument( "Array:+=  dimensions does not match: " + printArray(dimensions(),"dims")
                                                 + printArray(rhs.dimensions(),"  rhsdims") );
                }
                return *this;
            }
            
            template <typename U>
            Array<T> operator-( const Array<U>& rhs ) const {
                Array<T> tmp;
                this->copy(tmp);
                return std::move(tmp-=rhs);
            }
            template <typename U>
            const Array<T>& operator-=( const Array<U>& rhs ) {
                if( this->sameSize( rhs ) ) {
                    if(dense_ && rhs.dense_) {
                        std::transform(get()+begin_, get()+end_, rhs.get()+rhs.begin_, get()+begin_, redux::math::subtract<T,U>());
                    } else {
                        std::transform(begin(), end(), rhs.begin(), begin(), redux::math::subtract<T,U>());
                    }
                }
                else {
                    throw std::invalid_argument( "Array:-=  dimensions does not match: " + printArray(dimensions(),"dims")
                                                 + printArray(rhs.dimensions(),"  rhsdims") );
                }
                return *this;
            }

            template <typename U>
            Array<T> operator*( const Array<U>& rhs ) const {
                Array<T> tmp;
                this->copy(tmp);
                return std::move(tmp*=rhs);
            }
            template <typename U>
            const Array<T>& operator*=( const Array<U>& rhs ) {
                if( this->sameSize( rhs ) ) {
                    if(dense_ && rhs.dense_) {
                        std::transform(get()+begin_, get()+end_, rhs.get()+rhs.begin_, get()+begin_, redux::math::multiply<T,U>());
                    } else {
                        std::transform(begin(), end(), rhs.begin(), begin(), redux::math::multiply<T,U>());
                    }
                }
                else {
                    throw std::invalid_argument( "Array:*=  dimensions does not match: " + printArray(dimensions(),"dims")
                                                 + printArray(rhs.dimensions(),"  rhsdims") );
                }
                return *this;
            }

            template <typename U>
            Array<T> operator/( const Array<U>& rhs ) const {
                Array<T> tmp;
                this->copy(tmp);
                return std::move(tmp/=rhs);
            }
            template <typename U>
            const Array<T>& operator/=( const Array<U>& rhs ) {
                if( this->sameSize( rhs ) ) {
                    if(dense_ && rhs.dense_) {
                        std::transform(get()+begin_, get()+end_, rhs.get()+rhs.begin_, get()+begin_, redux::math::divide<T,U>());
                    } else {
                        std::transform(begin(), end(), rhs.begin(), begin(), redux::math::divide<T,U>());
                    }
                }
                else {
                    throw std::invalid_argument( "Array:/=  dimensions does not match: " + printArray(dimensions(),"dims")
                                                 + printArray(rhs.dimensions(),"  rhsdims") );
                }
                return *this;
            }
            
            
            
            virtual void zero( void ) { if( nElements_ ) memset( datablock.get(), 0, nElements_ * sizeof( T ) ); }


            template <typename U, typename V>
            const Array<T>& add( const Array<U>& rhs, const Array<V>& weight ) {
                if( this->sameSize( rhs ) && this->sameSize( weight )) {
                    typename Array<U>::const_iterator rhsit = rhs.begin();
                    typename Array<V>::const_iterator wit = weight.begin();
                    for( auto & it : *this ) it += (*rhsit++ * *wit++);
                }
                else {
                    throw std::invalid_argument( "Array::add  dimensions does not match." );
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
                    throw std::invalid_argument( "Array::add<U,V>  dimensions does not match." );
                }
                return *this;
            }

            template <typename U, typename V>
            void addTo( Array<U>& rhs, V weight ) const {
                if( this->sameSize( rhs ) ) {
                    typename Array<U>::const_iterator it = begin();
                    for( auto & rhsit : rhs ) rhsit += (*it++ * weight);
                }
                else {
                    throw std::invalid_argument( "Array::addTo  dimensions does not match." );
                }
            }

            template <typename U, typename V>
            const Array<T>& subtract( const Array<U>& rhs, V weight ) {
                if( this->sameSize( rhs ) ) {
                    if(dense_ && rhs.dense_) {
                        std::transform(get()+begin_, get()+end_, rhs.get()+rhs.begin_, get()+begin_, [weight](const T& a, const U& b) { return a-b*weight; } );
                    } else {
                        std::transform(begin(), end(), rhs.begin(), begin(), [weight](const T& a, const U& b) { return a-b*weight; } );
                    }
                }
                else {
                    throw std::invalid_argument( "Array::subtract  dimensions does not match: " + printArray(dimensions(),"dims")
                                                 + printArray(rhs.dimensions(),"  rhsdims") );
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
                    throw std::invalid_argument( "Array::mult  dimensions does not match: " + printArray(dimensions(),"dims")
                                                 + printArray(rhs.dimensions(),"  rhsdims") );
                }
                return *this;
            }

            
            bool operator==( const Array<T>& rhs ) const {
                if( ! sameSize( rhs ) ) {
                    return false;
                }
                if( get() == rhs.get() ) {      // shared data, check if it is the same sub-array.
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

            template <typename U>
            iterator pos( const std::vector<U>& indices ) {
                return iterator( *this, getOffset<U>( indices, dimFirst ), begin_, end_ );
            }
            template <typename ...S>
            iterator pos( S ...indices ) {
                return iterator( *this, getOffset<int64_t>( {static_cast<int64_t>( indices )...}, dimFirst ), begin_, end_ );
            }
            template <typename ...S>
            const_iterator pos( S ...indices ) const { return const_cast<Array<T>*>( this )->pos( indices... ); }
        //private:
            
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
        public:

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
            auto reshape(S ...s) -> std::shared_ptr<typename redux::util::detail::Dummy<T,S...>::dataType> {
                return redux::util::reshapeArray(datablock.get(),s...);
            }
            template <typename ...S>
            auto reshape(S ...s) const -> std::shared_ptr<const typename redux::util::detail::Dummy<T,S...>::dataType> {
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
                    if( ( currentSizes[i] ) != ( rhs.currentSizes[i] ) ) {
                        return false;
                    }
                }
                return true;
            }
            
            template <typename U>
            void setLimits( const std::vector<U>& first, const std::vector<U>& last) {
                if( first.size() != nDims_ || last.size() != nDims_ ) {
                    throw std::logic_error("Array::setLimits: Dimensions does not match." );
                }
                nElements_ = 1;
                for( size_t i = 0; i < nDims_; ++i ) {
                    //int64_t tmpLast = redux::util::bound_cast<int64_t>( dimFirst[i] + last[i], 0, dimLast[i] );
                    //int64_t tmpFirst = redux::util::bound_cast<int64_t>( dimFirst[i] + first[i], 0, dimLast[i] );
                    int64_t tmpLast  = std::min( std::max(int64_t(dimFirst[i] + last[i]), 0L), dimLast[i] );
                    int64_t tmpFirst = std::min( std::max(int64_t(dimFirst[i] + first[i]), 0L), dimLast[i] );
                    dimLast[i] = tmpLast;
                    dimFirst[i] = tmpFirst;
                    if( dimFirst[i] > dimLast[i] ) std::swap( dimFirst[i], dimLast[i] );
                    currentSizes[i] = ( dimLast[i] - dimFirst[i] + 1 );
                   // if( i > 0 && ( dimStrides[i - 1] > currentSizes[i] ) ) dense_ = false;
                    nElements_ *= currentSizes[i];
                }
                begin_ = getOffset( dimFirst );
                auto it = const_iterator( *this, getOffset( dimLast ), begin_, getOffset( dimLast ) );
                end_ = it.pos()+1;  // 1 step after last element;
                dense_ = (int64_t(begin_+nElements_) == end_);
                //std::cout << "Arr::SetLimits:  begin=" << begin_ << "  end=" << end_ << "  nEl=" << nElements_
                //<< " lastOff=" << getOffset( dimLast ) << " dense=" << dense_ << std::endl;
            }
            
            template <typename U>
            void setLimits( const std::vector<U>& limits ) {
                std::vector<U> first(nDims_);
                std::vector<U> last(nDims_);
                if( limits.size() == nDims_ ) {
                    for( size_t i = 0; i < nDims_; ++i ) {
                        last[i] = dimFirst[i]+limits[i]-1;
                        first[i] = dimFirst[i];
                    }
                } else if (limits.size() == 2*nDims_) {
                    for( size_t i = 0; i < nDims_; ++i ) {
                        last[i] = limits[2 * i + 1];
                        first[i] = limits[2 * i];
                    }
                } else {  // odd number of indices, or too many indices
                    throw std::logic_error("Array::setLimits: Dimensions does not match.");
                }
                setLimits(first,last);
            }
            template <typename ...S> void setLimits( S ...s ) { setLimits<int64_t>({static_cast<int64_t>(s)...}); }
            
            void resetLimits( void ) {
                nElements_ = 1;
                for( size_t i = 1; i < nDims_; ++i ) {
                    currentSizes[i] = dimStrides[i-1];
                    dimFirst[i] = 0;
                    dimLast[i] = currentSizes[i]-1;
                    nElements_ *= currentSizes[i];
                }
                currentSizes[0] = dataSize/nElements_;
                dimFirst[0] = 0;
                dimLast[0] = currentSizes[0]-1;
                nElements_ *= currentSizes[0];
                begin_ = getOffset( dimFirst );
                auto it = const_iterator( *this, getOffset( dimLast ), begin_, getOffset( dimLast ) + nElements_ );
                end_ = it.pos()+1;  // 1 step after last element;
                dense_ = (int64_t(begin_+nElements_) == end_);
            }
            

            int shift( size_t dim, int n ) {

                if( dim >= nDims_ || n == 0 ) {
                    return 0;   // no shift performed
                }
                
              
                // if we would step out of bounds, go as far as possible and return the actual shift.
                if( n+dimFirst[dim] < 0 ) {
                    n = -dimFirst[dim];
                }
                
                if(abs(n+dimLast[dim]) >= dimSizes[dim]) {
                    n = dimSizes[dim]-dimLast[dim]-1;
                }
              
                dimFirst[dim] += n;
                dimLast[dim] += n;
                
                begin_ = getOffset( dimFirst );
                auto it = const_iterator( *this, getOffset( dimLast ), begin_, getOffset( dimLast ) + nElements_ );
                end_ = it.pos()+1;  // 1 step after last element;
                //dense_ = (int64_t(begin_+nElements_) == end_);         // should not change with shift, but it's a cheap calculation anyway...
               
                return n;   // return 
            }

            template <typename U>
            int64_t getOffset( const std::vector<U>& indices, const std::vector<int64_t>& offsets ) const {
                int64_t offset = 0;
                if( indices.size() > offsets.size() ) {
                    return getOffset( indices );
                }
                else if( indices.size() > currentSizes.size() ) {
                    throw std::logic_error( "Array::getOffset() called with more indices than dimensions." );
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
                if( indices.size() > currentSizes.size() ) {
                    throw std::logic_error( "Array::getOffset() called with more indices than dimensions." );
                }
                else {
                    size_t dimDiff = nDims_ - indices.size();
                    for( size_t i = 0; i < indices.size(); ++i ) {
                        offset += indices[i] * dimStrides[dimDiff + i];
                    }
                }
                return offset;
            }
            
            std::shared_ptr<T>& getData(void) { return datablock; }

        private:
            void setSizes( const std::vector<size_t>& sizes ) {
                begin_ = end_ = nElements_ = dataSize = 0;
                dimSizes = sizes;
                dimSizes.erase(std::remove(dimSizes.begin(), dimSizes.end(), 0), dimSizes.end());
                nDims_ = dimSizes.size();
                currentSizes = dimSizes;
                dimFirst.resize( nDims_ );
                if( nDims_ > 0 ) {
                    dimLast.assign(dimSizes.begin(), dimSizes.end());
                    dense_ = true;
                    dataSize = 1;
                    for( auto & it : dimLast ) {
                        dataSize *= it;
                        it--;
                    }
                    end_ = nElements_ = dataSize;
                }
            }
            template <typename ...S> void setSizes( S ...sizes ) { setSizes( {static_cast<size_t>( sizes )...} );  }

            void create( void ) {
                if( dataSize ) {
                    datablock = sharedArray<T>(dataSize);
                }
                else {
                    datablock.reset();
                }
            }
            
            void setStrides(void) {
                if( nDims_ > 0 ) {
                    dimStrides.resize( nDims_ );
                    dimStrides.back() = 1;
                    for( int i = nDims_-1; i > 0; --i ) {
                        dimStrides[i - 1] = dimSizes[i] * dimStrides[i];
                    }
                }
            }

            void countElements( void ) {
                nElements_ = 0;
                if( currentSizes.size() > 0 ) {
                    nElements_ = 1;
                    for( int i = currentSizes.size() - 1; i > 0; --i ) {
                        nElements_ *= ( dimLast[i] - dimFirst[i] + 1 );
                    }
                    nElements_ *= ( dimLast[0] - dimFirst[0] + 1 );
                }
            }
            
            // These are modified on construct/resize and then constant
            std::vector<size_t> dimSizes;
            std::vector<size_t> dimStrides;
            size_t nDims_;
            size_t dataSize;

            // These are modified with sub-array operations
            std::vector<int64_t> dimFirst;
            std::vector<int64_t> dimLast;
            std::vector<size_t> currentSizes;
            size_t nElements_;
            bool dense_;
            int64_t begin_, end_;
            
            
        protected:
            std::shared_ptr<T> datablock;

            template<typename U> friend class Array;
            friend class const_iterator;
            friend class iterator;

        };
        
        template <> template <typename U>
        void Array<complex_t>::copy( Array<U>& out ) const {
            if( !sameSizes( out ) ) out.resize( currentSizes );
            const_iterator cit = begin();
            for( auto & it : out ) it = (*cit++).real();
        }
    
        template <> template <>
        inline const Array<double>& Array<double>::operator=( const Array<redux::complex_t>& rhs ) {
            if( !sameSize( rhs ) ) this->resize( rhs.currentSizes );
            if(dense_ && rhs.dense_) {
                std::transform(rhs.get()+rhs.begin_, rhs.get()+rhs.end_, get()+begin_, [](complex_t d) { return std::real(d); } );
            } else {
                std::transform(rhs.begin(), rhs.end(), begin(), [](complex_t d) { return std::real(d); } );
            }
            return *this;
        }

        template <> template <>
        inline const Array<float>& Array<float>::operator=( const Array<redux::complex_t>& rhs ) {
            if( !sameSize( rhs ) ) this->resize( rhs.currentSizes );
            if(dense_ && rhs.dense_) {
                std::transform(rhs.get()+rhs.begin_, rhs.get()+rhs.end_, get()+begin_, [](complex_t d) { return std::real(d); } );
            } else {
                std::transform(rhs.begin(), rhs.end(), begin(), [](complex_t d) { return std::real(d); } );
            }
            return *this;
        }
        

        /*! @} */


    }

}

#endif // REDUX_UTIL_ARRAY_HPP

