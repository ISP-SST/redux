#ifndef REDUX_MATH_INTERVAL_HPP
#define REDUX_MATH_INTERVAL_HPP

#include "redux/math/helpers.hpp"

#include <functional>
#include <limits>


namespace redux {

    namespace math {

        /*!  @ingroup math
         *  @{
         */

        namespace detail {

            enum IntervalType { UNDEFINEDINTERVAL = 0, OPENOPEN, OPENCLOSED, CLOSEDOPEN, CLOSEDCLOSED };

        }


        /**
         *   @brief     Class representing a continuous range of values (of type T)
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2013
         */
        template< typename T, detail::IntervalType TT = detail::CLOSEDCLOSED >
        class Interval {
            
            static T minStep;

        public:
            typedef T value_type;

            Interval(T minVal = std::numeric_limits<T>::min(), T maxVal = std::numeric_limits<T>::max()) {
                set(minVal, maxVal);
            };

            T start(void) const { return m_Start; };
            T end(void) const { return m_End; };
            T mid(void) const { return (m_End+m_Start)/2; };

            void set(T minVal, T maxVal) {
                if(minVal > maxVal) {
                    std::swap(minVal, maxVal);
                }
                m_Start = minVal;
                m_End = maxVal;
            }
            
            void setEnd(T maxVal) {
                m_End = maxVal;
            }
            void setStart(T minVal) {
                m_Start = minVal;
            }

            virtual Interval& expand(const T& a, const T& b) {
                m_Start = std::min(m_Start, std::min(a, b));
                m_End = std::max(m_End, std::max(a, b));
                return *this;
            };

            virtual Interval& expand(const Interval& rhs) {
                m_Start = std::min(m_Start, rhs.m_Start);
                m_End = std::max(m_End, rhs.m_End);
                return *this;
            };
            //std::function<int(int)> 
            //Interval& trim( bool (*skip)(T), T step = minStep ) {
            virtual Interval& trim( std::function<bool(T)>& skip, T step = minStep ) {
                //std::cout << "trim    m_Min = "<< m_Min << "   m_Max = " << m_Max << "   step = " << step << std::endl;
                while( m_Start < m_End && skip(m_Start) ) {
                    //std::cout << ".";
                    m_Start += step;
                }
                while( m_Start < m_End && skip(m_End) ) {
                    //std::cout << "|";
                    m_End -= step;
                }
                //std::cout << std::endl << "trim2   m_Min = "<< m_Min << "   m_Max = " << m_Max << "   step = " << step << std::endl;
                return *this;
            }

            virtual Interval from(const T& val) const {
                if(contains(val)) {
                    return Interval<T, TT>(val, m_End);
                }
                throw std::logic_error("Value not in range.");
            };

            virtual Interval until(const T& val) const {
                if(contains(val)) {
                    return Interval<T, TT>(m_Start, val);
                }
                throw std::logic_error("Value not in range.");
            };

            virtual Interval& join(const Interval& rhs) {
                if(!intersects(rhs)) {
                    throw std::logic_error("Intervals do not overlap.");
                }
                m_Start = std::min(m_Start, rhs.m_Start);
                m_End = std::max(m_End, rhs.m_End);
                return *this;
            };

            bool intersects(const Interval& rhs) const {
                return surrounds(rhs) || inside(rhs) ||
                       ((m_Start <= rhs.m_End && m_End > rhs.m_End) || (m_End >= rhs.m_Start && m_Start < rhs.m_Start));
            };

            virtual Interval intersection(const Interval& rhs) const {
                if(surrounds(rhs)) {
                    return rhs;
                }
                else if(inside(rhs)) {
                    return *this;
                }
                else if((m_Start <= rhs.m_End && m_End > rhs.m_End) || (m_End >= rhs.m_Start && m_Start < rhs.m_Start)) {
                    return Interval<T, TT>(std::max(m_Start, rhs.m_Start), std::min(m_End, rhs.m_End));
                }
                else {
                    throw std::logic_error("Intervals do not overlap.");
                }
            };

            Interval joined(const Interval& rhs) const {
                if(!intersects(rhs)) {
                    throw std::logic_error("Intervals do not overlap.");
                }
                return Interval<T, TT>(std::min(m_Start, rhs.m_Start), std::max(m_End, rhs.m_End));
            };

            bool contains(const T& val) const {
                return (m_Start <= val && m_End >= val);
            };

            bool surrounds(const Interval& rhs) const {
                return (m_Start <= rhs.m_Start && m_End >= rhs.m_End);
            };

            bool inside(const Interval& rhs) const {
                return (m_Start >= rhs.m_Start && m_End <= rhs.m_End);
            };

            bool operator==(const Interval& rhs) const {
                return (inside(rhs) && surrounds(rhs));
                //return intersects(rhs);
            };

            bool operator!=(const Interval& rhs) const {
                return !(inside(rhs) && surrounds(rhs));
            };

            bool operator<(const Interval& rhs) const {
                return (m_Start < rhs.m_Start); // && !intersects(rhs);
            };

            bool operator>(const Interval& rhs) const {
                return (m_End > rhs.m_End); // && !intersects(rhs);
            };

        protected:

            T m_Start;
            T m_End;

        };

        template<typename T,detail::IntervalType TT>
        T Interval<T, TT>::minStep = std::numeric_limits<T>::epsilon();

        
        /*! @} */


    }
}

#endif  // REDUX_MATH_INTERVAL_HPP
