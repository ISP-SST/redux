#ifndef REDUX_UTIL_PROJECTIVE_HPP
#define REDUX_UTIL_PROJECTIVE_HPP

#include "redux/util/point.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
namespace bnu = boost::numeric::ublas;


namespace redux {

    namespace util {


        /*!  @ingroup util
         *  @{
         */

        /*! @class     Projective projective.hpp
         *  @brief     Utilities for performing projections/maps between cameras.
         *  @author    Tomas Hillberg (hillberg@astro.su.se)
         *  @date      2016
         */
        
        class ProjectivePoints : public bnu::matrix<double> {

        public:

            explicit ProjectivePoints( size_t nPoints=1 );
            ProjectivePoints( double x, double y, size_t nPoints=1 );
            template <typename T> ProjectivePoints( const PointType<T>& pt, size_t nPoints=1 ) : matrix<double>(3, nPoints) {
                for (size_t i=0; i<nPoints; ++ i) {
                    (*this)( 0, i ) = pt.x;
                    (*this)( 1, i ) = pt.y;
                    (*this)( 2, i ) = 1;
                }
            }
            template <typename T> ProjectivePoints( const std::vector<PointType<T>>& pts ) : matrix<double>(3,1) {
                size_t nPoints = pts.size();
                if(nPoints) {
                    resize( nPoints );
                    for( size_t i=0; i<nPoints; ++i ) {
                        (*this)( 0, i ) = pts[i].x;
                        (*this)( 1, i ) = pts[i].y;
                        (*this)( 2, i ) = 1;
                    }
                }
            }
            
            void restrict(double minValue=0, double maxValue=1024);
            
            template <typename T> operator std::vector<PointType<T>>() {
                std::vector<PointType<T>> ret;
                size_t nPoints = size2();
                for( size_t i=0; i<nPoints; ++i ) {
                    ret.push_back(PointType<T>((*this)( 1, i ),(*this)( 0, i )));
                }
                return ret;
            }

            template <typename T> operator PointType<T>() {
                return PointType<T>((*this)( 1, 0 ),(*this)( 0, 0 ));
            }

            
            void resize(size_t);

            
        private:
            

        };

        
        class ProjectiveMap : public bnu::matrix<double> {
            
            typedef bnu::matrix<double> base;

        public:
            ProjectiveMap(void);
            explicit ProjectiveMap(const bnu::matrix<double>&);
            template <typename T> ProjectiveMap( const std::vector<T>& in ) : matrix<double>(bnu::identity_matrix<double>(3)) {
                size_t n = std::min<size_t>( in.size(), 9 );
                std::copy( in.begin(), in.begin()+n, data().begin() );
            }
            virtual ~ProjectiveMap(void);

            virtual ProjectiveMap inverse(void);
            
            void restrict(std::vector<PointD>& pts, double minValue=0, double maxValue=1024);

            bool operator==(const ProjectiveMap& rhs);

            ProjectiveMap& operator=(const ProjectiveMap& rhs);
            
            ProjectiveMap operator*(const ProjectiveMap& rhs) const;
            ProjectivePoints operator*(const ProjectivePoints& rhs) const;
            template <typename T> PointD operator*( const PointType<T>& rhs ) const {
                ProjectivePoints tmp(rhs);
                tmp = *this * tmp;
                return tmp;
            }
            template <typename T> std::vector<PointD> operator*( const std::vector<PointType<T>>& rhs ) const {
                ProjectivePoints tmp(rhs);
                tmp = *this * tmp;
                return tmp;
            }


        private:

            std::string m_Label;

        };
        
        

        /*! @} */


    }

}

#endif // REDUX_UTIL_PROJECTIVE_HPP
