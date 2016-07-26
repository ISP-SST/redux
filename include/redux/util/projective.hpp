#ifndef REDUX_UTIL_PROJECTIVE_HPP
#define REDUX_UTIL_PROJECTIVE_HPP

#include <redux/util/point.hpp>

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

            ProjectivePoints( size_t nPoints=1 );
            ProjectivePoints( double x, double y, size_t nPoints=1 );
            ProjectivePoints( const PointD&, size_t nPoints=1 );
            ProjectivePoints( const std::vector<PointD>& );
            
            void restrict(double minValue=0, double maxValue=1024);
            
            operator std::vector<PointD>();
            
            void resize(size_t);

            
        private:
            

        };

        
        class ProjectiveMap : public bnu::matrix<double> {
            
            typedef bnu::matrix<double> base;

        public:
            ProjectiveMap(void);
            ProjectiveMap(const bnu::matrix<double>&);
            virtual ~ProjectiveMap(void);

            virtual ProjectiveMap inverse(void);
            
            void restrict(std::vector<PointD>& pts, double minValue=0, double maxValue=1024);

            bool operator==(const ProjectiveMap& rhs);

            ProjectiveMap& operator=(const ProjectiveMap& rhs);
            
            ProjectiveMap operator*(const ProjectiveMap& rhs) const;
            ProjectivePoints operator*(const ProjectivePoints& rhs) const;
            std::vector<PointD> operator*(const std::vector<PointD>& rhs) const;


        private:

            std::string m_Label;

        };
        
        

        /*! @} */


    }

}

#endif // REDUX_UTIL_PROJECTIVE_HPP
