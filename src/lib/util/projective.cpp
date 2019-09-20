#include "redux/util/projective.hpp"


#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace redux::util;
using namespace redux;

using namespace boost::numeric::ublas;
using namespace std;


ProjectivePoints::ProjectivePoints( size_t nPoints )
    : matrix<double>(zero_matrix<double>(3, nPoints)) {

    for (size_t i=0; i<nPoints; ++ i) {
        (*this)( 2, i ) = 1;
    }
    
}


ProjectivePoints::ProjectivePoints(double x, double y, size_t nPoints)
    : matrix<double>(3, nPoints) {

    for (size_t i=0; i<nPoints; ++ i) {
        (*this)( 0, i ) = x;
        (*this)( 1, i ) = y;
        (*this)( 2, i ) = 1;
    }
    
}


void ProjectivePoints::restrict(double minValue, double maxValue) {
    
    size_t nPoints = size2();
    for( size_t i=0; i<nPoints; ++i ) {
        (*this)(0, i) = std::min(std::max((*this)(0, i), minValue), maxValue);
        (*this)(1, i) = std::min(std::max((*this)(1, i), minValue), maxValue);
    }

}

  
void ProjectivePoints::resize(size_t nPoints) {

    matrix<double>::resize(3, nPoints, true);
    for (size_t i=0; i<nPoints; ++ i) {
        (*this)( 2, i ) = 1;
    }

}


ProjectiveMap::ProjectiveMap(void) : matrix<double>(identity_matrix<double>(3)) {

}


ProjectiveMap::ProjectiveMap( const matrix<double>& rhs) {
    if(rhs.size1() != 3 || rhs.size2() != 3) {
        throw std::logic_error("Attempt to copy-construct ProjectiveMap from matrix with rows/cols != 3");
    }
    this->assign(rhs);
}


ProjectiveMap::~ProjectiveMap(void) {


}


ProjectiveMap ProjectiveMap::inverse(void) {

    ProjectiveMap tmp(*this);

    // do LU-factorization
    permutation_matrix<std::size_t> pm(3);
    if(lu_factorize(tmp, pm)) {
        throw std::logic_error("Failed to LU factorize matrix.");
    }

    ProjectiveMap inv;
    inv.assign(identity_matrix<double>(3));

    // backsubstitute
    lu_substitute(tmp, pm, inv);

    return inv;
}


void ProjectiveMap::restrict(std::vector<PointD>& pts, double minValue, double maxValue) {
    
    ProjectivePoints tmp(pts);
    
    tmp = *this * tmp;
    tmp.restrict( minValue, maxValue );
    tmp = this->inverse() * tmp;

    pts = tmp;
    
}


bool ProjectiveMap::operator==(const ProjectiveMap& rhs) {
    bool ret(true);
    for (size_t i=0; ret && i<4; ++i) {
        for (size_t j=0; ret && j<4; ++j) {
            ret &= (fabs(operator()(i,j)-rhs(i,j)) < 1E-9);
        }
    }
    return ret;
}


ProjectiveMap& ProjectiveMap::operator=(const ProjectiveMap& rhs) {
    base::operator=(rhs);
    return *this;
}


ProjectiveMap ProjectiveMap::operator*(const ProjectiveMap& rhs) const {
    ProjectiveMap ret;
    ret.assign(prod(*this,rhs));
    return ret;
}


ProjectivePoints ProjectiveMap::operator*(const ProjectivePoints& rhs) const {
    
    size_t nPoints = rhs.size2();
    ProjectivePoints ret( nPoints );
    ret.assign(prod(*this,rhs));
    for (size_t i=0; i<nPoints; ++ i) {
        double& w = ret( 2, i );
        if( w && (w != 1.0) ) {
            ret( 0, i ) /= w;
            ret( 1, i ) /= w;
            w = 1;
        }
    }
    return ret;
}


