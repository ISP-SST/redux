#include "envt.hpp"

using namespace std;

template< typename T>
BaseGDL* two_fun_template( BaseGDL* p0)
{
  T* p0C = static_cast<T*>( p0);
  T* res = new T( p0C->Dim(), BaseGDL::NOZERO );
  SizeT nEl = p0->N_Elements();

  for( SizeT ii=0; ii<nEl; ++ii) {
      (*res)[ii] = ((*p0C)[ii])*2.0;
    }

  return res;
}

extern "C" BaseGDL* two_fun( EnvT* e)
{
  SizeT nParam=e->NParam();

  if (nParam != 1) {
    cout << "TWO: Improper Number of Variables" << endl;
    return new DLongGDL( -1);
  }

  BaseGDL* p0 = e->GetPar( 0);//, "TWO");
  if( p0 ) {
    if( p0->Type() == GDL_DOUBLE) {
        return two_fun_template< DDoubleGDL>( p0);
    } else if( p0->Type() == GDL_FLOAT) {
        return two_fun_template< DFloatGDL>( p0);
    } 
  }
      return new DLongGDL( -1);

}

