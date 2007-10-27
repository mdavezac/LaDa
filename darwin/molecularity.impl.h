//
//  Version: $Id$
//
#ifndef _MOLECULARITY_IMPL_H_
#define _MOLECULARITY_IMPL_H_

#include <algorithm>
#include <functional>
#include <ext/algorithm>
#ifndef __PGI
  #include<ext/functional>
  using __gnu_cxx::compose1;
#else
  #include<functional>
  using std::compose1;
#endif

#include "print/xmg.h"
#include "print/stdout.h"
#include "functors.h"

namespace Molecularity
{
  inline void Evaluator::object_to_quantities( t_Individual & _indiv )
  {
    typedef t_Individual::t_IndivTraits::t_Object t_Object;
    const t_Object &object = (const t_Object&) _indiv;

    _indiv.quantities().clear();

    _indiv.quantities().push_back( inplane_stress( object.stress, direction ) );
    _indiv.quantities().push_back( object.cbm - object.vbm );
  }
  
  inline void Evaluator :: evaluate()
  {
    // relax structure
    vff();
    // Load relaxed structure into pescan
    pescan << vff; 

#ifdef _NOLAUNCH
    typedef t_Individual :: t_IndivTraits :: t_FourierRtoK t_Fourier;
    Layered::Fourier<2>( structure.atoms.begin(), structure.atoms.end(),
                         structure.k_vecs.begin(), structure.k_vecs.end() );
#endif

    // get band gap
    pescan( *current_object );
  
    // set quantity
    object_to_quantities( *current_individual );
  }

  inline eoF<bool>* Evaluator :: LoadContinue(const TiXmlElement &_el )
  {
    return new GA::mem_zerop_t<Pescan::Darwin>( pescan,
                                                &Pescan::Darwin::Continue,
                                                "Pescan::Continue"         );     
  }

  inline types::t_real inplane_stress( const atat::rMatrix3d &_stress,
                                       const atat::rVector3d &_dir     )
  {
    types::t_real norm = atat::norm2(_dir);
    types::t_real trace = _stress(0,0) + _stress(1,1) + _stress(2,2);
    atat::rVector3d axial = _stress * _dir;
    return ( trace - (   axial(0) * _dir(0)
                       + axial(1) * _dir(1)
                       + axial(2) * _dir(2) ) / norm ) * 0.5;
  }



  inline std::ostream& operator<<(std::ostream &_stream, const Object &_o)
  {
    return _stream << (const Layered::Object<>&) _o << " " 
                   << (const Pescan::Keeper&)  _o << " ";
  }
} // namespace Molecularity


#ifdef _MPI
namespace mpi
{
  template<>
  inline bool mpi::BroadCast::serialize<Molecularity::Object>
                                       ( Molecularity::Object & _object )
  {
    return     serialize<Pescan::Keeper>( _object )
           and serialize<Vff::Keeper>( _object )
           and _object.broadcast( *this );
  }
}
#endif

#endif // _TWOSITES_IMPL_H_
