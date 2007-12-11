//
//  Version: $Id$
//
#ifndef _MOLECULARITY_IMPL_H_
#define _MOLECULARITY_IMPL_H_

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
    Ising_CE :: Structure structure0 = structure;
    // relax structure
    vff( *current_object );
    // Load relaxed structure into bandgap
    bandgap << vff; 
    structure = structure0;

    // get band gap
    bandgap( *current_object );
  
    // set quantity
    object_to_quantities( *current_individual );
  }

  inline eoF<bool>* Evaluator :: LoadContinue(const TiXmlElement &_el )
  {
    return new GA::mem_zerop_t<BandGap::Darwin>( bandgap,
                                                &BandGap::Darwin::Continue,
                                                "BandGap::Continue"         );     
  }

  inline types::t_real inplane_stress( const atat::rMatrix3d &_stress,
                                       const atat::rVector3d &_dir     )
  {
    types::t_real norm = atat::norm2(_dir);
    types::t_real trace = _stress(0,0) + _stress(1,1) + _stress(2,2);
    types::t_real axial = (_dir * (_stress * _dir) ) / norm;
    return ( trace - axial ) * 0.5;
  }



  inline std::ostream& operator<<(std::ostream &_stream, const Object &_o)
  {
    return _stream << (const Layered::Object<>&) _o << " "
                   << (const BandGap::Keeper&)  _o << " ";
  }
} // namespace Molecularity


#ifdef _MPI
namespace mpi
{
  template<>
  inline bool mpi::BroadCast::serialize<Molecularity::Object>
                                       ( Molecularity::Object & _object )
  {
    return     serialize<BandGap::Keeper>( _object )
           and serialize<Vff::Keeper>( _object )
           and _object.broadcast( *this );
  }
}
#endif

#endif // _TWOSITES_IMPL_H_
