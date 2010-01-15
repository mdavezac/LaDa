//
//  Version: $Id$
//
#ifndef _DARWIN_EMASS_IMPL_H_
#define _DARWIN_EMASS_IMPL_H_

#include <lapack/lapack.h>

namespace LaDa
{
  namespace eMassSL
  {
    inline bool Object::Load( const TiXmlElement &_node )
    { 
      if( not t_VffBase::Load(_node) ) return false;
      if( not t_eMassSLBase::Load(_node) )  return false;

      types::t_real eigs[3];
      Eigen::Matrix3d vecs;
      Lapack::eigen( emass, vecs, eigs );

      n = 2;
      Q = eigs[0] * eigs[1] / eigs[2] ; 
      types::t_real q =  eigs[1] * eigs[2] / eigs[0] ; 
      if( q > Q ) { Q = q; n = 0; }
      q =  eigs[2] * eigs[0] / eigs[1] ; 
      if( q > Q ) { Q = q; n = 1; }

      return true;
    }
    inline bool Object::Save( TiXmlElement &_node ) const
    {
      if( not t_VffBase::Save(_node) ) return false;
      if( not t_eMassSLBase::Save(_node) )  return false;

      types::t_real eigs[3];
      Eigen::Matrix3d vecs_;
      Lapack::eigen( emass, vecs_, eigs );

      std::ostringstream dir; 
      dir << vecs_(0,n) << " " << vecs_(1,n) << " "  << vecs_(2,n);
      _node.SetAttribute( "direction", dir.str().c_str() );

      return true;
    }
    inline void Evaluator::object_to_quantities( t_Individual & _indiv )
    {
      typedef t_Individual::t_IndivTraits::t_Object t_Object;
      const t_Object &object = (const t_Object&) _indiv;

      _indiv.quantities().clear();

      _indiv.quantities().push_back( Vff::inplane_stress( object.stress, direction ) );

      _indiv.quantities().push_back( object.Q );
    }
    
    inline void Evaluator :: evaluate()
    {
      Crystal :: Structure structure0 = structure;
      // relax structure
      vff( *current_object );
      // Load relaxed structure into emass
      emass << vff; 
      structure = structure0;

      // get effective mass
      emass( *current_object );
    
      // get best Q
      types::t_real eigs[3];
      Eigen::Matrix3d vecs;
      Lapack::eigen( current_object->emass, vecs, eigs );

      current_object->n = 2;
      current_object->Q = eigs[0] * eigs[1] / eigs[2] ; 
      types::t_real q = eigs[1] * eigs[2] / eigs[0] ; 
      if( q > current_object->Q ) { current_object->Q = q; current_object->n = 0; }
      q = eigs[2] * eigs[0] / eigs[1] ; 
      if( q > current_object->Q ) { current_object->Q = q; current_object->n = 1; }

      // set quantity
      object_to_quantities( *current_individual );
    }

    inline eoF<bool>* Evaluator :: LoadContinue(const TiXmlElement &_el )
    {
      return new GA::mem_zerop_t<eMassSL::Darwin>( emass,
                                                  &eMassSL::Darwin::Continue,
                                                  "eMassSL::Darwin::Continue"         );     
    }

    inline types::t_real inplane_stress( const Eigen::Matrix3d &_stress,
                                         const Eigen::Vector3d &_dir     )
    {
      types::t_real norm = _dir->squaredNorm();
      types::t_real trace = _stress(0,0) + _stress(1,1) + _stress(2,2);
      types::t_real axial = (_dir * (_stress * _dir) ) / norm;
      return ( trace - axial ) * 0.5;
    }



    inline std::ostream& operator<<(std::ostream &_stream, const Object &_o)
    {

      types::t_real eigs[3];
      Eigen::Matrix3d vecs;
      Lapack::eigen( _o.emass, vecs, eigs );

      return _stream << (const Layered::Object<>&) _o << " "
                     << " Q=" << _o.Q
                     << " dir: " << vecs(0,_o.n) << " " << vecs(1,_o.n) << " " << vecs(2,_o.n);
    }
  } // namespace Molecularity
} // namespace LaDa


#endif // _TWOSITES_IMPL_H_
