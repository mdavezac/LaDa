//
//  Version: $Id$
//
#ifndef _MOLECULARITY_IMPL_H_
#define _MOLECULARITY_IMPL_H_

namespace LaDa
{

  namespace Molecularity
  {
    inline void Evaluator::object_to_quantities( t_Individual & _indiv )
    {
      typedef t_Individual::t_IndivTraits::t_Object t_Object;
      const t_Object &object = (const t_Object&) _indiv;

      _indiv.quantities().clear();

      _indiv.quantities().push_back(
          std::abs(Vff::inplane_stress( object.stress, direction )) );
      _indiv.quantities().push_back( object.cbm - object.vbm );
    }
    
    inline void Evaluator :: evaluate()
    {
      Crystal::Structure copy_structure = structure;
      concentration.get( *current_object );
      current_object->x = concentration.x;
      // get band gap
#     ifdef _NOLAUNCH
        typedef t_Individual :: t_IndivTraits :: t_FourierRtoK t_Fourier;
        t_Fourier( structure.atoms.begin(), structure.atoms.end(),
                   structure.k_vecs.begin(), structure.k_vecs.end() );
#     endif
      bandgap( *current_object );

      // set quantity
      object_to_quantities( *current_individual );

      structure = copy_structure;
    }

    inline eoF<bool>* Evaluator :: LoadContinue(const TiXmlElement &_el )
    {
      return new GA::mem_zerop_t<t_Functional>( bandgap,
                                                &t_Functional::Continue,
                                                "BandGap::Continue"         );     
    }


    inline std::ostream& operator<<(std::ostream &_stream, const Object &_o)
    {
      return _stream << (const Layered::Object<>&) _o << " "
                     << (const GA::Keepers::ConcOne&)  _o << " "
                     << (const GA::Keepers::BandGap&)  _o << " ";
    }
  } // namespace Molecularity
} // namespace LaDa

#endif // _TWOSITES_IMPL_H_
