//
//  Version: $Id$
//
#ifndef _VFF_FUNCTIONAL_LAYERED_H_
#define _VFF_FUNCTIONAL_LAYERED_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "functional.h"

namespace Vff
{


  //! \brief Valence Force Field for "layered" structures
  //! \details In this Vff implementation, strain is only allowed in one
  //!          epitaxial growth direction, Layered::direction. In practice,
  //!          this means a few changes to variable packing and unpacking (at
  //!          least were the strain/stress is concerned), as well as
  //!          redefining member functions which make use of packing and
  //!          unpacking.  It is expected that the first unit-cell vector of
  //!          structure (from input) is the direction in which relaxation is
  //!          allowed.
  //! \see Mostly, this class is meant to work with epitaxial structure
  //!      optimization as implemented in Darwin::Molecularity. 
  class Layered : public Vff::Functional
  {
#ifdef _MPI
    /// \cond
    friend bool mpi::BroadCast::serialize<Vff::Layered> ( Vff::Layered& );
    /// \endcond
#endif
    typedef Vff::Functional t_Base; //!< Base Class
    public:
      typedef t_Base::t_Type t_Type;            //!< see Functional::Base
      typedef t_Base::t_Container t_Container;  //!< see Functional::Base
      typedef t_Container :: iterator iterator; //!< see Functional::Base
      typedef t_Container :: const_iterator const_iterator; //!< see Functional::Base

    protected:
      //! Direction in which to allow lattice-cell relaxation
      atat::rVector3d direction; 
      //! Direction in which to allow lattice-cell relaxation, normalized
      atat::rVector3d u; 
      /** \brief The Strain \f$\hat{S}\f$, as defined in Vff::Functional is
       * \f$\hat{S} = \hat{1}+\epsilon \hat{S}'\f$, with
       * \f$\hat{S}'\f$ the template strain. */
      atat::rMatrix3d  template_strain; 
      
    public:
      //! \brief Constructor and Initializer
      //! \param _str structure for which to compute energy and stress
      Layered   ( Ising_CE :: Structure &_str )
              : t_Base( _str ), direction(0,0,0), u(0,0,0) 
        { template_strain.zero(); }
      //! \brief Copy Constructor
      Layered   ( const Vff::Layered &_c )
              : t_Base( _c ), direction( _c.direction ), u(_c.u) 
        { template_strain.zero(); }
      //! \brief Destructor
      ~Layered() {}

      //! \brief unpacks function::Base::variables, then calls energy
      //! \details This function is redeclared, so that it correctly calls the
      //!          correct unpack_variables() member function from this class,
      //!          and not its base. The alternative is to make pack and unpack
      //!          virtual.
      //! \sa function::Base, function::Base::evaluate
      types::t_real evaluate(); 
      //! \brief Evaluates gradients only
      //! \sa function::Base, function::Base::evaluate_gradient
      template< typename t_grad_iterator>
        void evaluate_gradient( t_grad_iterator const &_i_grad )
          { evaluate_with_gradient( _i_grad ); }
      //! \brief Evaluates gradients only
      //! \sa function::Base, function::Base::evaluate_gradient
      void evaluate_gradient( t_Type * const _i_grad )
        { evaluate_with_gradient<t_Type*>( _i_grad ); }  
      //! \brief Evaluates gradients and energy
      //! \sa function::Base, function::Base::evaluate_with_gradient
      template< typename t_grad_iterator>
        t_Type evaluate_with_gradient( t_grad_iterator const &_i_grad );
      //! \brief Evaluates gradients and energy
      //! \sa function::Base, function::Base::evaluate_with_gradient
      t_Type evaluate_with_gradient( t_Type * const _i_grad )
        { return evaluate_with_gradient<t_Type*>( _i_grad ); }  
      //! \brief Evaluates gradient in one direction only
      //! \todo Vff::Functional::implement evaluate_one_gradient
      //! \sa function::Base, function::Base::evaluate_one_gradient, minimizer::VA
      t_Type evaluate_one_gradient( types::t_unsigned _pos) {return 0;}; 
      //! \brief initializes stuff before minimization
      //! \details Defines the packing and unpacking process, such that only unfrozen
      //! degrees of liberty are known to the minimizer
      //! \sa function::Base, minimizer::Base
      bool init();
      
    protected:
      //! \brief unpacks variables from minimizer
      //! \details Functional knows about Functional::Structure, whereas minizers now
      //! about function::Base, this function does the interface between the two
      void unpack_variables(atat::rMatrix3d& strain);
      //! \brief packs variables from minimizer
      //! \details Functional knows about Functional::Structure, whereas minizers now
      //! about function::Base, this function does the interface between the two
      void pack_variables(const atat::rMatrix3d& _strain);
      //! \brief packs variables from minimizer
      //! \details Functional knows about Functional::Structure, whereas
      //! minizers now about function::Base, this function does the interface
      //! between the two
      template< typename t_grad_iterator>
      void pack_gradients(const atat::rMatrix3d& _stress, t_grad_iterator const &_grad) const;


      //! Initializes Layered::u and Layered::template_strain
      void create_template_strain();
  };

  template< typename t_grad_iterator>
  void Vff::Layered :: pack_gradients(const atat::rMatrix3d& _stress, 
                                      t_grad_iterator const &_grad) const
  {
    t_grad_iterator i_grad(_grad);

    // first, external stuff
    *i_grad = u * ( _stress(0,0) * u );
    ++i_grad;

    // then atomic position stuff
    std::vector<Atomic_Center> :: const_iterator i_center = centers.begin();
    std::vector<Atomic_Center> :: const_iterator i_end = centers.end();
    std::vector<Ising_CE::Atom> :: const_iterator i_atom0 = structure0.atoms.begin();
    i_center = centers.begin();
    for (; i_center != i_end; ++i_center, ++i_atom0)
    {
      const atat::rVector3d& gradient = i_center->get_gradient();
      if ( not (i_atom0->freeze & Ising_CE::Atom::FREEZE_X) ) 
        *i_grad = gradient[0], ++i_grad;
      if ( not (i_atom0->freeze & Ising_CE::Atom::FREEZE_Y) ) 
        *i_grad = gradient[1], ++i_grad;
      if ( not (i_atom0->freeze & Ising_CE::Atom::FREEZE_Z) ) 
        *i_grad = gradient[2], ++i_grad;
    }
  }

  template< typename t_grad_iterator>
  types::t_real Vff::Layered :: evaluate_with_gradient( t_grad_iterator const &_i_grad )
  {
    t_Type energy = 0;
    std::for_each( centers.begin(), centers.end(),
                   std::mem_fun_ref(&Atomic_Center::reset_gradient) );

    // unpacks variables into vff atomic_center and strain format
    unpack_variables(strain);

    // computes K0
    atat::rMatrix3d K0 = (!(~strain));

    // computes energy and gradient
    std::vector<Atomic_Center> :: iterator i_center = centers.begin();
    std::vector<Atomic_Center> :: iterator i_end = centers.end();
    stress.zero();
    for (; i_center != i_end; ++i_center)
      energy += functionals[i_center->kind()].
                      evaluate_with_gradient( *i_center, strain, stress, K0 );

    // now repacks into function::Base format
    pack_gradients(stress, _i_grad);

    return energy;
  }

} // namespace vff 

#ifdef _MPI
namespace mpi {
  /** \ingroup MPI
  * \brief Serializes Vff::Layered. 
  * \details It serializes Vff::Layered::direction, as well as everything else
  *          BroadCast::serialize<Vff::Functional>() did. 
  */
  template<>
  inline bool BroadCast::serialize<Vff::Layered> ( Vff::Layered& _vff )
  {
    bool result =     BroadCast::serialize<Vff::Functional>( _vff )
                  and BroadCast::serialize( _vff.direction );

    if ( stage == COPYING_FROM_HERE ) 
      _vff.create_template_strain();

    return result;
  }
}
#endif
#endif // _VFF_FUNCTIONAL_H_
