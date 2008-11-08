//
//  Version: $Id$
//
#ifndef _VFF_FUNCTIONAL_LAYERED_H_
#define _VFF_FUNCTIONAL_LAYERED_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "functional.h"


namespace LaDa
{ 
  namespace Vff
  {
    //! Computes in-plane stress from stress matrix \a _stress and plane \a _dir.
    types::t_real inplane_stress( const atat::rMatrix3d &_stress,
                                  const atat::rVector3d &_dir );


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
      //! The type of the atom  
      typedef Crystal::Structure::t_Atom  t_Atom;
      //! The type of the atom container
      typedef Crystal::Structure::t_Atoms t_Atoms;
      typedef Vff::Functional t_Base; //!< Base Class
      public:
        typedef t_Base::t_Type t_Type;            //!< see Functional::Base
        typedef t_Base::t_Container t_Container;  //!< see Functional::Base
        typedef t_Container :: iterator iterator; //!< see Functional::Base
        //! see Functional::Base
        typedef t_Container :: const_iterator const_iterator;

      protected:
        //! Type of the container holding the atomic centers
        typedef t_Base::t_Centers t_Centers;  
        //! Type of the atomic centers
        typedef t_Centers::value_type t_Center;  
        //! Type of the container holding the atomic functionals
        typedef t_Base::t_AtomicFunctionals t_AtomicFunctionals;  
        //! Type of the atomic functionals
        typedef t_AtomicFunctionals::value_type t_AtomicFunctional;  

      protected:
        //! Direction in which to allow lattice-cell relaxation
        atat::rVector3d direction; 
        //! Direction in which to allow lattice-cell relaxation, normalized
        atat::rVector3d u; 
        /** \brief The Strain \f$\hat{S}\f$, as defined in Vff::Functional is
         * \f$\hat{S} = \hat{1}+\epsilon \hat{S}'\f$, with
         * \f$\hat{S}'\f$ the template strain. */
        atat::rMatrix3d  template_strain; 
        //! Wether epitaxial direction if fixed by input or  structure cell
        bool is_fixed_by_input;
        
      public:
        //! \brief Constructor and Initializer
        //! \param _str structure for which to compute energy and stress
        Layered   ( Crystal :: Structure &_str )
                : Vff::Functional( _str ), direction(0,0,0), u(0,0,0),
                  template_strain(), is_fixed_by_input(false)
          { template_strain.zero(); }
        //! \brief Copy Constructor
        Layered   ( const Vff::Layered &_c )
                : t_Base( _c ), direction( _c.direction ), u(_c.u),
                  template_strain(_c.template_strain), 
                  is_fixed_by_input(_c.is_fixed_by_input) {}
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
        //! \sa function::Base, function::Base::evaluate_one_gradient, Minimizer::VA
        t_Type evaluate_one_gradient( types::t_unsigned _pos) {return 0;}; 
        //! \brief initializes stuff before minimization
        //! \details Defines the packing and unpacking process, such that only unfrozen
        //! degrees of liberty are known to the minimizer
        //! \sa function::Base, Minimizer::Base
        bool init();

        //! Prints functional to \a stream.
        void print_out( std::ostream &stream ) const;
        
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

        //! Load from XML
        bool Load( const TiXmlElement &_node );

        //! \brief Loads Functional directly from \a _node.
        //! \details If \a _node is not the correct node, the results are undefined.
        bool Load_( const TiXmlElement &_node );
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
      t_Centers :: const_iterator i_center = centers.begin();
      t_Centers :: const_iterator i_end = centers.end();
      t_Atoms :: const_iterator i_atom0 = structure0.atoms.begin();
      i_center = centers.begin();
      for (; i_center != i_end; ++i_center, ++i_atom0)
      {
        const atat::rVector3d& gradient = i_center->get_gradient();
        if ( not (i_atom0->freeze & t_Atom::FREEZE_X) ) 
          *i_grad = gradient[0], ++i_grad;
        if ( not (i_atom0->freeze & t_Atom::FREEZE_Y) ) 
          *i_grad = gradient[1], ++i_grad;
        if ( not (i_atom0->freeze & t_Atom::FREEZE_Z) ) 
          *i_grad = gradient[2], ++i_grad;
      }
    }

    template< typename t_grad_iterator>
    types::t_real Vff::Layered :: evaluate_with_gradient( t_grad_iterator const &_i_grad )
    {
      t_Type energy = 0;
      std::for_each( centers.begin(), centers.end(),
                     std::mem_fun_ref(&t_Center::reset_gradient) );

      // unpacks variables into vff atomic_center and strain format
      unpack_variables(strain);

      // computes K0
      atat::rMatrix3d K0 = (!(~strain));

      // computes energy and gradient
      t_Centers :: iterator i_center = centers.begin();
      t_Centers :: iterator i_end = centers.end();
      stress.zero();
      for (; i_center != i_end; ++i_center)
        energy += functionals[i_center->kind()].
                        evaluate_with_gradient( *i_center, strain, stress, K0 );

      // now repacks into function::Base format
      pack_gradients(stress, _i_grad);

      return energy;
    }


    inline types::t_real inplane_stress( const atat::rMatrix3d &_stress,
                                         const atat::rVector3d &_dir     )
    {
      types::t_real norm = atat::norm2(_dir);
      types::t_real trace = _stress(0,0) + _stress(1,1) + _stress(2,2);
      types::t_real axial = (_dir * (_stress * _dir) ) / norm;
      return ( trace - axial ) * 0.5;
    }

  } // namespace vff 
} // namespace LaDa

#ifdef _DOFORTRAN
  //! Creates an instance of a typical Minimizer::Frpr "C" function for calling Vff::Layered
  extern "C" inline double layeredvff_frprfun(double* _x, double* _y)
    { return Minimizer::typical_frprfun<LaDa::Vff::Layered>( _x, _y); }
  //! \brief returns a pointer to the correct extern "C" evaluation function
  //!        for Minimizer::Frpr.
  //! \details This routine allows for a standard for Vff::VA to intialize
  //!          Minimizer::Frpr.
  template<> inline t_FrprFunction choose_frpr_function<LaDa::Vff::Layered>() 
    { return layeredvff_frprfun; }
#elif defined(_DONAG)
#include <nag.h>
#include <nage04.h>
  //! Creates an instance of a typical NAG "C" function for calling Vff::Layered
  extern "C" inline void layeredvff_nagfun(int _n, double* _x, double* _r, double* _g, Nag_Comm* _p)
    { Minimizer::typical_nagfun<LaDa::Vff::Layered>( _n, _x, _r, _g, _p ); }
  //! \brief returns a pointer to the correct extern "C" evaluation function
  //!        for Minimizer::Nag.
  //! \details This routine allows for a standard for Vff::VA to intialize
  //!          Minimizer::Nag.
  template<>
  inline Minimizer::t_NagFunction choose_nag_function<LaDa::Vff::Layered>()
    { return layeredvaff_nagfun; }
#endif

#endif // _VFF_FUNCTIONAL_H_
