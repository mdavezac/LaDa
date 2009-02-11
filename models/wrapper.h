//
//  Version: $Id$
//
#ifndef _LADA_MODELS_WRAPPER_H_
#define _LADA_MODELS_WRAPPER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector> 

#include <boost/function.hpp> 

#include <crystal/atom.h>
#include <crystal/structure.h>
#include <opt/types.h>
#include <opt/debug.h>
#include <minimizer/cgs.h>
#include <minimizer/interpolated_gradient.h>


namespace LaDa
{
  //! Small model functionals and adapters
  namespace Models
  {
    //! The adapter for Coulomb + Lennard-Jones fortran functional itself.
    template< class T_POLICY >
      class Wrapper : public T_POLICY
      {
        public:
          //! Type of the policy.
          typedef T_POLICY t_Policy;
          //! Type of the return.
          typedef typename t_Policy :: t_Type t_Return;
          //! Type of the argument.
          typedef std::vector< t_Return > t_Arg;

          //! Constructor and Initializer
          Wrapper() : t_Policy () {}
          //! Copy Constructor
          Wrapper( const Wrapper &_c ) : t_Policy( _c ) {}
          //! \brief Destructor
          ~Wrapper() {}
  
          //! Returns the functional value.
          t_Return operator()( const t_Arg& _arg ) const;
          //! Returns the gradient.l
          t_Return gradient( const t_Arg& _arg, t_Return *const _grads ) const;

          //! Creates vectorial args from pointers.
          void create_arguments( t_Arg& _arg,
                                 size_t _Natoms,
                                 const typename t_Arg::value_type *const _cell,
                                 const typename t_Arg::value_type *const _positions ) const;
          //! Computes all at point \a _arg.
          t_Return results( const t_Arg& _arg,
                            typename t_Arg :: value_type *const _cell,
                            typename t_Arg :: value_type *const _positions,
                            t_Return *const _stress,
                            t_Return *const _forces ) const;
      };

    namespace Policy
    {
      class Dynamic
      {
        public:
          //! Initializes the dynamic functional wrapper.
          template< class T_FUNCTION > void init( const T_FUNCTION& _f ) { functional_ = _f; }

        protected: 
          //! Type of the argument and return type.
          typedef double t_Type;

          //! Constructor.
          Dynamic() {};
          //! Copy Constructor.
          Dynamic( const Dynamic& _c ) : functional_(_c.functional_) {}
          //! Destructor.
          ~Dynamic() {}

          //! Launches functional.
          t_Type operator()( const t_Type *const _cell,
                             const t_Type *const _positions,
                             t_Type *const _stress,
                             t_Type* const _forces ) const
            { return functional_( _cell, _positions, _stress, _stress ); }
          
          boost::function
          <
            double
            (
              const t_Type* const, // cell
              const t_Type* const, // positions
              t_Type* const,  // stress
              t_Type* const  // forces
            )
          > functional_;
      };
    }

    template< class T_POLICY >
      typename Wrapper<T_POLICY>::t_Return Wrapper<T_POLICY> :: operator()( const t_Arg& _arg ) const
      {
        std::vector< t_Return > grads( _arg.size(), 0 );
        const t_Return *const i_cell = &_arg[0];
        const t_Return *const i_pos = i_cell + 9;
        t_Return *const i_force = &grads[0] + 9;
        return t_Policy :: operator()( i_cell, i_pos, &grads[0], i_force );
      }
    template< class T_POLICY >
      typename Wrapper<T_POLICY> :: t_Return 
        Wrapper<T_POLICY> :: gradient( const t_Arg& _arg, t_Return *const _grads ) const
        {
          Fitting::Cgs cgs;
          cgs.verbose = false;
          cgs.tolerance = 1e-12;
          Minimizer::interpolated_gradient( *this, _arg, cgs, _grads );
          const typename t_Arg :: value_type result = operator()( _arg );
          return result;

//         const t_Return *const i_cell = &_arg[0];
//         const t_Return *const i_pos = i_cell + 9;
//         t_Return *const i_force = _grads + 9;
//         return t_Policy :: operator()( i_cell, i_pos, _grads, i_force );
        }

    template< class T_POLICY >
      void Wrapper<T_POLICY> :: create_arguments
      (
        t_Arg& _arg,
        const size_t _Natoms,
        const typename t_Arg::value_type *const _cell,
        const typename t_Arg::value_type *const _positions 
      ) const
      {
        _arg.resize( 9 + _Natoms * 3 );
        std::copy( _cell, _cell + 9, _arg.begin() );
        std::copy( _positions, _positions + _Natoms * 3, _arg.begin() + 9 );
      }

    template< class T_POLICY >
      typename Wrapper<T_POLICY> :: t_Return 
        Wrapper<T_POLICY> :: results( const t_Arg& _arg,
                                      typename t_Arg :: value_type *const _cell,
                                      typename t_Arg :: value_type *const _positions,
                                      t_Return *const _stress,
                                      t_Return *const _forces ) const
        {
          std::vector< t_Return > grads( _arg.size(), 0 );
          const t_Return result = gradient( _arg, &grads[0] );
          std::copy( _arg.begin(), _arg.begin() + 9, _cell  );
          std::copy( grads.begin(), grads.begin() + 9, _stress  );
          std::copy( _arg.begin() + 9, _arg.end(), _positions  );
          std::copy( grads.begin() + 9, grads.end(), _forces  );
          return result;
        }

  } // namespace CLJ.
} // namespace LaDa

#endif // _VFF_FUNCTIONAL_H_
