//
//  Version: $Id$
//
#ifndef _LADA_MODELS_FORTRAN_H_
#define _LADA_MODELS_FORTRAN_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <minimizer/variant.h>
#include <boost/mpl/vector.hpp>
#include <minimizer/minuit2.h>
#include <minimizer/frprmn.h>
#include <minimizer/gsl_mins.h>

#include "wrapper.h"

namespace LaDa
{
  namespace Models
  {
    //! Relaxes the fortran lennard-jones + coulomb functional.
    class Relaxer
    {
        //! Type of the functional.
        typedef Wrapper< Policy::Dynamic > t_Wrapper;
        //! Binds the fortran subroutine.
        class Binder_;
      public:
        //! Type of the fortran functional.
        typedef void (*t_Functional)
                (
                  const int* const, // Natoms
                  const double* const,    // Cell
                  const int* const, // occupations
                  const double* const,    // positions
                  double* const,    // forces
                  double* const,    // stresses
                  double*           // energy
                ); 
        //! Internal type for return and argument.
        typedef t_Wrapper :: t_Return t_Type;

        //! Constructor.
        Relaxer() {}
        //! Copy Constructor.
        Relaxer   ( const Relaxer& _c ) 
                : minimizer_( _c.minimizer_ ),
                  wrapper_( _c.wrapper_ ),
                  functional_ ( _c.functional_ ) {}
        //! Destructor.
        ~Relaxer() {}
       
        //! Runs relaxation.
        t_Type operator()( const size_t _Natoms,
                           const int* const _occupation,
                           t_Type *const _cell, 
                           t_Type *const _positions, 
                           t_Type *const _stress, 
                           t_Type *const _forces );

        
        //! Reads minimizer info from a file.
        void read_info( const std::string &_filename );
        //! Sends function to functional.
        void init( const t_Functional &_func ) { functional_ = _func; }

      protected:
        //! Type of the minimizer.
        typedef LaDa::Minimizer::Variant
                <
                  boost::mpl::vector
                  <
                    LaDa::Minimizer::Frpr, 
                    LaDa::Minimizer::Gsl, 
                    LaDa::Minimizer::Minuit2
                  >
                > t_Minimizer;
        //! The minimizer.
        t_Minimizer minimizer_;
        //! A wrapper around the functional.
        t_Wrapper wrapper_;
        //! Pointer to the fortran subroutine itself.
        t_Functional functional_;
    };

    class Relaxer :: Binder_
    {
      public:
        Binder_   ( t_Functional &_functional,
                     const int _natoms,
                     const int* const _occupations )
                : functional_(_functional), natoms_(_natoms), occupations_(_occupations) {}
        Binder_   ( const Binder_ &_c )
                : functional_(_c.functional_), natoms_(_c.natoms_), occupations_(_c.occupations_) {}

        double operator()( const double* const _cell, const double* const _positions,
                           double* const _stress, double* const _forces ) const 
        { 
          double energy(0);
          (*functional_)( &natoms_, _cell, occupations_, _positions, _forces, _stress, &energy );
          return energy;
        }

      protected:
        //! Number of atoms;
        const int natoms_;
        //! Pointer to the occupations.
        const int* const occupations_;
        //! A pointer to the functional.
        t_Functional functional_;
    };

  } // namespace CLJ.
} // namespace LaDa


extern "C" void FC_FUNC_(create_relaxer, CREATE_RELAXER)
           ( 
             int* _handle,
             const int* _fsize,
             const char* _filename,
             LaDa :: CLJ :: Relaxer :: t_Functional
           );
extern "C" void FC_FUNC_(release_relaxer, RELEASE_RELAXER)( int* const _handle );
extern "C" void FC_FUNC_(call_relaxer, CALL_RELAXER)
           ( 
             int* const,       // Relaxer handle
             const int* const, // Natoms
             double* const,    // Cell
             const int* const, // occupations
             double* const,    // positions
             double* const,    // forces
             double* const,    // stresses
             double*           // energy
           );
extern "C" double FC_FUNC_(boost_erfc, BOOST_ERFC)( const double *const );
#endif // _VFF_FUNCTIONAL_H_

