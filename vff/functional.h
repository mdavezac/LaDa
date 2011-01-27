#ifndef _LADA_VFF_ORIGINAL_H_
#define _LADA_VFF_ORIGINAL_H_

#include "LaDaConfig.h"

#include "vff.h"

namespace LaDa
{
  //! \brief Reimplements the Valence Force Field %Functional in c++
  //! \details Vff, or Valence Force Field functional, is an empirical functional which
  //! attempts to model the strain energy of a material from at most three-body
  //! interactions. The there body interactions which are considered are
  //! bond-stretching, change in bond-angles, and a combination of these two.
  //! 
  //! The implementation relies on a body-centered paradigm. In other words,
  //! four classes have been created:
  //!   - Vff::Functional is wrapper class and interface to Vff
  //!   - Vff::AtomicCenter represent a single atom and lists its first neighbor relationships
  //!   - Vff::AtomicCenter::const_iterator allows coders to travel
  //!   through a collection of Vff::Atomic_Centera along first neighbor
  //!   relationships.
  //!   - Vff::AtomicFunctional computes the strain energy of one single atom, eg
  //!   all three body terms in which a particular atom takes part.
  //!   .
  namespace vff
  {
    //! Vff functional called using a string of numbers rather than a complete structure.
    class Functional : public Vff
    {
      public:
        //! Type of the return.
        typedef types :: t_real t_Return;
        //! Type of the container.
        typedef std::vector< types :: t_real > t_Arg;
        //! Type of the gradient argument.
        typedef types :: t_real* t_GradientArg;

        
      public:
        //! Constructor and Initializer
        Functional() : fixed_index(-1,-1,-1)
          { stress = math::rMatrix3d::Zero(); }
        //! \brief Constructor and Initializer
        //! \param _str structure for which to compute energy and stress
        Functional   (Vff::t_Arg const &_str)
                   : Vff( _str ), structure0(_str), fixed_index(-1,-1,-1)
          { stress = math::rMatrix3d::Zero(); }
        //! \brief Copy Constructor
        Functional   ( const Functional &_c )
                   : Vff( _c ), structure0( _c.structure0 ), stress( _c.stress ),
                     fixed_index( _c.fixed_index ) {}
        //! \brief Destructor
        ~Functional() {}

#       ifdef LADA_MPI
          //! Evaluates a functional after unpacking it.
          t_Return operator()( const t_Arg& _arg, boost::mpi::communicator const &_comm)
            { comm = _comm; return operator()(_arg); }
  
          //! Evaluates a gradient.
          t_Return gradient( const t_Arg& _arg, t_GradientArg _grad,
                             boost::mpi::communicator const &_comm)
            { comm = _comm; return gradient(_arg, _grad); }
#       endif
        //! Evaluates a functional after unpacking it.
        t_Return operator()(const t_Arg& _arg) const; 

        //! Evaluates a gradient.
        t_Return gradient(const t_Arg& _arg, t_GradientArg _grad) const;

        //! \brief initializes stuff before minimization
        //! \details Defines the packing and unpacking process, such that only unfrozen
        //! degrees of liberty are known to the minimizer
        //! \sa function::Base, Minimizer::Base
        bool init( t_Arg& _arg );
        //! \brief Returns a reference to the computed stress
        //! \sa Functional::stress
        const math::rMatrix3d& get_stress() const { return stress; }

        //! Copies parameters from argument.
        void copy_parameters(Functional const &_f)
        {
          bond_cutoff = _f.bond_cutoff;
          bonds_params = _f.bonds_params;
          angles_params = _f.angles_params;
        }

      protected:
        //! \brief unpacks variables from minimizer
        //! \details Functional knows about Functional::Structure, whereas minizers now
        //! about function::Base, this function does the interface between the two
        void unpack_variables(const t_Arg& _arg, math::rMatrix3d& strain) const;
        //! Unpacks position variables only.
        void unpack_positions( t_Arg::const_iterator& _i_x, math::rMatrix3d& strain) const;
        //! \brief packs variables from minimizer
        //! \details Functional knows about Functional::Structure, whereas minizers now
        //! about function::Base, this function does the interface between the two
        void pack_variables( t_Arg& _arg, const math::rMatrix3d& _strain ) const;
        //! Packs position variables only.
        void pack_positions( t_Arg :: iterator & _i_x) const;
        //! Counts positional degrees of freedom.
        types::t_unsigned posdofs();
        //! \brief packs variables from minimizer
        //! \details Functional knows about Functional::Structure, whereas
        //! minizers now about function::Base, this function does the interface
        //! between the two
        void pack_gradients(const math::rMatrix3d& _stress, t_GradientArg _grad) const;
      
        //! original structure,  needed for gradients
        Vff::t_Arg structure0;
        //! stores stress in Functional::structure after computation
        mutable math::rMatrix3d stress;
        //! Index of the first atoms with fixed x, y, z;
        math::iVector3d fixed_index; 
    };

  } // namespace vff 
} // namespace LaDa


#endif // _VFF_FUNCTIONAL_H_
