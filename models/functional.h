//
//  Version: $Id$
//
#ifndef _LADA_MODEL_FUNCTIONAL_H_
#define _LADA_MODEL_FUNCTIONAL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "clj.h"

namespace LaDa
{
  namespace Models
  {
    //! Creates a functional type interface from clj.
    class Functional : protected Clj
    {
      public:
        //! "Namespace" for relaxation type.
        struct Relaxation
        {
          enum type
          {
            default_ = 0,
            volume   = 1
          };
        };
        //! Type of the return.
        typedef types :: t_real t_Return;
        //! Type of the container.
        typedef std::vector< types :: t_real > t_Arg;
        //! Type of the gradient argument.
        typedef types :: t_real* t_GradientArg;

        //! Output structure with forces, stress, energy.
        mutable Clj::t_Arg forces;
        //! Output/input structure with atomic positions.
        mutable Clj::t_Arg structure;
        //! Type of the relaxation.
        Relaxation::type relaxation;


        
        //! \brief Constructor and Initializer
        //! \param _str structure for which to compute energy and stress
        Functional() : relaxation(Relaxation::default_) {}
        //! Constructor
        Functional   ( const Clj &_c )
                   : Clj(_c), relaxation(Relaxation::default_) {}
        //! Copy Constructor
        Functional   ( const Functional &_c )
                   : Clj(_c), forces(_c.forces), structure(_c.structure),
                     relaxation(_c.relaxation), cell0_( _c.cell0_ ), scaling_(_c.scaling_) {}
        //! \brief Destructor
        ~Functional() {}

        //! \brief Loads input to functional from  xml 
        //! \param _element should point to an xml node which is the functional data
        //! or contains a node to the funtional data as a direct child
        bool Load( const TiXmlElement &_element ) { return Clj::Load( _element ); } 

        //! Evaluates a functional after unpacking it.
        t_Return operator()( const t_Arg& _arg ) const; 

        //! Evaluates a gradient.
        void gradient( const t_Arg& _arg, t_GradientArg _grad ) const;

        //! Initializes array.
        bool init( t_Arg &_arg );


      protected:
        //! \brief unpacks variables from minimizer
        //! \details Functional knows about Functional::Structure, whereas minizers now
        //! about function::Base, this function does the interface between the two
        void unpack_variables(t_Arg const& _arg) const;
        //! \brief packs variables from minimizer
        //! \details Functional knows about Functional::Structure, whereas minizers now
        //! about function::Base, this function does the interface between the two
        void pack_variables(t_Arg& _arg) const;
        //! \brief packs variables from minimizer
        //! \details Functional knows about Functional::Structure, whereas
        //! minizers now about function::Base, this function does the interface
        //! between the two
        void pack_gradients(t_GradientArg _grad) const;
        //! Original cell, for volume relaxation.
        math::rMatrix3d cell0_;
        //! Scaling for volume relaxation.
        mutable types::t_real scaling_;
    };

  } // namespace vff 
} // namespace LaDa


#endif // _VFF_FUNCTIONAL_H_
