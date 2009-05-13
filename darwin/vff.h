//
//  Version: $Id$
//
#ifndef _VFF_H_
#define _VFF_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>

#include <tinyxml/tinyxml.h>

#include <vff/functional.h>
#include <vff/va.h>
#include <crystal/structure.h>
#include <opt/function_base.h>
#include <opt/types.h>
#include <atat/vectmac.h>


namespace LaDa
{
  namespace Vff
  {
    /** \ingroup Genetic 
     * @{ */
    //! \brief Stub class for a %GA "Object" working with Vff.
    //! \details You can derive your object from this stub. It will then inherit
    //!          usefull Keeper::Load() and Keeper::Save() member routines for
    //!          reading from/writing to XML, and work well with the evaluator
    //!          stub Vff::Darwin.
    //! \xmlrestart This stub save its info, Keeper::energy and Keeper::stress, in
    //! the following format:
    //! \code
    // <VffResult energy="0.0" xx="-0.1" xy="0.0" xz="-0.1" yx="0.0" yy="-0.1" yz="0.0" zx="-0.1" zy="0.0" zz="0.0" />
    //! \endcode
    //! Yes, its ugly... but all the info is there, so stop complaining. 
    struct Keeper 
    {
      types::t_real energy; //!< The strain energy computed by %Vff.
      atat::rMatrix3d stress; //!< The stress computed by %Vff.

      //! Constructor
      Keeper() : energy(0) { stress.zero(); }
      //! Copy Constructor
      Keeper(const Keeper &_c) : energy(_c.energy), stress(_c.stress) {};
      //! Destructor
      ~Keeper() {};

      //! Loads Keeper::energy and Keeper::stress from XML.
      bool Load( const TiXmlElement &_node );
      //! Saves Keeper::energy and Keeper::stress to XML.
      bool Save( TiXmlElement &_node ) const;
      //! Serializes a vff Keeper.
      template<class Archive> void serialize(Archive & _ar, const unsigned int _version)
        { _ar & energy; _ar & stress; }
    };
    /* @} */

    /** \ingroup Genetic 
     * @{ */
    //! \brief Evaluator stub class for Valence Force Field (Vff)
    //! \details This class simply defines the necessary members to minimize a
    //!          structure using Vff. It is templated so that any Vff derivative
    //!          can be used. Whatever the flavor of Vff, the structure is
    //!          minimized using a minimizer.
    //!          Communication between an Evaluator class and this class are done
    //!          \e via the reference to the Crystal::Structure
    //!          Vff::Functional::structure which is (should) fed to the
    //!          constructor. This is not quite obvious until you see the
    //!          constructor code of BandGap::Evaluator and
    //!          Molecularity::Evaluator, so go and check it out.
    //! \param T_FUNCTIONAL Vff::Functional (default) or derived class
    template< class T_FUNCTIONAL = LaDa::Vff::Functional >
    class Darwin : public LaDa::Vff::VABase< T_FUNCTIONAL >
    {
      public:
        typedef T_FUNCTIONAL t_Functional; //!< The base class
        typedef LaDa::Vff::VABase<T_FUNCTIONAL> t_Base; //!< The base class

      public:
        //! \brief Constructor.
        //! \details  Should be used referencing a temporary instance of the
        //!           base class.
        Darwin ( const t_Base &_func ) : t_Base( _func ) {} 
        //! Copy Constructor
        Darwin ( const Darwin &_d ) : t_Base(_d) {}
        //! Destructor
        ~Darwin() {};

        //! Load t_Base and the minimizer from XML
        bool Load( const TiXmlElement &_node );
        //! Minimizes the structure
        void operator()() { t_Functional :: structure.energy = t_Base::evaluate(); }
        //! Initializes the base functional
        void init() { t_Base::init(); } 
        //! Minimizes the structure and stores the results in \a _keeper
        void operator()( Keeper &_keeper );
#       ifdef _MPI
          //! Sets communicator and suffix for mpi stuff.
          void set_mpi( boost::mpi::communicator *_comm, const std::string &_suffix )
            { t_Base::set_mpi( _comm ); }
#       endif
    };
    /* @} */



    template< class T_FUNCTIONAL>
    bool Darwin<T_FUNCTIONAL> :: Load( const TiXmlElement &_node )
    {
      if ( not t_Base::Load( _node ) )
      {
        std::cerr << " Could not load vff input!! " << std::endl; 
        return false;
      }
      if ( not t_Base::initialize_centers() )
      {
        std::cerr << " Could not initialize Atomic_Center list in vff!! " << std::endl
                  << " Are you sure the lattice and the structure correspond? " << std::endl; 
        return false;
      }

      return true;
    }
    template< class T_FUNCTIONAL>
    inline void Darwin<T_FUNCTIONAL> :: operator()( Keeper &_keeper )
    {
      operator()();
      _keeper.energy = t_Functional::structure.energy;
      _keeper.stress = t_Base::stress;
    }

  } // namespace BandGap

} // namespace LaDa
#endif // _BANDGAP_H_
