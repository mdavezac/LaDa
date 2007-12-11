//
//  Version: $Id$
//
#ifndef _PESCAN_H_
#define _PESCAN_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <tinyxml/tinyxml.h>

#include <pescan_interface/emass.h>
#include <pescan_interface/bandgap.h>
#include <lamarck/structure.h>
#include <opt/types.h>

#include "vff.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

namespace BandGap
{

  //! \brief Object stub which keeps track of CBM and VBM.
  //! \details This stub is supposed to be used a base class of a %GA object.
  //!          It keeps track of the conduction band minimum and of the valence
  //!          band maximum. Very few routines are implemented. Basically, it
  //!          can serialize itself for mpi purposes, it can be dumped to a
  //!          stream, and it can load/write istelf from/to and XML element. 
  //! \sa BandGap::Object, Molecularity::Object
  //! \xmlinput This stub can load \a cbm and \a vm attributes from "a" node. 
  //! \code
  //   <SOMETAG cvm="?" vbm="?"/>
  //! \endcode
  //! \a SOMETAG will generally be an %Individual.
  struct Keeper
  {
#ifdef _MPI
    //! \cond
    friend bool mpi::BroadCast::serialize<BandGap::Keeper>(BandGap::Keeper &);
    //! \endcond
#endif
    types::t_real cbm; //!< Conduction Band Minimum
    types::t_real vbm; //!< Valence Band Maximum

    //! Constructor
    Keeper() : cbm(0), vbm(0) {}
    //! Copy Constructor
    Keeper(const Keeper &_c) : cbm(_c.cbm), vbm(_c.vbm) {};
    //! Detructor
    ~Keeper() {};

    //! Loads Keeper::vbm and Keeper::cvm from attributes of \a _node.
    bool Load( const TiXmlElement &_node );
    //! Saves Keeper::vbm and Keeper::cvm as attributes of \a _node.
    bool Save( TiXmlElement &_node ) const;
  };
  //! Dumps a BandGap::Keeper to a stream.
  std::ostream& operator<<(std::ostream &_stream, const Keeper &_o);


  //! \brief GA::Evaluator stub for the band-gap (pescan) functional.
  //! \details Implements all the functionalities necessary to get and keep
  //!          pescan running. This include loading the pescan interface from
  //!          xml (see Pescan::Interface and Pescan::BandGap), reading and writing to a file
  //!          containing the reference energies, extracting the band-gap and
  //!          making sure is positive definit. If the band-gap is not positive
  //!          definit, then an all-electron calculation is performed
  //!          automatically. An all-electron calculation can also be performed
  //!          periodically, if required.
  //!          In practive an evaluator class should be derived from
  //!          GA::Evaluator and a BandGap::Darwin object used for obtaining band-gaps.
  //! \see BandGap::Evaluator, Molecularity::Evaluator
  class Darwin 
  {
    public:
      //! \brief File from which to read the atomic configuration.
      //! \details This is supposed to be given to/by vff. It should be
      //!          processor specific when executing in MPI.
      //! \see Darwin::operator<<()
      std::string atomicconfig;
      
    protected:
      //! \brief Structure on which to perform calculations
      //! \details The structure which is referenced should be updated prior to
      //!          calling Darwin::operator()(). Although much of the
      //!          information should come from Darwin::atomicconfig, such
      //!          things as the number of electrons (for all-electron
      //!          calculations) does not.
      Ising_CE::Structure &structure;
      //! The pescan interface
      Pescan::BandGap bandgap;
      //! The filename to/from are written/read the energy references for the
      //! CBM and VBM
      std::string references_filename;
      //! \brief The directory name where to perform bandgap calculations.
      //! \details Should be processor specific in the case of MPI execution.
      //!          See Pescan::Darwin::operator()() implementation for details.
      std::string dirname;
      //! Tracks the number of bandgap evaluations
      types::t_int nbeval;
      //! \brief Tracks the number of generations.
      //! \details Actually, it tracks the number of calls to
      //!          Darwin::Continue(), which should be the same thing.
      types::t_int age;
      //! How often all-electron calculations should be performed
      types::t_int check_ref_every;

    public:
      //! Constructor and Initializer
      Darwin   ( Ising_CE::Structure &_s )
             : structure(_s), references_filename("BandEdge"), 
               nbeval(0), age(0), check_ref_every(-1) {}
      //! Copy Constructor
      Darwin   ( const Darwin &_b ) 
             : structure(_b.structure), bandgap( _b.bandgap ),
               references_filename(_b.references_filename),
               nbeval(_b.nbeval), age(_b.age), check_ref_every(_b.check_ref_every) {}
      //! Destructor
      ~Darwin() {};

      //! Load Pescan::BandGap from XML
      bool Load( const TiXmlElement &_node );
      //! \brief Checks whether to perform all-electron calculations.
      //! \details if Darwin::age \% Darwin::check_ref_every is false, then the
      //! first calculations of the root node will be an all-electron on the next
      //! turrn. Darwin::check_ref_every == -1 turns this behavior off.
      bool Continue();
      //! Computed the band-gap of Darwin::structure.
      void operator()();
      //! Computes the band-gap of Darwin::structure and stores the results in
      //! \a _keeper.
      void operator()( Keeper &_keeper );
      //! \brief Correlates bandgap to  %Vff results.
      //! \details More accurately, \a _vff is asked to print atomic
      //!          configurations to Darwin::atomicconfig. That's all folks.
      //!          The templating makes this routine functional for both
      //!          Vff::Functional and its derived class Vff::Layered. In
      //!          practice, the compiler should find the type out for itself.
      template<class T_BASE>
      void operator<<( const Vff::Darwin<T_BASE> &_vff );

    protected:
      //! \brief Sets the next computation to be all-electron
      //! \details The next call to Darwin::operator()() will automatically
      //!          return the setting to folded spectra calculations.
      void set_all_electron() { bandgap.set_method( Pescan::Interface::ALL_ELECTRON ); }
      //! Reads Folded Spectra reference energies from file
      void read_references();
      //! Writes Folded Spectra reference energies to file 
      void write_references();
  };

  inline std::ostream& operator<<(std::ostream &_stream, const Keeper &_o)
  { 
    _stream << " CBM " << std::fixed << std::setw(12) << std::setprecision(6) << _o.cbm 
            << "  --  VBM " << std::fixed << std::setw(12) << std::setprecision(6) << _o.vbm; 
    return _stream; 
  } 

  inline void Darwin :: operator()( Keeper &_keeper )
  {
    Darwin::operator()();
    // copies band edges into object
    _keeper.vbm = bandgap.bands.vbm; 
    _keeper.cbm = bandgap.bands.cbm;
  }
  template <class T_BASE> 
  inline void Darwin::operator<<( const Vff::Darwin<T_BASE> &_vff )
  {
    // creates an mpi aware file name for atomic configurations
    std::ostringstream  sstr;
    sstr << "atom_config";
#ifdef _MPI
    sstr << "." << mpi::main.rank();
#endif
    // prints atomic configurations
    _vff.print_escan_input(sstr.str());
    // tells bandgap where to find atomic configurations
    atomicconfig = sstr.str();
  }


} // namespace BandGap

//! Optimizes effective elecronic mass for SL
namespace eMassSL
{
  //! \brief Object stub which keeps track of effective electronic mass
  //! \details This stub is supposed to be used a base class of a %GA object.
  //           It keeps track of the conduction band minimum and of electronic
  //           the effective mass tensor.  Very few routines are implemented.
  //           Basically, it can serialize itself for mpi purposes, it can be
  //           dumped to a stream, and it can load/write istelf from/to and XML
  //           element. 
  //! \sa BandGap::Object, Molecularity::Object
  //! \xmlinput This stub can load \a cbm and \a vm attributes from "a" node. 
  //! \code
  //   <SOMETAG cvm="?" vbm="?"/>
  //! \endcode
  //! \a SOMETAG will generally be an %Individual.
  struct Keeper
  {
#ifdef _MPI
    //! \cond
    friend bool mpi::BroadCast::serialize<eMassSL::Keeper>(eMassSL:escan::Keeper &);
    //! \endcond
#endif
    types::t_real cbm; //!< Conduction Band Minimum
    atat::rMatrix3d emass;

    //! Constructor
    Keeper() : cbm(0) { emass.zero(); }
    //! Copy Constructor
    Keeper(const Keeper &_c) : cbm(_c.cbm), emass(_c.emass) {};
    //! Detructor
    ~Keeper() {};

    //! Loads Keeper::vbm and Keeper::cvm from attributes of \a _node.
    bool Load( const TiXmlElement &_node );
    //! Saves Keeper::vbm and Keeper::cvm as attributes of \a _node.
    bool Save( TiXmlElement &_node ) const;
  };
  //! Dumps a BandGap::Keeper to a stream.
  std::ostream& operator<<(std::ostream &_stream, const Keeper &_o);


  //! \brief GA::Evaluator stub for the effective electronic mass (pescan) functional.
  //! \details Implements all the functionalities necessary to get and keep
  //!          electronic mass running. This include loading the pescan
  //!          interface from xml (see Pescan::Interface and Pescan::eMassSL),
  //!          reading and writing to a file containing the reference energy,
  //!          extracting the CBM and the effective mass. An all-electron
  //!          calculation can be performed periodically, if required.
  //!          In practive an evaluator class should be derived from
  //!          GA::Evaluator and a eMassSL::Darwin object used for obtaining
  //!          electronic effective masses.
  //! \see BandGap::Evaluator, Molecularity::Evaluator
  class Darwin
  {
    public:
      //! \brief File from which to read the atomic configuration.
      //! \details This is supposed to be given to/by vff. It should be
      //!          processor specific when executing in MPI.
      //! \see Darwin::operator<<()
      std::string atomicconfig;
      
    protected:
      //! \brief Structure on which to perform calculations
      //! \details The structure which is referenced should be updated prior to
      //!          calling Darwin::operator()(). Although much of the
      //!          information should come from Darwin::atomicconfig, such
      //!          things as the number of electrons (for all-electron
      //!          calculations) does not.
      Ising_CE::Structure &structure;
      //! The pescan interface
      Pescan::eMassSL emass;
      //! The filename to/from are written/read the energy references for the
      //! CBM and VBM
      std::string references_filename;
      //! \brief The directory name where to perform bandgap calculations.
      //! \details Should be processor specific in the case of MPI execution.
      //!          See Pescan::Darwin::operator()() implementation for details.
      std::string dirname;
      //! Tracks the number of bandgap evaluations
      types::t_int nbeval;
      //! \brief Tracks the number of generations.
      //! \details Actually, it tracks the number of calls to
      //!          Darwin::Continue(), which should be the same thing.
      types::t_int age;
      //! How often all-electron calculations should be performed
      types::t_int check_ref_every;

    public:
      //! Constructor and Initializer
      Darwin   ( Ising_CE::Structure &_s )
             : structure(_s), references_filename("BandEdge"), 
               nbeval(0), age(0), check_ref_every(-1) {}
      //! Copy Constructor
      Darwin   ( const Darwin &_b ) 
             : structure(_b.structure), emass( _b.emass ),
               references_filename(_b.references_filename),
               nbeval(_b.nbeval), age(_b.age), check_ref_every(_b.check_ref_every) {}
      //! Destructor
      ~Darwin() {};

      //! Load Pescan::eMassSL from XML
      bool Load( const TiXmlElement &_node );
      //! \brief Checks whether to perform all-electron calculations.
      //! \details if Darwin::age \% Darwin::check_ref_every is false, then the
      //! first calculations of the root node will be an all-electron on the next
      //! turrn. Darwin::check_ref_every == -1 turns this behavior off.
      bool Continue();
      //! Computed the band-gap of Darwin::structure.
      void operator()();
      //! Computes the band-gap of Darwin::structure and stores the results in
      //! \a _keeper.
      void operator()( Keeper &_keeper );
      //! \brief Correlates bandgap to  %Vff results.
      //! \details More accurately, \a _vff is asked to print atomic
      //!          configurations to Darwin::atomicconfig. That's all folks.
      //!          The templating makes this routine functional for both
      //!          Vff::Functional and its derived class Vff::Layered. In
      //!          practice, the compiler should find the type out for itself.
      template<class T_BASE>
      void operator<<( const Vff::Darwin<T_BASE> &_vff );

    protected:
      //! \brief Sets the next computation to be all-electron
      //! \details The next call to Darwin::operator()() will automatically
      //!          return the setting to folded spectra calculations.
      void set_all_electron() { emass.set_method( Pescan::Interface::ALL_ELECTRON ); }
      //! Reads Folded Spectra reference energies from file
      void read_references();
      //! Writes Folded Spectra reference energies to file 
      void write_references();
  };
}

#ifdef _MPI
namespace mpi
{
    /** \ingroup MPI
     * \brief Serializes BandGap::Keeper class for mpi purposes.
     * \details It serializes Keeper::cbm and Keeper::vbm.    */
  template<>
  inline bool BroadCast::serialize<BandGap::Keeper>( BandGap::Keeper & _keeper )
  {
    return     serialize( _keeper.cbm ) 
           and serialize( _keeper.vbm );
  }
    /** \ingroup MPI
     * \brief Serializes eMassSL::Keeper class for mpi purposes.
     * \details It serializes Keeper::cbm and Keeper::emass.    */
  template<>
  inline bool BroadCast::serialize<eMassSL::Keeper>( eMassSL::Keeper & _keeper )
  {
    return     serialize( _keeper.cbm ) 
           and serialize( _keeper.emass );
  }
}
#endif

#endif // _PESCAN_H_
