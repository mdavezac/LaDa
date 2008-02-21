//
//  Version: $Id$
//
#ifndef _DARWIN_BANDGAP_STUBS_H_
#define _DARWIN_BANDGAP_STUBS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <tinyxml/tinyxml.h>

#include <pescan_interface/bandgap.h>
#include <lamarck/structure.h>
#include <opt/types.h>
#include <mpi/mpi_object.h>

#include "vff.h"


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
      //! \brief File in which energy references for the VBM/CBM are written and read.
      //! \details This file should contain only one line with two number: the
      //!          first number is the VBM and the second number is the CBM.
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
      //! \brief Reads energy references and checks whether to perform
      //!        all-electron calculations.
      //! \details The file Darwin::references_filename is read at the end of
      //!          each generation, whatever else happens. As such, any change
      //!          in the reference energies will only affect children in the
      //!          next generation. If the file exists from the outset \ref
      //!          TagReference will take precedence for computing the starting
      //!          population and the children of the first generation. 
      //!          Generally, this file exists only if it has been specifically
      //!          created by the user, or an all-electron calculation has been
      //!          performed. If the file exists from the outset and there are
      //!          no \ref TagReference in th input, then an all-electron
      //!          calculation is performed and the the file \e overwritten.
      //!
      //!          If Darwin::age \% Darwin::check_ref_every is false, then the
      //!          first calculations of the root node will be an all-electron
      //!          on the next turrn. Darwin::check_ref_every == -1 turns this
      //!          behavior off.
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
    sstr << "atom_config" __MPICODE(<< "." << mpi::main.rank());
    // prints atomic configurations
    _vff.print_escan_input(sstr.str());
    // tells bandgap where to find atomic configurations
    atomicconfig = sstr.str();
  }


} // namespace BandGap


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
}
#endif

#endif // _PESCAN_H_
