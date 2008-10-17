//
//  Version: $Id$
//
#ifndef _DARWIN_BANDGAP_STUBS_H_
#define _DARWIN_BANDGAP_STUBS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <tinyxml/tinyxml.h>

#include <pescan_interface/va.h>
#include <crystal/structure.h>
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
  struct Keeper: public Vff::Keeper
  {
    types::t_real cbm; //!< Conduction Band Minimum
    types::t_real vbm; //!< Valence Band Maximum

    //! Constructor
    Keeper() : cbm(0), vbm(0) {}
    //! Copy Constructor
    Keeper(const Keeper &_c) : Vff::Keeper( _c ), cbm(_c.cbm), vbm(_c.vbm) {};
    //! Detructor
    ~Keeper() {};

    //! Loads Keeper::vbm and Keeper::cvm from attributes of \a _node.
    bool Load( const TiXmlElement &_node );
    //! Saves Keeper::vbm and Keeper::cvm as attributes of \a _node.
    bool Save( TiXmlElement &_node ) const;
    //! Serializes a scalar individual.
    template<class Archive> void serialize(Archive & _ar, const unsigned int _version);
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
  class Darwin __MPICODE( : public MPI_COMMDEC )
  {
    protected:
      //! \brief Structure on which to perform calculations
      //! \details The structure which is referenced should be updated prior to
      //!          calling Darwin::operator()(). Although much of the
      //!          information should come from Darwin::atomicconfig, such
      //!          things as the number of electrons (for all-electron
      //!          calculations) does not.
      Crystal::Structure &structure;
      //! The pescan interface
      Pescan::VirtualAtom bandgap;
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
      __MPICODE
      (
        //! mpi suffix to add to calculation files.
        std::string suffix;
      )

    public:
      //! Constructor and Initializer
      Darwin   ( Crystal::Structure &_s )
             : __MPICODE( MPI_COMMCOPY( *::mpi::main ) __COMMA__ )
               structure(_s), bandgap( _s ), references_filename("BandEdge"), 
               nbeval(0), age(0), check_ref_every(-1) 
               __MPICODE( __COMMA__ suffix("") ) {}
      //! Copy Constructor
      Darwin   ( const Darwin &_b ) 
             : __MPICODE( MPI_COMMCOPY( _b ) __COMMA__ )
               structure(_b.structure), bandgap( _b.bandgap ),
               references_filename(_b.references_filename),
               nbeval(_b.nbeval), age(_b.age), check_ref_every(_b.check_ref_every) 
               __MPICODE( __COMMA__ suffix( _b.suffix ) ) {}
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
#ifdef _MPI
      //! Sets communicator and suffix for mpi stuff.
      void set_mpi( boost::mpi::communicator *_comm, const std::string &_suffix )
        { suffix = _suffix; MPI_COMMDEC::set_mpi( _comm ); bandgap.set_mpi( _comm, _suffix ); }
#endif
      //! initializes the vff part
      bool init( bool _redocenters = false )
        { return bandgap.init( _redocenters ); }

    protected:
      //! \brief Sets the next computation to be all-electron
      //! \details The next call to Darwin::operator()() will automatically
      //!          return the setting to folded spectra calculations.
      void set_all_electron()
        { bandgap.BandGap().set_method( Pescan::Interface::ALL_ELECTRON ); }
      //! Reads Folded Spectra reference energies from file
      void read_references();
      //! Writes Folded Spectra reference energies to file 
      void write_references();
  };

  inline std::ostream& operator<<(std::ostream &_stream, const Keeper &_o)
  { 
    _stream << " CBM " 
            << std::fixed << std::setw(12) << std::setprecision(6) << _o.cbm 
            << "  --  VBM "
            << std::fixed << std::setw(12) << std::setprecision(6) << _o.vbm; 
    return _stream; 
  } 

  inline void Darwin :: operator()( Keeper &_keeper )
  {
    Darwin::operator()();
    // copies band edges into object
    _keeper.vbm = bandgap.BandGap().bands.vbm; 
    _keeper.cbm = bandgap.BandGap().bands.cbm;
    _keeper.energy = structure.energy;
    bandgap.get_stress(_keeper.stress);
  }

  template<class Archive>
    void Keeper :: serialize(Archive & _ar, const unsigned int _version)
    {
      _ar & boost::serialization::base_object< Vff::Keeper >(*this); 
      _ar & vbm;
      _ar & cbm;
    }
} // namespace BandGap


#endif // _PESCAN_H_
