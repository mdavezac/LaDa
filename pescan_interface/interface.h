//
//  Version: $Id$
//
#ifndef _PESCAN_INTERFACE_H_
#define _PESCAN_INTERFACE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>

#include <string>
#include <vector>
#include <utility>
#include <fstream>

#include <tinyxml/tinyxml.h>

#include <physics/physics.h>
#include <opt/types.h>
#include <opt/tuple_serialize.h>
#include <crystal/structure.h>
#include <mpi/mpi_object.h>
#include <print/stdout.h>
#include <math/serialize.h>

#ifdef _DIRECTIAGA
# define __DIAGA( code ) code
# define __DIAGASUFFIX( code ) \
         t_Path( code.string() + "." \
       + boost::lexical_cast<std::string>(boost::mpi::communicator().rank()) )
# define __IIAGA( code ) 
# else
# define __DIAGA( code ) 
# define __DIAGASUFFIX( code ) code
# define __IIAGA( code ) code
#endif

namespace LaDa
{

  //! \brief Holds everything pescan
  //! \details Escan is a semi-empirical potential code capable of computing
  //!          eigenvalues and related quantities for structure containing
  //!          thousands of atoms. The Hamiltonian is constructed from
  //!          pseudo-potentials which are fitted to LDA values and corrected for
  //!          the infamous LDA errors. It reproduces total energies to a feww
  //!          meVs and wavefunction with mean overlap with LDA wavefunctions of
  //!          99%. For more information, see <A
  //!          HREF="http://dx.doi.org/10.1103/PhysRevB.51.17398"> L-W Wang and
  //!          A.  Zunger PRB \b 51, 17398 (1995) </A>, and <A
  //!          HREF="http://dx.doi.org/10.1103/PhysRevB.66.045208"> K. Kim. P. R.
  //!          C. Kent, Alex Zunger and C. B. Geller, PRB \b 66, 045208 (2002)
  //!          </A>. 
  //!
  //!          The interface comes in one of two flavors: pescan (and genpot) can
  //!          be called as an external program specified on input, or it comes
  //!          as a library which is directly integrated at compile time. In the
  //!          former case, pescan itself can only be called in serial by
  //!          subsequent programs ( eg genetic algorithms ). In the second case,
  //!          it computes the eigenvalues using all the procs available to the
  //!          Pescan::Interface::comm communicator. This communicator must be
  //!          set prior to any computation through a call to
  //!          Pescan::Interface::set_mpi(). This last member is not available
  //!          when the interface is compiled for use with external programs. You
  //!          can use the __DIAGA and __IAGA macros to include/exclude code
  //!          depending on the type of call.
  namespace Pescan
  {

    //! Number of valence states in crystal structures.
    types::t_unsigned nb_valence_states( const Crystal::Structure &_str );
    //! Number of valence states in crystal structures.
    types::t_unsigned nb_valence_states( const Crystal::TStructure<std::string> &_str );

    //! \brief Defines an interface for the nanopse pescan program.
    //! \details Mostly, this class writes the input and recovers the output eigenvalues.
    class Interface __DIAGA( : public MPI_COMMDEC )
    {
      friend class boost::serialization::access;
      public:
        //! Path type.
        typedef boost::filesystem::path t_Path;
 
        //! Method for solving the eigenvalue problem
        enum t_method
        { 
          NOMET, //!< No method, placeholder.
          FOLDED_SPECTRUM, //!< Folded spectrum method.
          ALL_ELECTRON  //!< Full diagonalization method.
        };
 
        //! Contains parameters necessary for generating the potential.
        struct GenPot;
        //! Parameters for spin-orbit hamiltonian
        struct SpinOrbit;
        //! Parameters for nanopse's escan program.
        struct Escan;
        //! Whether to use verbose escan output.
        bool verbose;
 
        // Now includes declaration of Genpot, SpinOrbit, Escan classes.
#       include "interface.localclasses.h"
    
        //! Standard output filename
        t_Path stdout_file;
        //! Standard error filename
        t_Path stderr_file;
        //! Filename of the atomic configuation input
        t_Path atom_input;
        //! All potential generation parameters
        GenPot genpot;
        //! All escan-specific parameters
        Escan escan;
        //! Name of the maskr file.
        boost::filesystem::path maskr;
        //! Stores results
        std::vector<types::t_real> eigenvalues;
        //! Whether to delete directory where computations are being performed.
        bool do_destroy_dir;
 
        //! Constructor
        Interface()
          : stdout_file("out"), stderr_file("err"), atom_input("atom_config"), genpot(), escan(),
            maskr("maskr"), dirname("ESCAN"), do_destroy_dir(true), verbose(true) {}
        //! Copy Constructor
        Interface   ( const Interface &_c )
                  : __DIAGA( MPI_COMMCOPY( _c ) __COMMA__ )
                    stdout_file(_c.stdout_file), stderr_file(_c.stderr_file),
                    atom_input( _c.atom_input ), genpot( _c.genpot ),
                    escan( _c.escan ), maskr( _c.maskr ), eigenvalues( _c.eigenvalues ),
                    dirname( _c.dirname ), do_destroy_dir( _c.do_destroy_dir ), verbose(false) {}
        //! Destructor
       ~Interface() {};
 
       //! Loads all parameters from XML.
       bool Load( const TiXmlElement &_node );
 
       //! Sets the name of the directory where to perform calculations.
       void set_dirname( const boost::filesystem::path &_dir ) { dirname = _dir; }
       //! Sets the name of the directory where to perform calculations.
       const t_Path& get_dirname() const { return dirname; }
       //! Sets the reference energies for folded spectrum calculations.
       void set_reference( types::t_real _val )  { escan.Eref = _val; }
       //! Stores the reference energies for folded spectrum calculations.
       types::t_real get_reference() const { return escan.Eref; }
       //! return eigenvalue method to used.
       t_method get_method() const { return escan.method; }
       //! Set eigenvalue method to used.
       void set_method( t_method _method = FOLDED_SPECTRUM )
         { escan.method = _method; }
       //! Sets the name of the atomic configuration input file from Vff.
       void set_atom_input( const t_Path &_str ) { atom_input = _str; }
       //! Destroys directory for computations.
       void destroy_directory();
       //! \brief Sets the scale to that of the structure
       void set_scale( const Crystal::Structure &_str )
         { escan.scale = _str.scale; }
       //! \brief Sets the scale to that of the structure
       types::t_real get_scale() const { return escan.scale; }
       //! sets the number of states to compute
       void set_nbstates( const types::t_unsigned _n )
         { escan.nbstates = std::max<types::t_unsigned>(_n, 1); }
   
       void check_existence() const;
 
#      ifdef _DIRECTIAGA
         //! Gives access to ::mpi::AddCommunicator members. 
         void set_mpi( boost::mpi::communicator* _c );
         //! Gives access to ::mpi::AddCommunicator members. 
         const boost::mpi::communicator &comm() const { return MPI_COMMDEC::comm(); } 
         //! Gives access to ::mpi::AddCommunicator members. 
         boost::mpi::communicator &comm() { return MPI_COMMDEC::comm(); } 
#      endif
 
      //! Launches a calculation 
      bool operator()(); 

    protected:
      //! launches pescan once potential has been created.
      types::t_real launch_pescan(); 
      //! Creates directory for computations.
      void create_directory();
      //! Destroys directory for computations, depending on Interface::do_destroy_dir
      void destroy_directory_() { if (do_destroy_dir) destroy_directory(); }
      //! Interfaces to potential creation
      void create_potential();
      //! Writes escan parameter input.
      void write_escan_input();
      //! Writes potential generation parameter input.
      void write_genpot_input();
      //! Reads results from escan output.
      bool read_result();
      //! Loads functional directly from \a _node
      bool Load_( const TiXmlElement &_node );
 
      //! Directory where to perform computations.
      boost::filesystem::path dirname;
 
      private:
        //! Serializes escan functional.
        template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version)
        {
          _ar & verbose; _ar & stdout_file; _ar & stderr_file; _ar & atom_input; 
          _ar & genpot; _ar & escan; _ar & maskr; _ar & eigenvalues; 
          _ar & do_destroy_dir; _ar & dirname;
        }
    };
      
    //! Prints out the genpot input file.
    std::ostream& operator<<( std::ostream& _stream, const Interface::GenPot &_g );
 
    inline void Interface::check_existence() const 
    { 
      namespace bfs = boost::filesystem; 
      bool docontinue(true);
 
#     ifdef __thistry__
#       error macro already defined
#     endif
#     ifdef __endthistry__
#       error macro already defined
#     endif
#     define __thistry__ try {
#     define __endthistry__ \
        } \
        catch( std::exception &_e ) \
        { \
          docontinue = false; \
          Print::out << _e.what() << Print::endl; \
          std::cerr << _e.what() << std::endl; \
        } 
 
#     ifndef _NOLAUNCH
        __IIAGA
        (
          __thistry__
          __DOASSERT( not bfs::exists( genpot.launch ),
                      genpot.launch << " does not exist.\n" );
          __endthistry__
          __thistry__
          __DOASSERT( not bfs::exists( escan.launch ),
                      escan.launch << " does not exist.\n" );
          __endthistry__
        )
#     endif
 
      __thistry__
      __DOASSERT( not bfs::exists( maskr ), maskr << " does not exist.\n" );
      __DOASSERT( not ( bfs::is_regular_file( maskr ) or bfs::is_symlink( maskr ) ), 
                  maskr << " is not a regular file nor a symlink.\n" );
      __endthistry__
 
      t_Path file;
      std::vector<t_Path> :: const_iterator i_ps = genpot.pseudos.begin();
      std::vector<t_Path> :: const_iterator i_ps_end = genpot.pseudos.end();
      for(; i_ps != i_ps_end; ++i_ps )
      {
        if ( file == *i_ps ) continue;
        file = *i_ps;
        __thistry__
        __DOASSERT( not bfs::exists( file ), file << " does not exist.\n" );
        __DOASSERT( not ( bfs::is_regular_file( file ) or bfs::is_symlink( file ) ), 
                    file << " is not a regular file nor a symlink.\n" );
        __endthistry__
      }
 
      if( escan.potential != Escan::SPINORBIT ) return;
      std::vector<SpinOrbit> :: const_iterator i_so = escan.spinorbit.begin();
      std::vector<SpinOrbit> :: const_iterator i_so_end = escan.spinorbit.end();
      for(; i_ps != i_ps_end; ++i_ps )
      {
        if ( file == *i_ps ) continue;
        file = *i_ps;
        __thistry__
        __DOASSERT( not bfs::exists( file ), file << " does not exist.\n" );
        __DOASSERT( not ( bfs::is_regular_file( file ) or bfs::is_symlink( file ) ), 
                    file << " is not a regular file nor a symlink.\n" );
        __endthistry__
      }
      __DOASSERT( docontinue != true, "" )
#     undef __thistry__
#     undef __endthistry__
    }

  } // namespace pescan_interface
} // namespace LaDa

#endif // _PESCAN_INTERFACE_H_
