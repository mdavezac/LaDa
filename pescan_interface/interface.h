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
#include <boost/lexical_cast.hpp>

#include <string>
#include <vector>
#include <utility>
#include <fstream>

#include <tinyxml/tinyxml.h>

#include <physics/physics.h>
#include <opt/types.h>
#include <atat/vectmac.h>
#include <crystal/structure.h>
#include <mpi/mpi_object.h>
#include <print/stdout.h>

#ifdef _DIRECTIAGA
# define __DIAGA( code ) code
# define __DIAGASUFFIX( code ) \
         t_Path( code.string() + "." + boost::lexical_cast<std::string>(MPI_COMM.rank()) )
# define __IIAGA( code ) 
# else
# define __DIAGA( code ) 
# define __DIAGASUFFIX( code ) code
# define __IIAGA( code ) code
#endif


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

  //! \brief Defines an interface for the nanopse pescan program.
  //! \details Mostly, this class writes the input and recovers the output eigenvalues.
  class Interface __DIAGA( : public MPI_COMMDEC )
  {
    protected:
      //! Path type.
      typedef boost::filesystem::path t_Path;
    public:
      //! Method for solving the eigenvalue problem
      enum t_method
      { 
        NOMET, //!< No method, placeholder.
        FOLDED_SPECTRUM, //!< Folded spectrum method.
        ALL_ELECTRON  //!< Full diagonalization method.
      };

    protected:

    //! Contains parameters necessary for generating the potential.
    struct GenPot
    {
      //! Filenam where to write the output.
      t_Path filename;
      types::t_int x; //!< Number of points in x coordinate of real space mesh
      types::t_int y; //!< Number of points in y coordinate of real space mesh
      types::t_int z; //!< Number of points in z coordinate of real space mesh
      types::t_real cutoff; //!< Plane-wave energy cutoff 
      t_Path output;   //!< File to which to write the potential
      t_Path launch;   //!< Command for launching pescan's getVLarg.
      //! Name of the files describing each pseudo-potentials
      std::vector<t_Path> pseudos;

      //! Constructor
      GenPot   () 
             : filename("pot.input"), x(0), y(0), z(0), 
               cutoff(0), output("pot.output"), launch("getVLarg") {}
      //! Copy Constructor
      GenPot   ( const GenPot & _c )
             : filename( _c.filename ), x( _c.x ), y( _c.y ), z( _c.z ),
               cutoff( _c.cutoff ), output( _c.output ), launch( _c.launch ),
               pseudos( _c.pseudos ) {}
      //! Some coherency checks.
      bool check() { return x && y && z && cutoff && ( pseudos.size() >= 2 ); }
    };
    //! Parameters for spin-orbit hamiltonian
    struct SpinOrbit
    {
      t_Path filename; //!< Filename of the spin-orbit empirical pseudo-potential
      std::string izz; //!< Don't know.
      types::t_real s;   //!< Don't know. 
      types::t_real p;   //!< Don't know. 
      types::t_real d;   //!< Don't know.   
      types::t_real pnl; //!< Don't know.
      types::t_real dnl; //!< Don't know.

      //! Constructor.
      SpinOrbit () : filename(""), izz(""), s(0), p(0), d(0), pnl(0), dnl(0) {};
      //! Copy Constructor.
      SpinOrbit   ( const SpinOrbit &_c)
                : filename(_c.filename), izz(_c.izz), s(_c.s), 
                  p(_c.p), d(_c.d), pnl(_c.pnl), dnl(_c.dnl) {};
    };
    //! Parameters for nanopse's escan program.
    struct Escan
    {
        //! Type of the potential
        enum t_potential
        { 
          NOPOT, //!< No potential, place-holder.
          LOCAL,  //!< Local potential only.
          NONLOCAL, //!< Non-local potential
          SPINORBIT //!< Include spin-orbit.
        };

      //! Name where to write input.
      t_Path filename;
      //! Name of the pescan output
      t_Path output;
      //! Stub name for output wavefunction files
      t_Path wavefunction_out;
      //! Stub name for input wavefunction files
      t_Path wavefunction_in;
      //! Wavefunctions to read in
      std::vector<types::t_unsigned> read_in;
      //! Eigenvalue solution method.
      t_method method;
      //! Reference energy for folded spectrum method
      types::t_real Eref;
      //! Kinetic scaling factor for plane-wave cutoff
      types::t_real kinscal;
      //! Smoothness fcator plane-wave cutoff
      types::t_real smooth;
      //! The number of states to compute.
      types::t_int nbstates;
      //! Should be the number of iterations. Right?
      types::t_int niter;
      //! God and Peter only know.
      types::t_int nlines;
      //! Tolerance for diagonalization convergence.
      types::t_real tolerance;
      //! Reciprocal-space vector for which to compute eigenvalue.
      atat::rVector3d kpoint;
      //! \brief Reciprocal-Space scale \f$\frac{2\pi}{a}\f$, with \e a the lattice
      //! constant Crystal::Structure::scale
      types::t_real scale;
      //! Type of Hamiltonian to solve.
      t_potential potential;
      //! Real-space cutoff?
      types::t_real rcut;
      //! System call to lauch nanopse's pescan
      boost::filesystem::path launch;
      //! Spin orbit parameters.
      std::vector<SpinOrbit> spinorbit;

      //! Constructor.
      Escan () : filename("escan.input"), output("escan.out"),
                 wavefunction_out("wavefunction"), wavefunction_in("wavefunction"), 
                 method(FOLDED_SPECTRUM), Eref(0), smooth(0.5), kinscal(0.0),
                 nbstates(3), niter(10), nlines(50),
                 tolerance(types::tolerance), kpoint(0,0,0), scale(0),
                 potential(LOCAL), rcut(0), launch("escanCNL")
        { read_in.clear(); read_in.reserve(nbstates); }
      //! Copy Constructor
      Escan   ( const Escan &_c)
            : filename(_c.filename), output(_c.output),
              wavefunction_out( _c.wavefunction_out),
              wavefunction_in( _c.wavefunction_in ), read_in( _c.read_in ),
              method(_c.method), Eref(_c.Eref), smooth( _c.smooth ),
              kinscal( _c.kinscal), nbstates(_c.nbstates), niter(_c.niter),
              nlines(_c.nlines), tolerance(_c.tolerance), kpoint(_c.kpoint),
              scale(_c.scale), potential(_c.potential), rcut(_c.rcut),
              launch(_c.launch) {}
    };
  
    protected:
      //! Filename of the atomic configuation input
      t_Path atom_input;
      //! All potential generation parameters
      GenPot genpot;
      //! All escan-specific parameters
      Escan escan;
      //! Directory where to perform computations.
      boost::filesystem::path dirname;
      //! Name of the maskr file.
      boost::filesystem::path maskr;
      //! Whether to delete directory where computations are being performed.
      bool do_destroy_dir;
     
    public:
      //! Stores results
      std::vector<types::t_real> eigenvalues;

    public:
      //! Constructor
      Interface()
        : __DIAGA( MPI_COMMCOPY( *::mpi::main ) MPI_COMMA )
          atom_input("atom.config"), genpot(), escan(),
          maskr("maskr"), dirname("ESCAN"), do_destroy_dir(true) {}
      //! Copy Constructor
      Interface   ( const Interface &_c )
                : __DIAGA( MPI_COMMCOPY( _c ) MPI_COMMA )
                  atom_input( _c.atom_input ), genpot( _c.genpot ),
                  escan( _c.escan ), maskr( _c.maskr ),
                  dirname( _c.dirname ), do_destroy_dir( _c.do_destroy_dir ),
                  eigenvalues( _c.eigenvalues ) {}
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
     //! \brief Tries to find a <Functional type="escan"> tag in \a _node.
     //! \details Checks wether \a _node or its immediate offpsrings are the
     //!          right functional node.
     const TiXmlElement* find_node (const TiXmlElement &_node );
     //! \brief Sets the scale to that of the structure
     void set_scale( const Crystal::Structure &_str )
       { escan.scale = _str.scale; }
     //! sets the number of states to compute
     void set_nbstates( const types::t_unsigned _n )
       { escan.nbstates = std::max<types::t_unsigned>(_n, 1); }
 
     void check_existence() const;

     __DIAGA
     (
       //! Allows derived classes to have access to ::mpi::AddCommunicator members. 
       void set_mpi( boost::mpi::communicator* _c );
     )

#    ifdef _DIRECTIAGA
       //! Allows derived classes to have access to ::mpi::AddCommunicator members. 
       const boost::mpi::communicator &comm() const { return MPI_COMMDEC::comm(); } 
#    endif

   protected:
    //! Launches a calculation 
    bool operator()(); 
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
#   ifdef _DIRECTIAGA
      //! Allows derived classes to have access to ::mpi::AddCommunicator members. 
      boost::mpi::communicator &comm() { return MPI_COMMDEC::comm(); } 
#   endif
  };


  inline bool Interface::Load( const TiXmlElement &_node )
  {
    const TiXmlElement *parent = find_node( _node );
    if ( parent ) return Load( *parent );
    
    std::cerr << "Could not find an <Functional type=\"escan\"> tag in input file" 
              << std::endl;
    return false;
  }
    

  inline void Interface::check_existence() const 
  { 
    namespace bfs = boost::filesystem; 
    bool docontinue(true);

#   ifdef __thistry__
#     error macro already defined
#   endif
#   ifdef __endthistry__
#     error macro already defined
#   endif
#   define __thistry__ try {
#   define __endthistry__ \
      } \
      catch( std::exception &_e ) \
      { \
        docontinue = false; \
        Print::out << _e.what() << Print::endl; \
        std::cerr << _e.what() << std::endl; \
      } 

#   ifndef _NOLAUNCH
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
#   endif

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
    __DOASSERT( docontinue == true, "" )
#   undef __thistry__
#   undef __endthistry__
  }

} // namespace pescan_interface

#endif // _PESCAN_INTERFACE_H_
