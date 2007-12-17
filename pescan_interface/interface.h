//
//  Version: $Id$
//
#ifndef _PESCAN_INTERFACE_H_
#define _PESCAN_INTERFACE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <string>
#include <vector>
#include <utility>

#include <physics/physics.h>
#include <opt/types.h>
#include <atat/vectmac.h>
#include <lamarck/structure.h>

#include <tinyxml/tinyxml.h>

#ifdef _MPI
#include <mpi/mpi_object.h>
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
namespace Pescan
{

  //! \brief Defines an interface for the nanopse pescan program.
  //! \details Mostly, this class writes the input and recovers the output eigenvalues.
  class Interface 
  {
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
      std::string filename;
      types::t_int x; //!< Number of points in x coordinate of real space mesh
      types::t_int y; //!< Number of points in y coordinate of real space mesh
      types::t_int z; //!< Number of points in z coordinate of real space mesh
      types::t_real cutoff; //!< Plane-wave energy cutoff 
      std::string output;   //!< File to which to write the potential
      std::string launch;   //!< Command for launching pescan's getVLarg.
      std::vector<std::string> pseudos; //!< Name of the files describing each pseudo-potentials

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
#ifdef _MPI
      //! \ingroup MPI
      //! \brief mpi broadcasting for this object
      //! \details There seems to be pbs with gcc 3.4 
      //!          when using friend template mpi::BroadCast::Serialize().
      bool broadcast( mpi::BroadCast& );
#endif 
    };
    //! Parameters for spin-orbit hamiltonian
    struct SpinOrbit
    {
      std::string filename; //!< Filename of the spin-orbit empirical pseudo-potential
      std::string izz; //!< Don't know.
      types::t_real s;   //!< Don't know. 
      types::t_real p;   //!< Don't know. 
      types::t_real d;   //!< Don't know.   
      types::t_real pnl; //!< Don't know.
      types::t_real dnl; //!< Don't know.

      //! Constructor.
      SpinOrbit () : filename(""), izz(""), s(0), p(0), d(0), pnl(0), dnl(0) {};
      //! Copy Constructor.
      SpinOrbit ( const SpinOrbit &_c) : filename(_c.filename), izz(_c.izz),
                                         s(_c.s), p(_c.p), d(_c.d), pnl(_c.pnl), dnl(_c.dnl) {};
#ifdef _MPI
      //! \ingroup MPI
      //! \brief mpi broadcasting for this object
      //! \details There seems to be pbs with gcc 3.4 
      //!          when using friend template mpi::BroadCast::Serialize().
      bool broadcast( mpi::BroadCast& _bc );
#endif 
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
      std::string filename;
      //! Name of the pescan output
      std::string output;
      //! Stub name for output wavefunction files
      std::string wavefunction_out;
      //! Stub name for input wavefunction files
      std::string wavefunction_in;
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
      //! constant Ising_CE::Structure::scale
      types::t_real scale;
      //! Type of Hamiltonian to solve.
      t_potential potential;
      //! Real-space cutoff?
      types::t_real rcut;
      //! System call to lauch nanopse's pescan
      std::string launch;
      //! Spin orbit parameters.
      std::vector<SpinOrbit> spinorbit;

      //! Constructor.
      Escan () : filename("escan.input"), output("escan.out"),
                 wavefunction_out("wavefunction"), wavefunction_in("wavefunction"), 
                 method(FOLDED_SPECTRUM), Eref(0), smooth(0.5), kinscal(0.0), nbstates(3),
                 niter(10), nlines(50), tolerance(types::tolerance),
                 kpoint(0,0,0), scale(0), potential(LOCAL), rcut(0), 
                 launch("escanCNL") { read_in.clear(); read_in.reserve(nbstates); }
      //! Copy Constructor
      Escan   ( const Escan &_c)
            : filename(_c.filename), output(_c.output), wavefunction_out( _c.wavefunction_out), 
              wavefunction_in( _c.wavefunction_in ), read_in( _c.read_in ), method(_c.method), 
              Eref(_c.Eref), smooth( _c.smooth ), kinscal( _c.kinscal), nbstates(_c.nbstates), 
              niter(_c.niter), nlines(_c.nlines), tolerance(_c.tolerance), kpoint(_c.kpoint), 
              scale(_c.scale), potential(_c.potential), rcut(_c.rcut), launch(_c.launch) {}
#ifdef _MPI
      //! \ingroup MPI
      //! \brief mpi broadcasting for this object
      //! \details There seems to be pbs with gcc 3.4 
      //!          when using friend template mpi::BroadCast::Serialize().
      bool broadcast( mpi::BroadCast& _bc );
#endif 
    };
#ifdef _MPI
    //! \cond
    friend bool mpi::BroadCast::serialize<Interface> ( Interface& );
    //! \endcond
#endif 
  
    protected:
      //! Filename of the atomic configuation input
      std::string atom_input;
      //! All potential generation parameters
      GenPot genpot;
      //! All escan-specific parameters
      Escan escan;
      //! Directory where to perform computations.
      std::string dirname;
      //! Whether to delete directory where computations are being performed.
      bool do_destroy_dir;
     
    public:
      //! Stores results
      std::vector<types::t_real> eigenvalues;

    public:
      //! Constructor
      Interface () : atom_input("atom.config"), genpot(), escan(),
                     dirname("escan"), do_destroy_dir(true) {}
      //! Copy Constructor
      Interface   ( const Interface &_c )
                : atom_input( _c.atom_input ), genpot( _c.genpot ),
                  escan( _c.escan ), 
                  dirname( _c.dirname ), do_destroy_dir( _c.do_destroy_dir ),
                  eigenvalues( _c.eigenvalues ) {}
      //! Destructor
     ~Interface() {};

     //! Loads all parameters from XML.
     bool Load( const TiXmlElement &_node );

     //! Sets the name of the directory where to perform calculations.
     void set_dirname( const std::string &_dir ) { dirname = _dir; }
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
     void set_atom_input( const std::string &_str ) { atom_input = _str; }
     //! Destroys directory for computations.
     void destroy_directory();
     //! \brief Tries to find a <Functional type="escan"> tag in \a _node.
     //! \details Checks wether \a _node or its immediate offpsrings are the
     //!          right functional node.
     const TiXmlElement* find_node (const TiXmlElement &_node );

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
     //! \brief Sets scale of reciprocal in mesh as \f$\frac{2\pi}{a}\f$, with \e a the
     //! lattice in constant in atomic units.
     void set_scale( const Ising_CE::Structure &_str )
       { escan.scale = _str.scale; }
  };


  inline bool Interface::Load( const TiXmlElement &_node )
  {
    const TiXmlElement *parent = find_node( _node );
    if ( parent ) return Load( *parent );
    
    std::cerr << "Could not find an <Functional type=\"escan\"> tag in input file" 
              << std::endl;
    return false;
  }
    

} // namespace pescan_interface

#endif // _PESCAN_INTERFACE_H_
