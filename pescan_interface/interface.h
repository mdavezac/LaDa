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

#include<lamarck/structure.h>

#include <opt/types.h>

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

  //! Keeps track of the HOMO and LUMO
  struct Bands
  { 
    types::t_real cbm;  //!< Conduction Band Minimum
    types::t_real vbm;  //!< Valence Band minimum
    Bands() : cbm(0), vbm(0) {}; //!< Constructor
    //! Constructor and Initializer
    Bands( const types :: t_real _cbm, types::t_real _vbm ) : cbm(_cbm), vbm(_vbm) {};
    //! Copy Constructor.
    Bands( const Bands &_bands ) : cbm(_bands.cbm), vbm(_bands.vbm) {};
    //! Returns the gap.
    types::t_real gap() const { return cbm - vbm; }
#ifdef _NOLAUNCH
    //! Debug: for _NOLAUNCH fake functional only
    void swap() { types::t_real a = cbm; cbm = vbm; vbm = a; }
#endif
  };

  //! \brief Defines an interface for the nanopse pescan program.
  //! \details Mostly, this class writes the input and recovers the output.
  //! \todo Isolate the pescan launcher itself from the band-gap computation?
  class Interface 
  {
    public:

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
      //! Some coherency checks.
      bool check() { return x && y && z && cutoff && ( pseudos.size() >= 2 ); }
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
    };
    //! Parameters for nanopse's escan program.
    struct Escan
    {
      public:
        //! Method for solving the eigenvalue problem
        enum t_method
        { 
          NOMET, //!< No method, placeholder.
          FOLDED_SPECTRUM, //!< Folded spectrum method.
          ALL_ELECTRON  //!< Full diagonalization method.
        };
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
      //! Eigenvalue solution method.
      t_method method;
      //! Reference energy for folded spectrum method
      Bands Eref;
      //! \brief Index of VBM and CBM wavefunctions
      std::pair<types::t_unsigned> wfn_index;
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
      //! Real-space scale 
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
                 method(FOLDED_SPECTRUM), Eref(0,0), smooth(0.5), kinscal(0.0), nbstates(3),
                 niter(10), nlines(50), tolerance(types::tolerance),
                 kpoint(0,0,0), scale(0), potential(LOCAL), rcut(0), 
                 launch("escanCNL") {}
    };
    //! For folded spectrum, whcih band is being computed.
    enum t_computation
    { 
      VBM,  //!< Valence Band Maximum is being computed.
      CBM   //!< Conduction Band Minimum is being computed.
    };

#ifdef _MPI
    //! \cond
    friend bool mpi::BroadCast::serialize<Interface> ( Interface& );
    //! \endcond
#endif 
#ifdef _NOLAUNCH
    //! \cond
    friend void nolaunch_functional( const Ising_CE::Structure &_str,
                                     Interface::Bands &bands );
    //! \endcond
#endif

  
    protected:
      //! Filename of the atomic configuation input
      std::string atom_input;
      //! All potential generation parameters
      GenPot genpot;
      //! All escan-specific parameters
      Escan escan;
      //! For folded spectrum, which computation is being performed.
      t_computation computation;
      //! Directory where to perform computations.
      std::string dirname;
      //! Stores results
      Bands bands;
      //! Whether to delete directory where computations are being performed.
      bool do_destroy_dir;
      //! Whether to input wavefunctions
      bool do_input_wavefunctions;

    public:
      //! Constructor
      Interface () : atom_input("atom.config"), genpot(), escan(),
                     computation(VBM), do_destroy_dir(true),
                     do_input_wavefunctions(false) {}
      //! Destructor
     ~Interface() {};

     //! Loads all parameters from XML.
     bool Load( const TiXmlElement &_node );

     //! Launches a band-gap calculation for structure \a _str
     types::t_real operator()( Ising_CE::Structure &_str ); 
     
     //! launches pescan once potential has been created.
     types::t_real launch_pescan( Ising_CE::Structure &_str ); 
     //! Sets the name of the directory where to perform calculations.
     void set_dirname( const std::string &_dir ) { dirname = _dir; }
     //! Sets the reference energies for folded spectrum calculations.
     void set_references( types::t_real _val, types::t_real _cond )
       { escan.Eref.vbm = _val; escan.Eref.cbm = _cond; }
     //! Stores the reference energies for folded spectrum calculations.
     void get_references( types::t_real &_val, types::t_real &_cond ) const
       { _val = escan.Eref.vbm; _cond = escan.Eref.cbm; }
     //! returns the band edges.
     void get_bands( types::t_real &_VBM, types::t_real &_CBM ) const
       { _CBM = bands.cbm; _VBM = bands.vbm; }
     //! Returns the current eigenvalue method being used.
     Escan::t_method get_method() const { return escan.method; }
     //! Set eigenvalue method to used.
     void set_method( Escan::t_method _method = Escan::FOLDED_SPECTRUM )
       { escan.method = _method; }
     //! Sets the name of the atomic configuration input file from Vff.
     void set_atom_input( const std::string &_str ) { atom_input = _str; }
     //! \brief sets and codes name of output wavefunctions according to calculation
     //! \details If \a _name is not empty, Escan::wavefunction_out is set to \a _name.
     //!          The function returns a name which codes for the kind of
     //!          computation done
     const std::string &wfn_name(const std::string _name);
     //! Destroys directory for computations.
     void destroy_directory();

    protected:
     //! Creates directory for computations.
     void create_directory();
     //! Destroys directory for computations, depending on Interface::do_destroy_dir
     void destroy_directory_() { if (do_destroy_dir) destroy_diretory(); }
     //! Interfaces to potential creation
     void create_potential();
     //! Writes escan parameter input.
     void write_escan_input( Ising_CE::Structure &_str );
     //! Writes potential generation parameter input.
     void write_genpot_input();
     //! Reads results from escan output.
     types::t_real read_result( Ising_CE::Structure &_str );
  };


  inline const std::string Interface :: wfn_name(std::string _name)
  {
    if( escan.method == FOLDED_SPECTRUM )
      _name +=  "." + (computation == VBM ? "vbm": "cbm" )
    else
      _name +=  ".ae";
    return _name;
  }

} // namespace pescan_interface

#endif // _PESCAN_INTERFACE_H_
