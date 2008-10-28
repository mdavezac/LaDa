//
//  Version: $Id$
//

// Convenience file to move local class declaration outside Pescan::Interface
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
    //! Possible real-space wavefunction output.
    enum t_rWfnOutput
    {
      NOOUTPUT, //!< No output.
      DENSITY, //!< Square of wavefunctions.
      DENSITY_AFTER_CALL, //!< Square of wavefunctions after calling escan.
      WFN, //!< Wavefunctions after calling escan.
      WFN_AFTER_CALL  //!< Wavefunctions after calling escan.
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
    __IIAGA
    (
      //! \brief System call to lauch nanopse's pescan.
      //! \details With __IIAGA defined only.
      boost::filesystem::path launch;
    )
    //! Spin orbit parameters.
    std::vector<SpinOrbit> spinorbit;
    //! real-space wavefunction  output.
    t_rWfnOutput rspace_output;
    //! Path to real-space wavefunctions.
    t_Path rspace_wfn;

    //! Constructor.
    Escan () : filename("escan.input"), output("escan.out"),
               wavefunction_out("wavefunction"), wavefunction_in("wavefunction"), 
               method(FOLDED_SPECTRUM), Eref(0), smooth(0.5), kinscal(0.0),
               nbstates(3), niter(10), nlines(50),
               tolerance(types::tolerance), kpoint(0,0,0), scale(0),
               potential(LOCAL), rcut(0), __IIAGA( launch("escanCNL") __COMMA__ )
               rspace_output(NOOUTPUT), rspace_wfn( "rzerowfn" )
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
            __IIAGA( launch(_c.launch) __COMMA__ )
            rspace_output( _c.rspace_output ), rspace_wfn( _c.rspace_wfn ) {}
  };

  struct GenPot
  {
    //! Filenam where to write the output.
    t_Path filename;
    types::t_int x; //!< Number of points in x coordinate of real space mesh
    types::t_int y; //!< Number of points in y coordinate of real space mesh
    types::t_int z; //!< Number of points in z coordinate of real space mesh
    types::t_real cutoff; //!< Plane-wave energy cutoff 
    t_Path output;   //!< File to which to write the potential
    __IIAGA
    ( 
      t_Path launch;   //!< Command for launching pescan's getVLarg. 
    )
    //! Name of the files describing each pseudo-potentials
    std::vector<t_Path> pseudos;

    //! Constructor
    GenPot   () 
           : filename("pot.input"), x(0), y(0), z(0), 
             cutoff(0), output("pot.output")
             __IIAGA( __COMMA__ launch("getVLarg") ) {}
    //! Copy Constructor
    GenPot   ( const GenPot & _c )
           : filename( _c.filename ), x( _c.x ), y( _c.y ), z( _c.z ),
             cutoff( _c.cutoff ), output( _c.output ),
             __IIAGA( launch( _c.launch ) __COMMA__ )
             pseudos( _c.pseudos ) {}
    //! Some coherency checks.
    bool check() { return x && y && z && cutoff && ( pseudos.size() >= 2 ); }
  };
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


