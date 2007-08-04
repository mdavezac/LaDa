#ifndef _PESCAN_INTERFACE_H_
#define _PESCAN_INTERFACE_H_

#include <string>
#include <vector>
#include <utility>

#include<lamarck/structure.h>

#include <opt/types.h>

#ifdef _MPI
#include <mpi/mpi_object.h>
#endif

namespace Pescan
{
  class Interface 
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<Interface> ( Interface& );
#endif 
  
    struct Bands
    { 
      types::t_real cbm, vbm; 
      Bands() : cbm(0), vbm(0) {};
      Bands( const types :: t_real _cbm, types::t_real _vbm ) : cbm(_cbm), vbm(_vbm) {};
      Bands( const Bands &_bands ) : cbm(_bands.cbm), vbm(_bands.vbm) {};
      types::t_real gap() const { return cbm - vbm; }
    };

    public:
    struct GenPot
    {
      std::string filename;
      types::t_int x, y, z; // real space mesh
      types::t_real cutoff; // Energy cutoff
      std::string output;
      std::string launch;
      std::vector<std::string> pseudos;

      GenPot   () 
             : filename("pot.input"), x(0), y(0), z(0), 
               cutoff(0), output("pot.output"), launch("getVLarg") {}
      bool check()
      {
        return x && y && z && cutoff && ( pseudos.size() >= 2 );
      }
    };
    struct SpinOrbit
    {
      std::string filename;
      std::string izz;
      types::t_real s, p, d, pnl, dnl;

      SpinOrbit () : filename(""), izz(""), s(0), p(0), d(0), pnl(0), dnl(0) {};
      SpinOrbit ( const SpinOrbit &_c) : filename(_c.filename), izz(_c.izz),
                                         s(_c.s), p(_c.p), d(_c.d), pnl(_c.pnl), dnl(_c.dnl) {};
    };
    struct Escan
    {
      public:
        enum t_method { NOMET, FOLDED_SPECTRUM, ALL_ELECTRON };
        enum t_potential { NOPOT, LOCAL, NONLOCAL, SPINORBIT };

      std::string filename;
      std::string output;
      std::string wavefunction;
      t_method method;
      Bands Eref;
      types::t_real smooth, kinscal;
      types::t_int nbstates;
      types::t_int niter, nlines;
      types::t_real tolerance;
      atat::rVector3d kpoint;
      types::t_real scale;
      t_potential potential;
      types::t_real rcut;
      std::string launch;
      std::vector<SpinOrbit> spinorbit;

      Escan () : filename("escan.input"), output("escan.out"), wavefunction("wg.cbm"), 
                 method(FOLDED_SPECTRUM), Eref(0,0), smooth(0.5), kinscal(0.0), nbstates(3),
                 niter(10), nlines(50), tolerance(types::tolerance),
                 kpoint(0,0,0), scale(0), potential(LOCAL), rcut(0), 
                 launch("escanCNL") {}
    };
    enum t_computation { VBM, CBM };
    protected:
      std::string atom_input;
      GenPot genpot;
      Escan escan;
      t_computation computation;
      std::string dirname;
      Bands bands;
      bool do_destroy_dir;

    public:
      Interface () : atom_input("atom.config"), genpot(), escan(), computation(VBM), do_destroy_dir(true) {}
     ~Interface() {};

     bool Load( const TiXmlElement &_node );

     types::t_real operator()( Ising_CE::Structure &_str ); 
     types::t_real launch_pescan( Ising_CE::Structure &_str ); 
     void set_dirname( const std::string &_dir ) { do_destroy_dir = false; dirname = _dir; }
     void set_references( types::t_real _val, types::t_real _cond )
       { escan.Eref.vbm = _val; escan.Eref.cbm = _cond; }
     void get_references( types::t_real &_val, types::t_real &_cond ) const
       { _val = escan.Eref.vbm; _cond = escan.Eref.cbm; }
     void get_bands( types::t_real &_VBM, types::t_real &_CBM ) const
       { _CBM = bands.cbm; _VBM = bands.vbm; }
     Escan::t_method get_method() const
       { return escan.method; }
     void set_method( Escan::t_method _method = Escan::FOLDED_SPECTRUM )
       { escan.method = _method; }

    protected:
     void create_directory();
     void destroy_directory();
     void create_potential();
     void write_escan_input( Ising_CE::Structure &_str );
     void write_genpot_input();
     types::t_real read_result( Ising_CE::Structure &_str );
  };

} // namespace pescan_interface

#endif // _PESCAN_INTERFACE_H_
