//
//  Version: $Id$
//
#ifndef _PESCAN_H_
#define _PESCAN_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <tinyxml/tinyxml.h>

#include <pescan_interface/interface.h>
#include <lamarck/structure.h>
#include <opt/types.h>

#include "vff.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

namespace Pescan
{

  struct Keeper
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<Pescan::Keeper>(Pescan::Keeper &);
#endif
    types::t_real cbm, vbm;

    Keeper() : cbm(0), vbm(0) {}
    Keeper(const Keeper &_c) : cbm(_c.cbm), vbm(_c.vbm) {};
    ~Keeper() {};

    bool Load( const TiXmlElement &_node );
    bool Save( TiXmlElement &_node ) const;
  };


  class Darwin 
  {
    public:
      std::string atomicconfig;
      
    protected:
      Ising_CE::Structure &structure;
      Pescan::Interface pescan;
      std::string references_filename;
      std::string dirname;
      types::t_int nbeval;
      types::t_int age, check_ref_every;

    public:
      Darwin   ( Ising_CE::Structure &_s )
             : structure(_s), references_filename("BandEdge"), 
               nbeval(0), age(0), check_ref_every(-1) {}
      Darwin   ( const Darwin &_b ) 
             : structure(_b.structure), references_filename(_b.references_filename),
               nbeval(_b.nbeval), age(_b.age), check_ref_every(_b.check_ref_every) {}
      ~Darwin() {};

      bool Load( const TiXmlElement &_node );
      bool Continue();
      void operator()();
      void operator()( Keeper &_keeper );
      void operator<<( const Vff::Darwin &_vff );

    protected:
      void set_all_electron() { pescan.set_method( Pescan::Interface::Escan::ALL_ELECTRON ); }
      void read_references();
      void write_references();
      void get_bands( types::t_real &_vbm, types::t_real &_cbm ) { pescan.get_bands( _vbm, _cbm); }
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
    get_bands( _keeper.vbm, _keeper.cbm );
  }

} // namespace Pescan

#ifdef _MPI
namespace mpi
{
  template<>
  inline bool mpi::BroadCast::serialize<Pescan::Keeper>( Pescan::Keeper & _keeper )
  {
    return     serialize( _keeper.cbm ) 
           and serialize( _keeper.vbm );
  }
}
#endif

#endif // _PESCAN_H_
