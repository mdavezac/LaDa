//
//  Version: $Id$
//
#ifndef _VFF_H_
#define _VFF_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>

#include <tinyxml/tinyxml.h>

#include <vff/functional.h>
#include <lamarck/structure.h>
#include <opt/opt_function_base.h>
#include <opt/opt_minimize_gsl.h>
#include <opt/types.h>
#include <atat/vectmac.h>

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

namespace Vff
{

  struct Keeper 
  {
    types::t_real energy;
    atat::rMatrix3d stress;

    Keeper() : energy(0), stress() {}
    Keeper(const Keeper &_c) : energy(_c.energy), stress(_c.stress) {};
    ~Keeper() {};

    bool Load( const TiXmlElement &_node );
    bool Save( TiXmlElement &_node ) const;
  };


  class Base : public Functional
  {
    public:
      Base( Ising_CE::Structure &_str ) : Functional( _str ) {} 
      ~Base() {};

      Load( const TiXmlElement &_node );
      bool evaluate( Keeper &_keeper )
      {
        minimizer.minimize();
        structure.energy = vff.energy();
        _keeper.energy = structure.energy;
        _keeper.stress = stress;
      }
  };

} // namespace BandGap

#ifdef _MPI
namespace mpi
{
  template<>
  inline bool mpi::BroadCast::serialize< Vff::Keeper >( Vff::Keeper & _k )
  {
    return     serialize( _k.energy ) 
           and serialize( _k.stress ); 
  }
}
#endif

#endif // _BANDGAP_H_
