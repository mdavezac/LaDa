//
//  Version: $Id$
//
#ifndef _LAYERED_H_
#define _LAYERED_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <ostream>

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include <lamarck/structure.h>
#include <opt/types.h>

#include "evaluator.h"
#include "individual.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

//! \brief allows the creation
namespace Layered
{
  template <types::t_unsigned _D>
  class Concentration
  {
    protected:
      types::t_real x0[_D];
      types::t_unsigned N;
      types::t_int Nfreeze[_D];
      bool single_c;

    public:
      types::t_real x[_D];

    public:
      Concentration ();
      Concentration ( const Concentration &_conc);
      ~Concentration() {}


      void operator()( Ising_CE::Structure &_str );
      void operator()( Object &_obj );
      void operator()( const Ising_CE::Structure &_str, Object &_object,
                       types::t_int _concx, types::t_int _concy );
      void set( const Ising_CE::Structure &_str);
      void set( const Object &_obj );
      void setfrozen ( const Ising_CE::Structure &_str );

      std::string print() const;

    protected:
      void normalize( Ising_CE::Structure &_str, 
                      types::t_real _tochange);

  }

  class Physics 
  {
    protected:
      Ising_CE :: Lattice lattice;
      Ising_CE :: Structure structure;
      iVector3d direction;
      types::t_unsigned multiplicity;

    public:
      Evaluator() {}
      Evaluator   ( const Evaluator &_c )
                : lattice(_c.lattice), structure( _c.structure),
                  direction(_c.direction), multipicity(_c.multiplicity) {}
      ~Evaluator() {}

      bool Load( const TiXmlElement &_node );

    protected:
      bool Load_Structure( const TiXmlElement &_node );
  };

  class Depth
  {
    protected:
      rVector3d depth;

    public:
      Depth( const rVector3d &_vec ) : depth(_vec) {}
      Depth( const Depth &_c) : depth(_c.vec) {}

      bool operator()(const rVector3d& _1, const rVector3d& _2 );
  };

} // namespace Layered


#include "layered.impl.h"

#endif // _LAYERED_H_
