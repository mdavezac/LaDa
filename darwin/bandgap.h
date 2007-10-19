//
//  Version: $Id$
//
#ifndef _BANDGAP_H_
#define _BANDGAP_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include <vff/functional.h>
#include <pescan_interface/interface.h>
#include <lamarck/structure.h>
#include <opt/opt_function_base.h>
#include <opt/opt_minimize_gsl.h>
#include <opt/types.h>

#include "two_sites.h"
#include "evaluator.h"
#include "concentration.h"
#include "functors.h"
#include "individual.h"
#include "vff.h"
#include "pescan.h"


#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

namespace BandGap
{

  struct Object : public TwoSites::Object, public Pescan::Keeper
  {
    friend std::ostream& operator<<(std::ostream &_stream, const Object &_o);
    typedef TwoSites::Object :: t_Container t_Container;
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<BandGap::Object>(BandGap::Object &);
#endif
    types::t_real x, y;


    Object() : TwoSites::Object(), Pescan::Keeper() {}
    Object   (const Object &_c)
           : TwoSites::Object(_c), Pescan::Keeper(_c),
             x(_c.x), y(_c.y) {};
    ~Object() {};
  };

  inline std::ostream& operator<<(std::ostream &_stream, const Object &_o)
  { 
    if( _o.Container().size() <= 30 )
      _stream << (const SingleSite::Object& ) _o << " ";
    _stream << "x=" << _o.x << " y="  << _o.y 
            << (const Pescan::Keeper&) _o;
    return _stream; 
  } 

  typedef Individual::Types< BandGap::Object, 
                             TwoSites::Concentration, 
                             TwoSites::Fourier        > :: Scalar t_Individual;

  class Evaluator : public TwoSites::Evaluator< t_Individual >
  {
    public:
      typedef BandGap::t_Individual t_Individual;
      typedef Traits::GA< Evaluator > t_GATraits;
    protected:
      typedef TwoSites::Evaluator< t_Individual > t_Base;
      typedef Ising_CE::Structure::t_kAtoms t_kvecs;
      typedef Ising_CE::Structure::t_Atoms t_rvecs;

    public:
      using t_Base :: Load;
    protected:
      using t_Base :: current_individual;
      using t_Base :: current_object;

      
    protected:
      Vff::Darwin vff;
      Pescan::Darwin pescan;

    public:
      Evaluator() : t_Base(), vff(structure), pescan(structure)  {}
      ~Evaluator() {};

      bool Load( const TiXmlElement &_node );
      bool Load ( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      bool Save ( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const;
      void LoadAttribute ( const TiXmlAttribute &_att ) {};
      eoF<bool>* LoadContinue(const TiXmlElement &_el )
        { return new GA::mem_zerop_t<Pescan::Darwin>( pescan, &Pescan::Darwin::Continue,
                                                      "Pescan::Continue" );     }

      void evaluate();
  };


} // namespace BandGap

#ifdef _MPI
namespace mpi
{
  template<>
  inline bool mpi::BroadCast::serialize<BandGap::Object>( BandGap::Object & _object )
  {
    return     serialize<Pescan::Keeper>( _object )
           and serialize( _object.x )
           and serialize( _object.y )
           and _object.broadcast( *this );
  }
}
#endif

#endif // _BANDGAP_H_
