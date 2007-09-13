//
//  Version: $Id$
//
#ifndef _PESCAN_H_
#define _PESCAN_H_

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

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

namespace BandGap
{

  struct Object : public TwoSites::Object
  {
    friend std::ostream& operator<<(std::ostream &_stream, const Object &_o);
    typedef TwoSites::Object :: t_Container t_Container;
    typedef types::t_real t_Quantity;
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<BandGap::Object>(BandGap::Object &);
#endif
    types::t_real cbm, vbm;

    Object() : TwoSites::Object(), cbm(0), vbm(0) {}
    Object(const Object &_c) : TwoSites::Object(_c), cbm(_c.cbm), vbm(_c.vbm) {};
    ~Object() {};
  };


  class Evaluator : public TwoSites::Evaluator< Individual::Types<BandGap::Object>::Scalar >
  {
    public:
      typedef  Individual::Types<BandGap::Object>::Scalar t_Individual;
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
      Vff::Functional vff;
      Pescan::Interface pescan;
      minimizer::GnuSL<Vff::Functional> vff_minimizer;
      std::string references_filename;
      types::t_int nbeval;
      types::t_int age, check_ref_every;

    public:
      Evaluator() : t_Base(), 
                    vff(structure), vff_minimizer(vff), references_filename("BandEdge"), 
                    nbeval(0), age(0), check_ref_every(-1) {}
      ~Evaluator() {};

      bool initialize( t_Individual &_indiv )
        { return t_Base::initialize( _indiv ); }
      void init( t_Individual &_indiv );
      bool Load( const TiXmlElement &_node );
      bool Load ( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      bool Save ( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const;
      void LoadAttribute ( const TiXmlAttribute &_att ) {};
      eoF<bool>* LoadContinue(const TiXmlElement &_el )
        { return new darwin::mem_zerop_t<Evaluator>( *this, &Evaluator::Continue,
                                                     "Evaluator::Continue" );     }

      bool Continue();
      void evaluate();

    protected:
      void set_all_electron() { pescan.set_method( Pescan::Interface::Escan::ALL_ELECTRON ); }
      void read_references();
      void write_references();
      void get_bands( types::t_real &_vbm, types::t_real &_cbm ) 
        { pescan.get_bands( _vbm, _cbm); }
  };

  void rearrange_structure(Ising_CE::Structure &);

  inline std::ostream& operator<<(std::ostream &_stream, const Object &_o)
  { 
    _stream << (SingleSite::Object& ) _o 
            << " CBM " << std::fixed << std::setw(12) << std::setprecision(6) << _o.cbm 
            << "  --  VBM " << std::fixed << std::setw(12) << std::setprecision(6) << _o.vbm; 
    return _stream; 
  } 

} // namespace BandGap

#ifdef _MPI
namespace mpi
{
  template<>
  inline bool mpi::BroadCast::serialize<BandGap::Object>( BandGap::Object & _object )
  {
    if( not serialize( _object.cbm ) ) return false;
    if( not serialize( _object.vbm ) ) return false;
    return serialize< TwoSites::Object >( _object );
  }
}
#endif

#endif // _PESCAN_H_
