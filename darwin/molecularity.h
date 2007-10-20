//
//  Version: $Id$
//
#ifndef _MOLECULARITY_H_
#define _MOLECULARITY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <ostream>

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include <vff/functional.h>
#include <pescan_interface/interface.h>
#include <lamarck/structure.h>
#include <opt/opt_function_base.h>
#include <opt/opt_minimize_gsl.h>
#include <opt/types.h>

#include "two_sites.h"
#include "pescan.h"
#include "vff.h"
#include "evaluator.h"
#include "individual.h"
#include "gaoperators.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

namespace Molecularity
{
  struct Fourier
  {
    template<class T_R_IT, class T_K_IT>
    Fourier( T_R_IT _rfirst, T_R_IT _rend,
             T_K_IT _kfirst, T_K_IT _kend );
    template<class T_R_IT, class T_K_IT, class T_O_IT >
    Fourier( T_R_IT _rfirst, T_R_IT _rend,
             T_K_IT _kfirst, T_K_IT _kend,
             T_O_IT _rout ); // sets rvector values from kspace values
  };

  // Exact same object as SingleSite::Object, except for the two friends
  class Object : public SingleSite::Object, public Vff::Keeper, public Pescan::Keeper
  {
    protected:
      typedef SingleSite :: Object t_Base;

    public:
      typedef t_Base :: t_Type t_Type;
      typedef t_Base :: t_Container t_Container;

    public:
      Object() {}
      Object(const Object &_c) : t_Base(_c), Vff::Keeper(_c), Pescan::Keeper(_c) {};
      ~Object() {};
      
      types::t_int get_concentration() const
      {
        t_Container :: const_iterator i_var = bitstring.begin();
        t_Container :: const_iterator i_var_end = bitstring.end();
        types::t_int result = 0;
        for(; i_var != i_var_end; ++i_var )
          result += *i_var > 0 ? 1: -1;
        return result;
      }

      bool Load( const TiXmlElement &_node )
        { return Vff::Keeper::Load(_node) and Pescan::Keeper::Load(_node); }
      bool Save( TiXmlElement &_node ) const
        { return Vff::Keeper::Save(_node) and Pescan::Keeper::Save(_node); }
  };

  class Concentration : public SingleSite::Concentration
  {
    typedef SingleSite::Concentration t_Base;
    public:
      Concentration () : t_Base() {}
      Concentration   ( const Concentration &_c) 
                    : t_Base(_c) {}

      void operator()( Ising_CE::Structure &_str);
      void operator()( Object &_obj) { t_Base::operator()( _obj ); }
      void set( const Ising_CE::Structure &_structure);
      void set( const Object &_obj) { t_Base::set( _obj ); }
      void setfrozen( const Ising_CE::Structure &_structure);

    protected:
      void normalize( Ising_CE::Structure &_str, types::t_real _tochange );
  };

  inline std::ostream& operator<<(std::ostream &_stream, const Object &_o)
    { return _stream << (const SingleSite::Object&) _o; }
  inline void operator<<(Object &_o, const std::string &_str )
    { (SingleSite::Object&) _o << _str; }
  inline void operator<<(std::string &_str, const Object &_o)
    { _str << (const SingleSite::Object&) _o; }

  void operator<<(Ising_CE::Structure &_str, const Object &_o);
  void operator<<(Object &_o, const Ising_CE::Structure &_str);


  typedef Individual::Types< Object, 
                             Concentration, 
                             Fourier        > :: Vector t_Individual;

  class Evaluator : public TwoSites::Evaluator< Molecularity::t_Individual >
  {
    public:
      typedef Molecularity::t_Individual t_Individual;
      typedef Traits::GA< Evaluator > t_GATraits;
    protected:
      typedef Evaluator t_This;
      typedef TwoSites::Evaluator< t_Individual > t_Base;
      typedef Ising_CE::Structure::t_kAtoms t_kvecs;
      typedef Ising_CE::Structure::t_Atoms t_rvecs;

    public:
      using t_Base :: Load;
      using t_Base :: Save;
    protected:
      using t_Base :: current_individual;
      using t_Base :: current_object;

    protected:
      Pescan::Darwin pescan;
      Vff::Darwin vff;

    public:
      Evaluator() : t_Base(), pescan(structure), vff(structure) {}
      Evaluator   ( const Evaluator &_c )
                : lattice( _c.lattice ), structure( _c.structure ),
                  concentration( _c.concentration ), pescan(_c.pescan), vff(_c.vff) {}
      ~Evaluator() {}

      bool Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const;
      bool Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      bool Load( const TiXmlElement &_node );
      void LoadAttribute ( const TiXmlAttribute &_att )
        { concentration.LoadAttribute( _att ); };

      void evaluate()
      {
        // relax structure
        vff();
        // Load relaxed structure into pescan
        pescan << vff; 
        // get band gap
        pescan( *current_object );
      
        // set quantity
        current_individual->quantities().front() =
           ( current_object->stress(0,0) + current_object->stress(1,1) ) * 0.5;
        current_individual->quantities().back() = current_object->cbm - current_object->vbm;
      }
    protected:
      bool consistency_check();
  };

  template<class T_R_IT, class T_K_IT>
  Fourier :: Fourier( T_R_IT _rfirst, T_R_IT _rend,
                      T_K_IT _kfirst, T_K_IT _kend )
  {
    const std::complex<types::t_real>
       imath(0, -2*3.1415926535897932384626433832795028841971693993751058208);
    
    for (; _kfirst != _kend; ++_kfirst)
    {
      _kfirst->type = std::complex<types::t_real>(0);
      for(T_R_IT i_r( _rfirst ); i_r != _rend; i_r += 2 )
      {
        _kfirst->type +=   exp( imath * ( i_r->pos[0] * _kfirst->pos[0] +
                                          i_r->pos[1] * _kfirst->pos[1] +
                                          i_r->pos[2] * _kfirst->pos[2] ) )
                         * i_r->type;
      }
    }
  }
  template<class T_R_IT, class T_K_IT, class T_O_IT >
  Fourier :: Fourier( T_R_IT _rfirst, T_R_IT _rend,
                      T_K_IT _kfirst, T_K_IT _kend,
                      T_O_IT _rout ) // sets rvector values from kspace values
  {
    const std::complex<types::t_real>
       imath(0, 2*3.1415926535897932384626433832795028841971693993751058208);
    for (; _rfirst != _rend; _rfirst+=2, ++_rout)
    {
      *_rout = 0.0;
      for(T_K_IT i_k=_kfirst; i_k != _kend; ++i_k)
      {
        *_rout +=   exp( imath * ( _rfirst->pos[0] * i_k->pos[0] +
                                   _rfirst->pos[1] * i_k->pos[1] +
                                   _rfirst->pos[2] * i_k->pos[2] ) )
                  * i_k->type;
      }
    }
  }
} // namespace BandGap


#ifdef _MPI
namespace mpi
{
  template<>
  inline bool mpi::BroadCast::serialize<Molecularity::Object>( Molecularity::Object & _object )
  {
    return     serialize<Pescan::Keeper>( _object )
           and serialize<Vff::Keeper>( _object )
           and _object.broadcast( *this );
  }
}
#endif

#endif // _MOLECULARITY_H_
