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
#include "bitstring.h"
#include "gaoperators.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

//! \brief allows the creation
namespace Layered
{
  template<types::t_unsigned _D>
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


  //! \brief Layered Object Type
  //! \detail Redefines BitString::Object with the sole purpose of implementing
  //!         overloaded operator<<() on and from Ising_CE::Structures.
  //!         Note that we could run into awkward redefinitions by using typedefs
  //!         rather than a complete redeclaration.
  template<class T_CONTAINER = std::vector<types::t_real> >
  class Object: public BitString :: Object<T_CONTAINER>
  {
    public:
      typedef T_CONTAINER t_Container;
      typedef typename t_Container :: value_type t_Type;
      typedef typename t_Container :: iterator iterator;
      typedef typename t_Container :: const_iterator const_iterator;
    protected:
      typedef BitString :: Object<T_CONTAINER> t_Base;

    public:
      Object() {}
      Object(const Object &_c) : t_Base(_c) {};
      Object(const t_Container &_c) : t_Base(_c) {};
      ~Object() {};
  };

  template<class T_CONTAINER>
    void operator<<(Ising_CE::Structure &_str, const Object<T_CONTAINER> &_o);
  template<class T_CONTAINER>
    void operator<<(Object<T_CONTAINER> &_o, const Ising_CE::Structure &_str);

  template <types::t_unsigned _D>
  class Concentration
  {
    public:
      const types::t_unsigned _d;

    protected:
      types::t_real x0;
      types::t_unsigned N;
      types::t_int Nfreeze;
      bool single_c;

    public:
      types::t_real x;

    public:
      Concentration() : _d(_D), x0(0.0), N(0), Nfreeze(0),
                        single_c(false), x(0) {}
      Concentration   ( const Concentration &_c)
                    : _d(_D), x0(_c.x0), N(_c.N), Nfreeze(_c.Nfreeze),
                      single_c(_c.single_c), x(_c.x) {}
      ~Concentration() {}


      void operator()( Ising_CE::Structure &_str );
      template<class T_CONT> void operator()( BitString::Object<T_CONT> &_obj );
      void set( const Ising_CE::Structure &_str);
      template<class T_CONT> void set( const BitString::Object<T_CONT> &_obj );
      void setfrozen ( const Ising_CE::Structure &_str );

      std::string print() const;

      void LoadAttribute ( const TiXmlAttribute &_att );

      bool is_single_c() const { return single_c; }

    protected:
      void normalize( Ising_CE::Structure &_str, 
                      types::t_real _tochange);

  };

  template< class T_INDIVIDUAL >
  class Evaluator : public GA::Evaluator< T_INDIVIDUAL >
  {
    public:
      typedef T_INDIVIDUAL t_Individual;

    protected:
      typedef typename t_Individual ::t_IndivTraits t_IndivTraits;
      typedef typename t_IndivTraits::t_Concentration t_Concentration;
      typedef typename t_IndivTraits::t_FourierRtoK t_FourierRtoK;
      typedef typename t_IndivTraits::t_FourierKtoR t_FourierKtoR;
      typedef GA::Evaluator<t_Individual> t_Base;
      typedef Evaluator<t_Individual>     t_This;

    protected:
      Ising_CE :: Lattice lattice;
      Ising_CE :: Structure structure;
      atat::rVector3d direction;
      types::t_unsigned multiplicity;
      t_Concentration concentration;

    public:
      using t_Base::Load;

    public:
      Evaluator() : t_Base() {}
      Evaluator   ( const Evaluator &_c )
                : t_Base(), lattice(_c.lattice), structure( _c.structure),
                  direction(_c.direction), multiplicity(_c.multiplicity),
                  concentration(_c.concentration) {}
      ~Evaluator() {}

      bool Load( const TiXmlElement &_node );
      bool Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      bool Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const;

      eoGenOp<t_Individual>* LoadGaOp(const TiXmlElement &_el )
       { return GA::LoadGaOp<t_Individual>( _el, structure, concentration ); }
      GA::Taboo_Base<t_Individual>* LoadTaboo(const TiXmlElement &_el );
      bool initialize( t_Individual &_indiv );
      void LoadAttribute ( const TiXmlAttribute &_att )
        { concentration.LoadAttribute( _att ); };

    protected:
      bool Load_Structure( const TiXmlElement &_node );
  };

  class Depth
  {
    protected:
      atat::rVector3d depth;

    public:
      Depth( const atat::rVector3d &_vec ) : depth(_vec) {}
      Depth( const Depth &_c) : depth(_c.depth) {}

      bool operator()(const atat::rVector3d& _1, const atat::rVector3d& _2 );
  };

} // namespace Layered


#include "layered.impl.h"

#endif // _LAYERED_H_
