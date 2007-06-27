#ifndef _CE_H_
#define _CE_H_

#include <string>
#include <algorithm>
#include <functional>
#include <string>

#include <vff/functional.h>
#include <pescan_interface/interface.h>
#include <lamarck/structure.h>
#include <opt/opt_function_base.h>
#include <opt/opt_minimize_gsl.h>
#include <opt/types.h>

#include <tinyxml/tinyxml.h>

#include <eo/eoOp.h>

#include "evaluator.h"
#include "concentration.h"
#include "functors.h"

#ifdef _MPI
#include<mpi/mpi_object.h>
#endif

namespace BandGap
{
  // actual vff cum pescan functional
  class Functional : public function :: Base<>
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<BandGap::Functional>( BandGap::Functional & );
#endif
    public:
      typedef types::t_real t_Type;
      typedef std::vector<t_Type> t_Container;
      typedef t_Container :: iterator iterator;
      typedef t_Container :: const_iterator const_iterator;

    protected:
      Ising_CE::Structure &structure;
      Ising_CE::Structure evaluation_structure;
      Vff::Functional vff;
      Pescan::Interface pescan;
      minimizer::GnuSL<Vff::Functional> vff_minimizer;
      std::string band_edge_filename;

    public:
      Functional   ( Ising_CE::Structure &_str, std::string _f = "BandEdge" ) 
                 : structure( _str ), vff(evaluation_structure),
                   vff_minimizer(vff), band_edge_filename(_f) {}
      Functional   ( const Functional &_func )
                 : structure(_func.structure), evaluation_structure(_func.evaluation_structure),
                   vff(_func.vff), pescan(_func.pescan), vff_minimizer(_func.vff_minimizer),
                   band_edge_filename(_func.band_edge_filename) {};
      ~Functional() {};
      bool Load( const TiXmlElement &_node );
      
      void set_filename( const std::string &_f )
        { band_edge_filename = _f; } 
      void read_band_edges();
      void write_band_edges();
      void set_all_electron() { pescan.set_method( Pescan::Interface::Escan::ALL_ELECTRON ); }
      t_Type evaluate()
      {
        // first minimizes strain
        vff_minimizer.minimize();
        structure.energy = vff.energy();
        vff.print_escan_input();

        // then evaluates band gap
        types::t_real result = pescan(structure);
        if ( pescan.get_method() != Pescan::Interface::Escan::ALL_ELECTRON )
          return result;
        
#ifdef _MPI
        if ( mpi::main.rank() != mpi::ROOT_NODE ) // not root no read write
          return result;
#endif 

        write_band_edges();
        pescan.set_method(); // resets to folded spectrum if necessary
        return result;
      } 
      void evaluate_gradient( t_Type* const _i_grad ) 
      { 
        std::cerr << "Not implemented !!" << std::endl;
        exit(0);
      }
      t_Type evaluate_one_gradient( types::t_unsigned _pos)
      { 
        std::cerr << "Not implemented !!" << std::endl;
        exit(0);
      }
      t_Type evaluate_with_gradient( t_Type* const _i_grad )
      { 
        std::cerr << "Not implemented !!" << std::endl;
        exit(0);
      }
  };

  struct Object 
  {
    typedef types::t_real t_Type;
    typedef std::vector<t_Type>  t_Container;
    typedef t_Container :: iterator iterator;
    typedef t_Container :: const_iterator const_iterator;
    friend void operator<<(std::string &_str, const Object &_o);
    friend void operator<<(Ising_CE::Structure &_str, const Object &_o);
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<BandGap::Object>(BandGap::Object &);
#endif
    t_Container bitstring;
    Object() {}
    Object(const Object &_c) : bitstring(_c.bitstring) {};
    bool operator<<(const Ising_CE::Structure &_c)
    {
      bitstring.clear(); bitstring.reserve( _c.atoms.size() );
      Ising_CE::Structure :: t_Atoms :: const_iterator i_atom = _c.atoms.begin();
      Ising_CE::Structure :: t_Atoms :: const_iterator i_end = _c.atoms.end();
      for(; i_atom != i_end; ++i_atom )
        if ( not (i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T) )
          bitstring.push_back( i_atom->type > 0 ? 1.0: -1.0 );
      return true;
    }
    bool operator<<(const std::string &_c)
    {
      types::t_unsigned size = _c.size();
      bitstring.resize( size );
      std::vector<types::t_real> :: iterator i_var = bitstring.begin();
      std::vector<types::t_real> :: iterator i_end = bitstring.end();
      for(types::t_unsigned n=0; i_var != i_end; ++i_var, ++n )
        *i_var = ( _c[n] == '1' ) ? 1.0: -1.0;
      return true;
    }
    ~Object() {};
    
    bool operator==( const Object &_c ) const
    {
      return std::equal( bitstring.begin(), bitstring.end(), 
                          _c.bitstring.begin() ); 
    }
    t_Container :: const_iterator begin() const
      { return bitstring.begin();  }
    t_Container:: const_iterator end() const
      { return bitstring.end();  }
    t_Container:: iterator begin() 
      { return bitstring.begin();  }
    t_Container:: iterator end()
      { return bitstring.end();  }

    types::t_real get_concentration() const
    {
      t_Container :: const_iterator i_var = bitstring.begin();
      t_Container :: const_iterator i_var_end = bitstring.end();
      types::t_real result = 0.0;
      for(; i_var != i_var_end; ++i_var )
        result += *i_var > 0 ? 1.0: -1.0;
      result /= static_cast<types::t_real>(bitstring.size());
      return result;
    }

    void print_out( std::ostream &_stream ) const
      { std::string str; str << *this; _stream << str; }

    void mask( types::t_unsigned _start, types::t_unsigned _end)
    {
      if ( _end > bitstring.size() )
        _end = bitstring.size();
      std::transform( bitstring.begin()+_start, bitstring.begin()+_end,
                      bitstring.begin()+_start, std::logical_not<t_Type>() );  
    }
    types::t_unsigned size()
      { return bitstring.size(); }
  };


  class Evaluator : public darwin::Evaluator<Object>
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<BandGap::Evaluator>( BandGap::Evaluator & );
#endif
    public:
      typedef function::Base<> t_Functional;
    public:
      using darwin::Evaluator<Object> :: Load;
      // typedef t_Object is inherited from darwin::Evaluator<...> 
    protected:
      typedef Ising_CE::Structure::t_kAtoms t_kvecs;
      typedef Ising_CE::Structure::t_Atoms t_rvecs;
      
    protected:
      Ising_CE::Lattice lattice;
      Ising_CE::Structure structure;
      Functional functional;
      types::t_real crossover_probability;
      types::t_real x, y;
      types::t_real lessthan, morethan;
      X_vs_y x_vs_y;
      bool single_concentration;
      types::t_int age, check_ref_every;

    public:
      Evaluator() : functional(structure), crossover_probability(0.5),
                    x(0), y(0), lessthan(1.0), morethan(-1.0),
                    single_concentration(false), age(0), check_ref_every(-1) {}; 
      ~Evaluator() {};

      void* const init( t_Object &_object );
      bool Load( const TiXmlElement &_node );
      bool Load ( t_Object &_indiv, const TiXmlElement &_node, bool _type );
      bool Save ( const t_Object &_indiv, TiXmlElement &_node, bool _type ) const;
      void LoadAttribute ( const TiXmlAttribute &_att ) {};
      eoOp<t_Object>* LoadGaOp(const TiXmlElement &_el );
      eoF<bool>* LoadContinue(const TiXmlElement &_el )
        { return new darwin::mem_zerop_t<Evaluator, t_Object>( *this, &Evaluator::Continue,
                                                               "Evaluator::Continue" );     }
      eoMonOp<const t_Object>* LoadTaboo(const TiXmlElement &_el );

      bool Krossover( t_Object  &_offspring, const t_Object &_parent,
                      bool _range = false );
      bool Crossover( t_Object &_obj1, const t_Object &_obj2 );
      bool Continue();

      bool initialize( t_Object &_object );
      void set_object( t_Object &_object, const void * const _f ) {}
      void* const LoadMinimizer(const TiXmlElement &_el );

    protected:
      void get_xy_concentrations( const Ising_CE::Structure &_str );
      bool consistency_check();
      void set_concentration( Ising_CE::Structure &_str );
      bool Taboo(const t_Object &_object );
  };

  void rearrange_structure(Ising_CE::Structure &);
  types::t_real set_concentration( Ising_CE::Structure &_str,
                                   types::t_int _site,
                                   types::t_real _target );

// class SiteIterator
// {
//   protected: 
//     Ising_CE::Structure &structure;
//     Ising_CE::Structure::t_Atoms::iterator current_iterator;
//     enum t_where { NOWHERE, BEGIN, END };
//     enum t_which { SITES0, SITES1, ALLSITES } which;
//     
//   public:
//     SiteIterator   ( const Isince_CE::Structure &_str, 
//                      t_where _where = NOWHERE, t_which _i = ALLSITES ) 
//                  : structure( _str )
//     {
//       which = _i;
//       if ( _where == BEGIN )
//         current_iterator = _str.atoms.begin();
//       if ( _where == END )
//         current_iterator = _str.atoms.end();
//     }
//     SiteIterator   ( const SiteIterator &_site )
//                  : structure( _site.structure ), current_iterator( _site.current_iterator ),
//                    which( _site.which ) {}
//     void which_sites ( t_which _w ) { which = _w; }
//     void operator++()
//     {
//       if ( which == ALLSITES )
//         { ++ current_iterator; return; }
//       
//       current_iterator += 2;
//     }
//     Ising_CE::Structure::t_Atom & operator*() 
//     {
//       if ( which == ALLSITES  or which == SITES0 )
//         return *current_iterator;
//       return *(current_iterator+1);
//     }
//     Ising_CE::Structure::t_Atom * operator->() 
//     {
//       if ( which == ALLSITES  or which == SITES0 )
//         return &(*current_iterator);
//       return &(*(current_iterator+1));
//     }
//     bool operator == ( const SiteIterator &_it )
//     {
//       return current_iterator == _it.current_iterator;
//     }
//     bool operator != ( const SiteIterator &_it )
//     {
//       return current_iterator != _it.current_iterator;
//     }
//     
// };

  template<class T_R_IT, class T_K_IT>
  void fourrier_to_kspace( T_R_IT _rfirst, T_R_IT _rend,
                           T_K_IT _kfirst, T_K_IT _kend ) // sets kvector values from rspace values
  {
    const std::complex<types::t_real> imath(0, -2*3.1415926535897932384626433832795028841971693993751058208);
    
    for (; _kfirst != _kend; ++_kfirst)
    {
      _kfirst->type = std::complex<types::t_real>(0);
      for(T_R_IT i_r( _rfirst ); i_r != _rend; ++i_r )
      {
        _kfirst->type +=   exp( imath * ( i_r->pos[0] * _kfirst->pos[0] +
                                          i_r->pos[1] * _kfirst->pos[1] +
                                          i_r->pos[2] * _kfirst->pos[2] ) )
                         * std::complex<types::t_real>(i_r->type, (++i_r)->type);
      }
    }
  }
  template<class T_R_IT, class T_K_IT, class T_O_IT >
  void fourrier_to_rspace( T_R_IT _rfirst, T_R_IT _rend,
                           T_K_IT _kfirst, T_K_IT _kend,
                           T_O_IT _rout ) // sets rvector values from kspace values
  {
    const std::complex<types::t_real> imath(0, 2*3.1415926535897932384626433832795028841971693993751058208);
    for (; _rfirst != _rend; ++_rfirst, ++_rout)
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
    bool consistency_check();
  }


} // namespace BandGap

types::t_real set_concentration( Ising_CE::Structure &_str,
                                 types::t_int _site,
                                 types::t_real _target);

#endif // _CE_OBJECT_H_
