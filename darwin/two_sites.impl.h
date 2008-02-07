//
//  Version: $Id$
//
#ifndef _TWOSITES_IMPL_H_
#define _TWOSITES_IMPL_H_

#include <algorithm>

#include <print/xmg.h>
#include <print/stdout.h>

namespace TwoSites
{
  template<class T_R_IT, class T_K_IT>
  Fourier :: Fourier( T_R_IT _rfirst, T_R_IT _rend,
                      T_K_IT _kfirst, T_K_IT _kend )
  {
    const std::complex<types::t_real>
       imath(0, -2*3.1415926535897932384626433832795028841971693993751058208);
    
    for (; _kfirst != _kend; ++_kfirst)
    {
      _kfirst->type = std::complex<types::t_real>(0);
      for(T_R_IT i_r( _rfirst ); i_r != _rend; ++i_r )
      {
        // used to be directly in std::complex constructor below...
        // but then it hit me: which constructor argument does c++ look at first?
        // this is safe now.
        types::t_real a = i_r->type;
        types::t_real b = (++i_r)->type;
        _kfirst->type +=   exp( imath * ( i_r->pos[0] * _kfirst->pos[0] +
                                          i_r->pos[1] * _kfirst->pos[1] +
                                          i_r->pos[2] * _kfirst->pos[2] ) )
                         * std::complex<types::t_real>(a, b);
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


  inline bool Concentration :: Load( const TiXmlElement &_node )
  {
    if( not X_vs_y::Load( _node ) ) return false;
    if( not single_c )  return true;
    x = get_x();  y = get_y();

    return true;
  }


  template<class T_INDIVIDUAL>
  bool Evaluator<T_INDIVIDUAL> :: Load( t_Individual &_indiv,
                                        const TiXmlElement &_node,
                                        bool _type )
  {
    if ( _type == GA::LOADSAVE_SHORT )
    {
      if( not _node.Attribute("string") )
        return false;
      (_indiv.Object()) << std::string(_node.Attribute("string"));
      return true;
    }

    Ising_CE::Structure s; 
    if ( not s.Load(_node) )
      return false;
    (_indiv.Object()) << s;
    return true;
  }

  template<class T_INDIVIDUAL>
  bool Evaluator<T_INDIVIDUAL> :: Save( const t_Individual &_indiv, 
                                        TiXmlElement &_node, 
                                        bool _type ) const
  {
    if ( _type == GA::LOADSAVE_SHORT )
    {
      std::string str; str << _indiv.Object();
      _node.SetAttribute("string", str.c_str());
      return true;
    }

    Ising_CE::Structure s = structure; 
    s << _indiv.Object();
    TwoSites::Fourier( s.atoms.begin(),  s.atoms.end(),
                       s.k_vecs.begin(), s.k_vecs.end() );
    s.print_xml(_node);

    return true;
  }

  template<class T_INDIVIDUAL>
  inline bool Evaluator<T_INDIVIDUAL> :: Load( const TiXmlElement &_node )
  {
    if ( not lattice.Load( _node ) )
    {
      std::cerr << " Could not load lattice type from input!! " << std::endl; 
      return false;
    }
    Ising_CE::Structure::lattice = &lattice;
    if ( not structure.Load( _node ) )
    {
      std::cerr << " Could not load input structure!! " << std::endl; 
      return false;
    }
    if ( not structure.set_site_indices() )
    {
      std::cerr << " Could not set atomic indices! " << std::endl; 
      return false;
    }
    rearrange_structure(structure);
    if ( not consistency_check() )  return false;

    if ( not concentration.Load( _node ) ) 
    {
      std::cerr << " Could not load Concentration input!! " << std::endl; 
      return false;
    }

    concentration.setfrozen( structure );
    Print::xmg << Print::Xmg::comment << concentration.print() << Print::endl;
    
    return true;
  }

  template<class T_INDIVIDUAL>
  bool Evaluator<T_INDIVIDUAL> :: consistency_check()
  {
    Ising_CE::Structure::t_Atoms :: iterator i_atom = structure.atoms.begin();
    Ising_CE::Structure::t_Atoms :: iterator i_atom_end = structure.atoms.end();
    for(; i_atom != i_atom_end; ++i_atom )
      if ( i_atom->site != 1 and i_atom->site != 0 )
        return false;
    i_atom = structure.atoms.begin();

    for(; i_atom != i_atom_end; ++i_atom )
      if ( i_atom->site == 0 and not (i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
        break;
    if (     i_atom == i_atom_end
         and not (lattice.sites[0].freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
      lattice.sites[0].freeze |=  Ising_CE::Structure::t_Atom::FREEZE_T;
    if (     i_atom != i_atom_end 
         and (lattice.sites[0].freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
      for(i_atom = structure.atoms.begin(); i_atom != i_atom_end; i_atom+=2 )
        i_atom->freeze |= Ising_CE::Structure::t_Atom::FREEZE_T;

    for(i_atom = structure.atoms.begin(); i_atom != i_atom_end; ++i_atom )
      if ( i_atom->site == 1 and not (i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
        break;
    if (     i_atom == i_atom_end
         and not (lattice.sites[1].freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
      lattice.sites[1].freeze |=  Ising_CE::Structure::t_Atom::FREEZE_T;
    if (     i_atom != i_atom_end
         and (lattice.sites[1].freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
      for(i_atom = structure.atoms.begin(); i_atom != i_atom_end; i_atom+=2 )
        (i_atom+1)->freeze |= Ising_CE::Structure::t_Atom::FREEZE_T;

    if (     ( lattice.sites[0].freeze & Ising_CE::Structure::t_Atom::FREEZE_T )
         and ( lattice.sites[1].freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
    {
      std::cerr << "No atoms to optimize !? " << std::endl;
      return false;
    }
    return true;
  }


  template<class T_INDIVIDUAL>
  void Evaluator<T_INDIVIDUAL> :: presubmit( std::list<t_Individual> &_pop )
  {
    t_Individual pure;
    initialize( pure );
    // Pure A
    std::fill( pure.Object().Container().begin(), pure.Object().Container().end(), 1.0 );
    _pop.push_back(pure);
    // Pure B
    std::fill( pure.Object().Container().begin(), pure.Object().Container().end(), -1.0 );
    _pop.push_back(pure);
  } 

  template<class T_INDIVIDUAL> inline GA::Taboo_Base<T_INDIVIDUAL>* 
    Evaluator<T_INDIVIDUAL> :: LoadTaboo(const TiXmlElement &_el )
    {
      if ( concentration.single_c ) return NULL;
      GA::xTaboo<t_Individual> *xtaboo =
          new GA::xTaboo< t_Individual >( concentration );
      if ( xtaboo and xtaboo->Load( _el ) )  return xtaboo;
      if ( xtaboo ) delete xtaboo;
      return NULL;
    }

  template<class T_INDIVIDUAL> 
    inline bool Evaluator<T_INDIVIDUAL> :: initialize( t_Individual &_indiv )
    {
      GA::Random< t_Individual > random( concentration, structure, _indiv );
      _indiv.invalidate(); return true;
    }
  template<class T_INDIVIDUAL> 
    inline void Evaluator<T_INDIVIDUAL> :: init( t_Individual &_indiv )
    {
      t_Base :: init( _indiv );
      // sets structure to this object 
      structure << *current_object;
    }  
} // namespace TwoSites
#endif // _TWOSITES_IMPL_H_
