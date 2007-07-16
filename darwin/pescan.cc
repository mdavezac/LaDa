#include <functional>
#include <algorithm>
#include <ext/algorithm>
#include <fstream>
#ifndef __PGI
  #include<ext/functional>
  using __gnu_cxx::compose1;
#else
  #include<functional>
  using std::compose1;
#endif

#include "pescan.h"
#include "print_xmgrace.h"
#include "lamarck/atom.h"
#include "opt/va_minimizer.h"

namespace BandGap
{
  void operator<<(std::string &_str, const Object &_o)
  {
    Object::t_Container :: const_iterator i_var = _o.bitstring.begin();
    Object::t_Container :: const_iterator i_end = _o.bitstring.end();
    std::ostringstream sstr;
    for(; i_var != i_end; ++i_var )
      sstr << ( *i_var > 0 ? '1' : '0' );
    _str = sstr.str();
  }
  void operator<<(Ising_CE::Structure &_str, const Object &_o)
  {
    Object::t_Container :: const_iterator i_var = _o.bitstring.begin();
    Object::t_Container :: const_iterator i_end = _o.bitstring.end();
    Ising_CE::Structure :: t_Atoms :: iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure :: t_Atoms :: iterator i_atom_end = _str.atoms.end();
    for(; i_var != i_end and i_atom != i_atom_end; ++i_atom )
    {
      if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) 
        continue;
      i_atom->type = *i_var > 0 ? 1.0 : -1.0;
      ++i_var;
    }
  }

  bool Evaluator :: Load( Object &_obj, const TiXmlElement &_node, bool _type )
  {
    if ( _type == darwin::LOADSAVE_SHORT )
    {
      if( not _node.Attribute("string") )
        return false;
      _obj << std::string(_node.Attribute("string"));
      return true;
    }

    Ising_CE::Structure s; 
    if ( not s.Load(_node) )
      return false;
    _obj << s;
    return true;
  }
  bool Evaluator :: Save( const Object &_obj, TiXmlElement &_node, bool _type ) const
  {
    if ( _type == darwin::LOADSAVE_SHORT )
    {
      std::string str; str << _obj;
      _node.SetAttribute("string", str.c_str());
      return true;
    }

    Ising_CE::Structure s = structure; 
    s << _obj;
    BandGap::fourrier_to_kspace( s.atoms.begin(),  s.atoms.end(),
                                 s.k_vecs.begin(), s.k_vecs.end() );
    s.print_xml(_node);

    return true;
  }
  bool Evaluator :: Load( const TiXmlElement &_node )
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
    if ( not consistency_check() )
      return false;
    if ( not x_vs_y.Load( _node ) )
    {
      std::cerr << " Could not load Concentration input!! " << std::endl; 
      return false;
    }
    if ( x_vs_y.is_singlec() )
      { x = x_vs_y.get_x(); y = x_vs_y.get_y(); }
    if ( not functional.Load( _node ) )
    {
      std::cerr << " Could not load functional input!! " << std::endl; 
      return false;
    }
    if (     _node.FirstChildElement("Filenames") 
         and _node.FirstChildElement("Filenames")->Attribute("BandEdge") )
      functional.set_filename( _node.FirstChildElement("Filenames")->Attribute("BandEdge") );

    return true;
  }

  bool Evaluator :: Crossover ( t_Object &_obj1, const t_Object &_obj2 )
  {
    Object::t_Container :: iterator i_var1 = _obj1.begin();
    Object::t_Container :: const_iterator i_var2 = _obj2.begin();
    Object::t_Container :: const_iterator i_var2_end = _obj2.end();
    for(; i_var2 != i_var2_end; ++i_var1, ++i_var2)
      if ( rng.flip(crossover_probability) ) 
        *i_var1 = *i_var2;
    structure << _obj1;
    Evaluator::set_concentration( structure );
    _obj1 << structure;
    return true;
  }

  // expects kspace value to exist!!
  bool Evaluator :: Krossover( Object  &_offspring, const Object &_parent,
                               bool _range )
  {
    Ising_CE::Structure str1 = structure, str2 = structure;
    str1 << _offspring; str2 << _parent;
    BandGap::fourrier_to_kspace( str1.atoms.begin(),  str1.atoms.end(),
                                 str1.k_vecs.begin(), str1.k_vecs.end() );
    BandGap::fourrier_to_kspace( str2.atoms.begin(),  str2.atoms.end(),
                                 str2.k_vecs.begin(), str2.k_vecs.end() );
    if ( _range and str1.k_vecs.size() > 2 ) // range crossover ... kvec should be oredered according to size
    {  
      types::t_unsigned n = (types::t_unsigned)
         std::floor(   (types::t_real) rng.random ( str1.k_vecs.size() - 1 ) 
                     * (types::t_real) crossover_probability );
      __gnu_cxx::copy_n( str2.k_vecs.begin(), n, str1.k_vecs.begin() );
    }
    else // every point crossover
    {
      Ising_CE::Structure::t_kAtoms :: const_iterator i_p = str2.k_vecs.begin();
      Ising_CE::Structure::t_kAtoms :: const_iterator i_p_end = str2.k_vecs.end();
      Ising_CE::Structure::t_kAtoms :: iterator i_o = str1.k_vecs.begin();
      for ( ; i_p != i_p_end; ++i_p, ++i_o)
        if ( rng.flip(crossover_probability) ) 
          i_o->type = i_p->type;
    }
  
    Evaluator::set_concentration( str1 );
    _offspring << str1;

    return true; // offspring has changed!
  }

  eoOp<Object>* Evaluator::LoadGaOp(const TiXmlElement &_el )
  {
    std::string value = _el.Value();

    if ( value.compare("Crossover") == 0 )
    {
      _el.Attribute( "prob", &crossover_probability );
      crossover_probability = crossover_probability > 0 ? std::abs(crossover_probability) : 0.5;
      if ( crossover_probability > 1 ) crossover_probability = 0.5;
      std::ostringstream sstr;
      sstr << "Crossover rate = " << crossover_probability;
      darwin::printxmg.add_comment(sstr.str());
      // pointer is owned by caller !!
      return new darwin::mem_binop_t<Evaluator, Object, void>( *this, &Evaluator::Crossover, 
                                                               std::string( "Crossover" ) );
    }
    else if ( value.compare( "Krossover" ) == 0 )
    {
      bool att = false;
      _el.Attribute( "prob", &crossover_probability );
      crossover_probability = crossover_probability > 0 ? std::abs(crossover_probability) : 0.5;
      if ( crossover_probability > 1 ) crossover_probability = 0.5;
      std::ostringstream sstr;
      sstr << "Krossover, rate = " << crossover_probability;
      darwin::printxmg.add_comment(sstr.str());
      if ( _el.Attribute("type") )
      {
        std::string str =  _el.Attribute("type");
        if ( str.compare("range") == 0 ) 
          { att = true; darwin::printxmg.add_to_last( ", Range = true" ); }
      }
      // pointer is owned by caller !!
      return new darwin::mem_binop_t<Evaluator, Object, bool>( *this, &Evaluator::Krossover, 
                                                               std::string( "Krossover" ), att);
    }

    return NULL;
  }

  bool Evaluator :: Taboo(const t_Object &_object )
  {
    if ( x_vs_y.is_singlec() )
      return false;
    structure << _object;
    get_xy_concentrations( structure );
    return x > lessthan or x < morethan; // if true, _object is taboo
  }
  
  eoMonOp<const Object>* Evaluator :: LoadTaboo(const TiXmlElement &_el )
  {
    if ( x_vs_y.is_singlec() )
      return NULL;
    const TiXmlElement *child = _el.FirstChildElement( "Concentration" );
    if ( not child )
      return NULL;
    double d;
    if ( child->Attribute( "lessthan" ) )
      child->Attribute( "xlessthan", &d );
    lessthan = ( d > 0 and d < 1 ) ? 2.0*d-1.0: 1.0;
    if ( child->Attribute( "morethan" ) )
      child->Attribute( "morethan", &d );
    morethan = ( d > 0 and d < 1 ) ? 2.0*d-1.0: -1.0;
    if ( lessthan < morethan )
      return NULL;
   
    std::ostringstream sstr;
    sstr << std::fixed << std::setprecision(3) << "Taboo x in [ " << 0.5*(morethan+1.0)
         << ", "  << 0.5*(lessthan+1.0) << "] ";
    darwin::printxmg.add_comment(sstr.str());
    // pointer is owned by caller !!
    return new darwin::const_mem_monop_t<Evaluator, Object>( *this, &Evaluator::Taboo,
                                                             "Taboo" );
  }

  bool Evaluator::initialize( Object &_object )
  {
    _object.bitstring.clear(); 
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom = structure.atoms.begin();
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom_end = structure.atoms.end();
    types::t_int concx = 0;
    types::t_int concy = 0;
    for(; i_atom != i_atom_end; ++i_atom )
    {
      bool flip = rng.flip();
      if ( i_atom->freeze & (i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
        flip = ( i_atom->type > 0 );
      _object.bitstring.push_back( flip ? 1.0: -1.0 );
      flip ? ++concx: --concx;

      ++i_atom;
      flip = rng.flip();
      if ( i_atom->freeze & (i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
        flip = ( i_atom->type > 0 );
      _object.bitstring.push_back( flip ? 1.0: -1.0 );
      flip ? ++concy: --concy;
    }
    if ( x_vs_y.is_singlec() )
    {
      types::t_unsigned N = structure.atoms.size() >> 1; 
      types::t_real xto_change = (types::t_real) N * x  - concx;
      types::t_real yto_change = (types::t_real) N * y  - concy;
      if (      xto_change > -1.0 and xto_change < 1.0 
           and  yto_change > -1.0 and xto_change < 1.0 ) return true;
      do
      {
        types::t_unsigned i = 2 * rng.random(N-1);
        if ( xto_change > 1.0 and _object.bitstring[i] < 0 )
          { _object.bitstring[i] = 1; xto_change-=2; }
        else if ( xto_change < -1.0 and _object.bitstring[i] > 0 )
          { _object.bitstring[i] = -1; xto_change+=2; }
        
        if ( yto_change > -1.0 and yto_change < 1.0 ) continue;
        i = 2 * rng.random(N-1) + 1;
        if ( yto_change > 1.0 and _object.bitstring[i] < 0 )
          { _object.bitstring[i] = 1; yto_change-=2; }
        else if ( yto_change < -1.0 and _object.bitstring[i] > 0 )
          { _object.bitstring[i] = -1; yto_change+=2; }

      } while (    xto_change < -1.0 or xto_change > 1.0
                or yto_change < -1.0 or yto_change > 1.0 );
    }
    return true;
  }

  void* const Evaluator::init( Object &_object )
  {
    structure << _object;
    functional.set_variables( &_object.bitstring );
    return &functional;
  }

  void* const Evaluator::LoadMinimizer(const TiXmlElement &_el )
  {
    std::string name = _el.Value();
    if ( name.compare("Minimizer") )
      return NULL;

    if ( not _el.Attribute("type") )
      return false;
    
    name = _el.Attribute("type");

    if ( name.compare("VA") == 0 )
    {
      darwin::printxmg.add_comment("VA optimizer");
      // pointer is owned by caller !!
      // don't deallocate
      return new minimizer::VA<t_Functional>( _el );
    }
    else if ( name.compare("SA") == 0 )
    {
      darwin::printxmg.add_comment("SA optimizer");
      // pointer is owned by caller !!
      // don't deallocate
      return new minimizer::VA<t_Functional>( _el );
    }
    else if ( name.compare("Beratan") == 0 )
    {
      darwin::printxmg.add_comment("SA optimizer");
      // pointer is owned by caller !!
      // don't deallocate
      return new minimizer::Beratan<t_Functional>( _el );
    }
    

    return NULL;
  }

  void Evaluator :: get_xy_concentrations( const Ising_CE::Structure &_str)
  {
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom = structure.atoms.begin();
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom_end = structure.atoms.end();
    types :: t_unsigned Nx = 0, Ny = 0;
    for(; i_atom != i_atom_end; ++i_atom )
      ( i_atom->site > 0 ) ? ++Ny, y += i_atom->type: 
                             ++Nx, x += i_atom->type;
    x /= (types::t_real) Nx;
    y /= (types::t_real) Ny;
  }

  bool Evaluator :: consistency_check()
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
    if (    ( lattice.sites[0].freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) 
         or ( lattice.sites[1].freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
    {
      get_xy_concentrations( structure );
      x_vs_y.set_xy( x, y );
      std::ostringstream sstr;
      sstr << " Setting Concentrations to x=" << std::fixed << std::setprecision(3 ) << x 
           << " and y=" << y;
      darwin::printxmg.add_comment( sstr.str() );
      if (    std::abs(x - x_vs_y.get_x(y)) > types::tolerance 
           or std::abs(y - x_vs_y.get_y(x)) > types::tolerance )
        darwin::printxmg.add_comment( " WARNING: x and y pair are strain mismatched!! " );
    }

    return true;
  }


  bool Evaluator :: Continue()
  {
    // on first iteration, writes band edges... then read them on following iterations
    ++age;
    functional.read_band_edges();

#ifdef _MPI
    if ( mpi::main.rank() != mpi::ROOT_NODE )
      return true;
#endif

    if ( check_ref_every != -1 ) // recomputes all electron band_edges every so often, if required
    {
      if ( not ( age % check_ref_every ) )
        functional.set_all_electron();
    }
    return true;
  }

  bool Functional::Load( const TiXmlElement &_node )
  { 
    if ( not vff.Load( _node ) )
    {
      std::cerr << " Could not load vff input!! " << std::endl; 
      return false;
    }
    if ( not vff.initialize_centers() )
    {
      std::cerr << " Could not initialize Atomic_Center list in vff!! " << std::endl
                << " Are you sure the lattice and the structure correspond? " << std::endl; 
      return false;
    }
    if ( not vff_minimizer.Load( _node ) )
    {
      std::cerr << " Could not load vff minimizer from input!! " << std::endl;
      return false;
    }
    pescan.set_band_edges( -666.666, 666.666 ); // wrong order! just to check whether Refs are read from input
    if ( not pescan.Load( _node ) )
    {
      std::cerr << " Could not load pescan interface from input!! " << std::endl; 
      return false;
    }
    types::t_real x, y;
    pescan.get_band_edges(x,y); // band edges have not been read if below is true
    if (     std::abs(x + 666.666 ) < types::tolerance 
         and std::abs(x - 666.666 ) < types::tolerance )
      set_all_electron();

    return true;
  }

  void Functional::write_band_edges()
  {
#ifdef _MPI 
    if ( mpi::main.rank() != mpi::ROOT_NODE )
      return;
#endif
    std::ofstream file( band_edge_filename.c_str(), std::ios_base::out | std::ios_base::trunc ); 
    if ( not file.is_open() )
      return;
    types :: t_real a, b;
    pescan.get_band_edges( a, b );
    file << a; if ( file.fail() ) return;
    file << b; if ( file.fail() ) return;
    file.close();
    return;
  }
  void Functional::read_band_edges()
  {
    types :: t_real a, b;
#ifdef _MPI 
    if ( mpi::main.rank() == mpi::ROOT_NODE )
    {
#endif
      std::ifstream file( band_edge_filename.c_str(), std::ios_base::in ); 
      if ( not file.is_open() )
        return;
      file >> a; if ( file.fail() ) return;
      file >> b; if ( file.fail() ) return;
      if ( a >= b )
       return;
      file.close();
#ifdef _MPI
    }
    mpi::BroadCast broadcast(mpi::main);
    broadcast.serialize(a);
    broadcast.serialize(b);
    broadcast.allocate_buffers();
    broadcast.serialize(a);
    broadcast.serialize(b);
    broadcast();
    broadcast.serialize(a);
    broadcast.serialize(b);
#endif 

    pescan.set_band_edges( a, b ); 
    return;
  }

  void rearrange_structure( Ising_CE::Structure &_str)
  {
    if ( not _str.lattice and _str.lattice->sites.size() != 2)
      return;

    std::list< Ising_CE::Structure::t_Atom > sites0;
    std::list< Ising_CE::Structure::t_Atom > sites1;
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom_end = _str.atoms.end();
    for(; i_atom != i_atom_end; ++i_atom )
      ( i_atom->site == 0 ) ?
        sites0.push_back( *i_atom ): sites1.push_back( *i_atom );

    std::list< Ising_CE::Structure::t_Atom > :: iterator i_0 = sites0.begin();
    std::list< Ising_CE::Structure::t_Atom > :: iterator i_end = sites0.end();
    std::list< Ising_CE::Structure::t_Atom > :: iterator i_1;
    atat::rVector3d translation = _str.lattice->sites[1].pos - _str.lattice->sites[0].pos; 
    types::t_real (*ptr_norm)(const atat::FixedVector<types::t_real, 3> &) = &atat::norm2;
    _str.atoms.clear();
    for(; i_0 != i_end; ++i_0 )
    {
      atat::rVector3d atom = i_0->pos + translation;
      i_1 = std::min_element ( sites1.begin(), sites1.end(), 
                      opt::ref_compose2( std::less<types::t_real>(),
                                         compose1( std::ptr_fun(ptr_norm),
                                                   bind2nd(std::minus<atat::rVector3d>(), atom) ),
                                         compose1( std::ptr_fun(ptr_norm),
                                                      bind2nd(std::minus<atat::rVector3d>(), atom) ) ) );
      _str.atoms.push_back( *i_0 ); 
      _str.atoms.push_back( *i_1 ); 
    }

  }



  //  sets concentration from k-space values.
  //  individual need not be physical (ie S_i=+/-1) when fourrier transformed to real-space
  //  a normalization procedure is applied which takes into account:
  //  (i) that this ga run is at set concentration (x =x0, y=y0)
  //  (ii) that x and y are only constrained by load-balancing
  void Evaluator :: set_concentration( Ising_CE::Structure &_str )
  {
    types::t_unsigned N = (types::t_int) _str.atoms.size(); N = N>>1;
    types::t_complex  *hold = new types::t_complex[ N ];
    if ( not hold )
    {
      std::cerr << " Could not allocate memory in set_concentration!! " << std::endl; 
      exit(0);
    }

    // creates individual with unnormalized occupation
    types::t_complex  *i_hold = hold;
    BandGap::fourrier_to_rspace( _str.atoms.begin(), _str.atoms.end(),
                                 _str.k_vecs.begin(), _str.k_vecs.end(),
                                 i_hold );

    Ising_CE::Structure::t_Atoms::iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms::iterator i_atom_end = _str.atoms.end();
    types :: t_int concx = 0, concy = 0; 
    i_hold = hold;
    for (; i_atom != i_atom_end; ++i_atom, ++i_hold)
    {
      if ( not i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T )
        i_atom->type = std::real(*i_hold);
      ( i_atom->type > 0 ) ? ++concx : --concx;

      ++i_atom;
      if ( not i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T )
        i_atom->type = std::imag(*i_hold);
      ( i_atom->type > 0 ) ? ++concy : --concy;
    }

    // then normalize it while setting correct concentrations
    if ( x_vs_y.is_singlec() )
    {
      normalize( _str, 0, (types::t_real) concx - ( (types::t_real) N ) * x_vs_y.get_x(x) );
      normalize( _str, 1, (types::t_real) concy - ( (types::t_real) N ) * x_vs_y.get_y(y) );
      delete[] hold;
      return;
    }
    
    // Concentration is not set, but still constrained by load balancing
    // hence will randomly set concentration to ( x and load balanced y )
    // or ( y and load balanced x). The latter case may not always be possible. 
    x = (double) concx / (double) N;
    if ( rng.flip() or not x_vs_y.can_inverse(x) )
    {  // x and load balanced y
      y = (double) concy / (double) N;
      x = x_vs_y.get_x( y );
      normalize( _str, 1, 0 );
      normalize( _str, 0, (types::t_real) concx - ( (types::t_real) N ) * x );
      delete[] hold;
      return;
    }
     
    // y and load balanced x
    y = x_vs_y.get_y( x );
    normalize( _str, 0, 0 );
    normalize( _str, 1, (types::t_real) concy - ( (types::t_real) N ) * y );
    delete[] hold;
  }


  // Takes an "unphysical" individual and set normalizes its sites _sites to +/-1,
  // after flipping the _tochange spins closest to zero.
  // ie sets the concentration
  void Evaluator :: normalize( Ising_CE::Structure &_str, 
                               const types::t_int _site, types::t_real _tochange) 
  {
    Ising_CE::Structure::t_Atoms::iterator i_end = _str.atoms.end();
    Ising_CE::Structure::t_Atoms::iterator i_which;
    Ising_CE::Structure::t_Atoms::iterator i_atom;
    while( _tochange < -1.0 or _tochange > 1.0 )
    {
      // first finds atom with type closest to zero from +1 or -1 side,
      // depending on _tochange
      i_atom = _str.atoms.begin();
      i_which = i_end;
      types::t_real minr = 0.0;
      for(; i_atom != i_end; ++i_atom )
      {
        if ( _site ) ++i_atom; 
        if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T )
          goto endofloop;
        if ( _tochange > 0 )
        {
          if ( i_atom->type < 0 )
            goto endofloop;
          if ( minr != 0.0 and i_atom->type > minr )
            goto endofloop;
        }
        else // ( _tochange < 0 )
        {
          if ( i_atom->type > 0 )
            goto endofloop;
          if ( minr != 0.0 and i_atom->type < minr )
            goto endofloop;
        }

        i_which = i_atom;
        minr = i_atom->type;

endofloop: 
        if ( not _site ) ++i_atom;
      }
      if ( i_which == i_end )
        throw std::runtime_error( "Error while normalizing x constituents\n" );

      i_which->type = ( _tochange > 0 ) ? -1.0: 1.0;
      _tochange += ( _tochange > 0 ) ? -2: 2;
    }

    // finally, normalizes _str
    i_atom = _str.atoms.begin();
    for(; i_atom != i_end; ++i_atom )
    {
      if ( _site ) ++i_atom;
      i_atom->type = ( i_atom->type > 0 ) ? 1.0: -1.0;
      if ( not _site ) ++i_atom;
    }
#ifdef _DEBUG
    types::t_real concx = 0;
    types::t_real concy = 0;
    types :: t_unsigned N = _str.atoms.size() >> 1;
    i_atom = _str.atoms.begin();
    for(; i_atom != i_end; ++i_atom )
    {
      i_atom->type > 0 ? ++concx: --concx;
      ++i_atom;
      i_atom->type > 0 ? ++concy: --concy;
    }
    types::t_real result =  _site ? (types::t_real ) concy / (types::t_real) N:
                                    (types::t_real ) concx / (types::t_real) N;
    types::t_real inv = 2.0 / (types::t_real) N;
    if ( std::abs( result - (_site ? y:x) ) > inv )
      throw std::runtime_error("Could not normalize site\n" );
#endif
  }

} // namespace pescan



#ifdef _MPI
namespace mpi
{
  template<>
  bool mpi::BroadCast::serialize<BandGap::Functional>( BandGap::Functional & _func )
  {
    if( not serialize( _func.vff ) ) return false;
    if( not _func.vff_minimizer.broadcast( *this ) ) return false;
    if( not serialize( _func.pescan ) ) return false;

    if ( stage == COPYING_FROM_HERE and rank() != mpi::ROOT_NODE )
      _func.vff.initialize_centers();

    return true;
  }
  template<>
  bool mpi::BroadCast::serialize<BandGap::Evaluator>( BandGap::Evaluator & _ev )
  {
    if( not serialize( _ev.lattice ) ) return false;
    if( not serialize( _ev.structure ) ) return false;
    if( not serialize( _ev.functional ) ) return false;
    if( not serialize( _ev.crossover_probability ) ) return false;
    if( not serialize( _ev.x ) ) return false;
    if( not serialize( _ev.y ) ) return false;
    if( not serialize( _ev.lessthan ) ) return false;
    if( not serialize( _ev.morethan ) ) return false;
    return serialize( _ev.x_vs_y );
  }
  template<>
  bool mpi::BroadCast::serialize<BandGap::Object>( BandGap::Object & _object )
  {
    return serialize( _object.bitstring );
  }
}
#endif
