#include <functional>
#include <algorithm>
#include <ext/algorithm>
#include <fstream>

#include "pescan.h"
#include "print_xmgrace.h"
#include <lamarck/atom.h>
#include <opt/va_minimizer.h>

#ifndef __PGI
  #include<ext/functional>
  using __gnu_cxx::compose1;
#else
  #include<functional>
  using std::compose1;
#endif
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
    structure.lattice = &lattice;
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
    if ( single_concentration )
      return false;
    structure << _object;
    get_xy_concentrations( structure );
    return x > lessthen or x < morethen;
  }
  
  eoMonOp<const Object>* Evaluator :: LoadTaboo(const TiXmlElement &_el )
  {
    if ( single_concentration )
      return NULL;
    const TiXmlElement *child = _el.FirstChildElement( "Concentration" );
    if ( not child )
      return NULL;
    double d;
    if ( child->Attribute( "lessthen" ) )
      child->Attribute( "xlessthan", &d );
    lessthen = ( d > 0 and d < 1 ) ? 2.0*d-1.0: 1.0;
    if ( child->Attribute( "morethen" ) )
      child->Attribute( "morethen", &d );
    morethen = ( d > 0 and d < 1 ) ? 2.0*d-1.0: -1.0;
    if ( lessthen < morethen )
      return NULL;
   
    std::ostringstream sstr;
    sstr << std::fixed << std::setprecision(3) << "Taboo x in [ " << 0.5*(morethen+1.0)
         << ", "  << 0.5*(lessthen+1.0) << "] ";
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
    types::t_int result = 0;
    for(; i_atom != i_atom_end; ++i_atom )
    {
      bool flip = rng.flip();
      _object.bitstring.push_back( flip ? 1.0: -1.0 );
      flip ? ++result: --result;
    }
    if ( single_concentration )
    {
      types::t_unsigned N = structure.atoms.size();
      types::t_int xto_change = static_cast<types::t_int>( (double) N * x ) - result;
      types::t_int yto_change = static_cast<types::t_int>( (double) N * y ) - result;
      if ( not ( xto_change or yto_change ) ) return true;
      do
      {

        types::t_unsigned i = rng.random(N-1);
        if ( i % 2 )
        {
          if ( xto_change > 0 and _object.bitstring[i] < 0 )
            { _object.bitstring[i] = 1; xto_change-=2; }
          else if ( xto_change < 0 and _object.bitstring[i] > 0 )
            { _object.bitstring[i] = -1; xto_change+=2; }
          continue;
        }
        if ( yto_change > 0 and _object.bitstring[i] < 0 )
          { _object.bitstring[i] = 1; yto_change-=2; }
        else if ( yto_change < 0 and _object.bitstring[i] > 0 )
          { _object.bitstring[i] = -1; yto_change+=2; }

      } while ( xto_change or yto_change );
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
      single_concentration = true; 
      get_xy_concentrations( structure );
      std::ostringstream sstr;
      sstr << " Setting Concentrations to x=" << std::fixed << std::setprecision(3 ) << x 
           << " and y=" << y;
      darwin::printxmg.add_comment( sstr.str() );
      if (    std::abs(x - x_vs_y.inverse(y)) > types::tolerance 
           or std::abs(y - x_vs_y.evaluate(x)) > types::tolerance )
        darwin::printxmg.add_comment( " WARNING: x and y pair are strain mismatched!! " );
    }

    return true;
  }

  void Evaluator :: set_concentration( Ising_CE::Structure &_str )
  {
    if ( single_concentration )
    {
      ::set_concentration( _str, 0, x );
      ::set_concentration( _str, 1, y );
      return;
    }
    
    if ( rng.flip() ) // set y concentration
    {
      x = ::set_concentration( _str, 0, -2 );
      y = x_vs_y.evaluate( x );
      ::set_concentration( _str, 1, y );
      return;
    }

    y = ::set_concentration( _str, 1, -2 );
    x = x_vs_y.inverse( y );
    ::set_concentration( _str, 0, x );
  }

  bool Evaluator :: Continue()
  {
    ++age;
    functional.read_band_edges();
    if ( check_ref_every != -1 and not ( age % check_ref_every ) )
      functional.set_all_electron();
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
    if ( not pescan.Load( _node ) )
    {
      std::cerr << " Could not load pescan interface from input!! " << std::endl; 
      return false;
    }
    return true;
  }

  void Functional::write_band_edges()
  {
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

    std::cout << "before rearranging" << std::endl;
    _str.print_out(std::cout); std::cout << std::endl;
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

    std::cout << "after rearranging" << std::endl;
    _str.print_out(std::cout); std::cout << std::endl;
  }



} // namespace CE

types::t_real set_concentration( Ising_CE::Structure &_str,
                                 types::t_int _site,
                                 types::t_real _target )
{
  types::t_unsigned N = (types::t_int) _str.atoms.size(); N = N>>1;
  types::t_complex  *hold = new types::t_complex[ N ];
  if ( not hold )
  {
    std::cerr << " Could not allocate memory in set_concentration!! " << std::endl; 
    exit(0);
  }
  types::t_complex  *i_hold = hold;
  BandGap::fourrier_to_rspace( _str.atoms.begin(), _str.atoms.end(),
                               _str.k_vecs.begin(), _str.k_vecs.end(),
                               i_hold );

  Ising_CE::Structure::t_Atoms::iterator i_atom = _str.atoms.begin();
  Ising_CE::Structure::t_Atoms::iterator i_atom_end = _str.atoms.end();
  types :: t_int result = 0; 
  i_hold = hold;
  for (; i_atom != i_atom_end; ++i_atom, ++i_hold)
  {
    if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T 
         or ( _site != -1 and _site != i_atom->site ) ) 
      ( i_atom->type > 0 ) ? ++result : --result;
    else 
    {
      ( std::real(*i_hold) > 0 ) ? ++result : --result;
      i_atom->type = ( std::real(*i_hold) > 0 ) ? 1.0: -1.0;
    }
  }

  if ( _target == -2.0 )
  {
    delete[] hold;
    return (double) result  / (double) N;
  }

  types::t_int to_change = static_cast<types::t_int>( (double) N * _target ) - result;
  if (  not to_change ) 
  {
    delete[] hold;
    return _target;
  }

  i_atom = _str.atoms.begin();
  i_hold = hold;
  for (; i_atom != i_atom_end; ++i_atom, ++i_hold)
    if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T 
         or ( _site != -1 and _site != i_atom->site ) ) 
      *i_hold = to_change;
    else if ( to_change > 0 and std::real( *i_hold ) > 0 )
      *i_hold = to_change;
    else if ( to_change < 0 and std::real( *i_hold ) < 0 )
      *i_hold = to_change;

  i_atom = _str.atoms.begin();
  types::t_int which;
  do
  {
    i_hold = hold;
    which = -1;
    types::t_real maxr = 1.0;
    if ( to_change > 0 )
    {
      for( types::t_unsigned i = 0; i < N; ++i, ++i_hold ) 
      {
        if ( _site != -1 and _site != i_atom->site )
          continue;
        if (     ( maxr == 1.0 and std::real( *i_hold ) < 0 )
             or  (  std::real( *i_hold ) < 0 and maxr <  std::real( *i_hold ) ) )
        {
          maxr = std::real( *i_hold );
          which = i;
        }
      }
      if ( which == -1 )
      {
        std::cerr << "Couldn't find a site to flip !? " << std::endl;
        exit(0);
      }
      ( i_atom + which )->type = 1.0;
      *( hold + which ) = to_change;
      result+=2; to_change-=2;
      continue;
    }
    maxr = -1.0;
    for( types::t_unsigned i = 0; i < N; ++i, ++i_hold ) 
    {
      if ( _site != -1 and _site != i_atom->site )
        continue;
      if (     ( maxr == -1.0 and std::real( *i_hold ) > 0 )
           or  (  std::real( *i_hold ) > 0 and maxr >  std::real( *i_hold ) ) )
      {
        maxr = std::real( *i_hold );
        which = i;
      }
    }
    ( i_atom + which )->type = -1.0;
    *( hold + which ) = to_change;
    result-=2; to_change+=2;

  } while (to_change); 
  
  delete[] hold;
  return (double) result  / (double) N;
}


#ifdef _MPI
namespace mpi
{
  template<>
  bool mpi::BroadCast::serialize<BandGap::Functional>( BandGap::Functional & _func )
  {
    if( not serialize( _func.vff ) ) return false;
    if( not _func.vff_minimizer.broadcast( *this ) ) return false;
    if( not serialize( _func.pescan ) ) return false;

    if ( stage == COPYING_FROM_HERE and rank() != ROOT_NODE )
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
    if( not serialize( _ev.lessthen ) ) return false;
    if( not serialize( _ev.morethen ) ) return false;
    if( not serialize( _ev.x_vs_y ) ) return false;

    return serialize( _ev.single_concentration );
  }
  template<>
  bool mpi::BroadCast::serialize<BandGap::Object>( BandGap::Object & _object )
  {
    return serialize( _object.bitstring );
  }
}
#endif
