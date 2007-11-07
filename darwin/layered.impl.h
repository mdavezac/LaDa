//
//  Version: $Id$
//

#ifndef _LAYERED_IMPL_H_
#define _LAYERED_IMPL_H_

#include <print/stdout.h>
#include <print/manip.h>
#include <lamarck/structure.h>
#include <lamarck/atom.h>

namespace Layered
{
  template<types::t_unsigned _D> template<class T_R_IT, class T_K_IT>
  Fourier<_D> :: Fourier ( T_R_IT _rfirst, T_R_IT _rend,
                           T_K_IT _kfirst, T_K_IT _kend )
  {
    const std::complex<types::t_real>
       imath(0, -2*3.1415926535897932384626433832795028841971693993751058208);
    
    for (; _kfirst != _kend; ++_kfirst)
    {
      _kfirst->type = std::complex<types::t_real>(0);
      for(T_R_IT i_r( _rfirst ); i_r != _rend; i_r += _D )
      {
        _kfirst->type +=   exp( imath * ( i_r->pos[0] * _kfirst->pos[0] +
                                          i_r->pos[1] * _kfirst->pos[1] +
                                          i_r->pos[2] * _kfirst->pos[2] ) )
                         * i_r->type;
      }
    }
  }
  template<types::t_unsigned _D> template<class T_R_IT, class T_K_IT, class T_O_IT >
  Fourier<_D> :: Fourier ( T_R_IT _rfirst, T_R_IT _rend,
                           T_K_IT _kfirst, T_K_IT _kend,
                           T_O_IT _rout )
  {
    const std::complex<types::t_real>
       imath(0, 2*3.1415926535897932384626433832795028841971693993751058208);
    for (; _rfirst != _rend; _rfirst+=_D, ++_rout)
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


  template<class T_CONTAINER>
  void operator<<(Ising_CE::Structure &_str, const Object<T_CONTAINER> &_o)
  {
    typedef typename Object<T_CONTAINER>::const_iterator const_iterator;
    const_iterator i_var = _o.bitstring.begin();
    const_iterator i_end = _o.bitstring.end();
    Ising_CE::Structure :: t_Atoms :: iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure :: t_Atoms :: iterator i_atom_end = _str.atoms.end();
    types :: t_unsigned _D = _str.lattice->sites.size();
    for(; i_var != i_end and i_atom != i_atom_end; )
    {
      if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) 
        { i_atom += _D; continue; }
      
      for(types::t_unsigned i=0; i < _D; ++i, ++i_atom )
        i_atom->type = BitString::spin_up(*i_var) ? 1.0 : -1.0;
      ++i_var; 
    }
  }
  template<class T_CONTAINER>
  void operator<<(Object<T_CONTAINER> &_o, const Ising_CE::Structure &_c)
  {
    _o.bitstring.clear(); _o.bitstring.reserve( _c.atoms.size() );
    Ising_CE::Structure :: t_Atoms :: const_iterator i_atom = _c.atoms.begin();
    Ising_CE::Structure :: t_Atoms :: const_iterator i_end = _c.atoms.end();
    types :: t_unsigned _D = _c.lattice->sites.size();
    for(; i_atom != i_end; i_atom += _D )
      if ( not (i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T) )
        _o.bitstring.push_back( BitString::spin_up(i_atom->type) ? 1.0: -1.0 );
  }



  template <types::t_unsigned _D>
  void Concentration<_D> :: operator()( Ising_CE::Structure &_str )
  {
    if ( _str.atoms.size() % _D != 0 )
      throw std::runtime_error( "Number of atoms is not a multiple of _D\n");

    if ( not single_c ) return;

    types::t_complex  *hold = new types::t_complex[ N ];
    if ( not hold )
    {
      std::cerr << " Could not allocate memory in set_concentration!! " << std::endl; 
      exit(0);
    }

    // creates individual with unnormalized occupation
    types::t_complex  *i_hold = hold;
    Fourier<_D>( _str.atoms.begin(), _str.atoms.end(),
                 _str.k_vecs.begin(), _str.k_vecs.end(),
                 i_hold );

    Ising_CE::Structure::t_Atoms::iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms::iterator i_atom_end = _str.atoms.end();
    types :: t_int concx = 0;
    i_hold = hold;
    if( not single_c )
    {
      for (; i_atom != i_atom_end; ++i_atom, ++i_hold)
      {
        if (i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T )  continue;
      
        if ( std::abs( std::real(*i_hold) ) < types::tolerance )
              i_atom->type = rng.flip() ? 1.0: -1.0;
        else  i_atom->type = BitString::spin_up( std::real(*i_hold) ) ? 1.0: -1.0;

        for(types::t_unsigned i = 1; i < _D; ++i, ++i_atom )
          (i_atom + 1)->type = BitString::spin_up(i_atom->type) ? 1.0: -1.0;
      }
      return;
    }

    for (; i_atom != i_atom_end; i_atom+=_D, ++i_hold)
    {
      if ( not ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
        i_atom->type = std::real(*i_hold);
      BitString::spin_up(i_atom->type) ? ++concx : --concx;
    }

    // then normalize it while setting correct concentrations
    normalize( _str, (types::t_real) concx - ( (types::t_real) N ) * x0 );
    delete[] hold;
    return;
  }



  // Takes an "unphysical" individual and normalizes its sites _sites to +/-1,
  // after flipping the _tochange spins closest to zero.
  // ie sets the concentration
  template <types::t_unsigned _D>
  void Concentration<_D> :: normalize( Ising_CE::Structure &_str,
                                       types::t_real _tochange ) 
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
      for(; i_atom != i_end; i_atom += _D )
      {
        if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) continue;
        if ( _tochange > 0 )
        {
          if ( i_atom->type < 0 )                    continue;
          if ( minr != 0.0 and i_atom->type > minr ) continue;
        }
        else // ( _tochange < 0 )
        {
          if ( i_atom->type > 0 )                    continue;
          if ( minr != 0.0 and i_atom->type < minr ) continue;
        }

        i_which = i_atom;
        minr = i_atom->type;
      }
      if ( i_which == i_end )
        throw std::runtime_error( "Error while normalizing x constituents\n" );

      i_which->type = ( _tochange > 0 ) ? -1.0: 1.0;
      _tochange += ( _tochange > 0 ) ? -2: 2;
    }

    // finally, normalizes _str
    i_atom = _str.atoms.begin();
    for(unsigned i=0; i_atom != i_end; ++i_atom, ++i )
    {
      i_atom->type = ( i_atom->type > 0 ) ? 1.0: -1.0;
      for(types::t_unsigned i = 1; i < _D; ++i, ++i_atom )
        (i_atom+1)->type = BitString::spin_up(i_atom->type) ? 1.0: -1.0;
    }
  }

  
  template <types::t_unsigned _D> template<class T_CONT>
  void Concentration<_D> :: operator()( BitString::Object<T_CONT> &_obj )
  {
    if ( not single_c ) return;

    // computes concentrations first
    get( _obj );

    types::t_real to_change = (types::t_real) N * ( x0  - x );
    // inline with non-use of fzzy math below...
    if ( to_change > -1.0 and to_change < 1.0 ) return;

    // Puts the positions which can be changed into a list
    std::vector<types::t_unsigned> pos;
    typedef typename BitString :: Object<T_CONT> :: const_iterator const_iterator;
    typedef typename BitString :: Object<T_CONT> :: t_Type t_Type;
    const_iterator i_bit = _obj.Container().begin();
    const_iterator i_bit_end = _obj.Container().end();
    for(types::t_unsigned i=0; i_bit != i_bit_end; ++i_bit, ++i)
      if( to_change > 0.0 and BitString::spin_down(*i_bit)  )
        pos.push_back( i );
      else if( to_change < 0.0 and BitString::spin_down(*i_bit)  )
        pos.push_back( i );

    // shuffle directions
    types::t_unsigned (*ptr_to_rng)(const types::t_unsigned& )
        = &eo::random<types::t_unsigned>;
    std::vector<types::t_unsigned> :: iterator i_pos = pos.begin();
    std::vector<types::t_unsigned> :: iterator i_pos_end = pos.end();
    std::random_shuffle(i_pos, i_pos_end, ptr_to_rng );
    i_pos = pos.begin();
    for(; i_pos != i_pos_end; ++i_pos)
    {
      BitString::flip<typename T_CONT::value_type>(_obj.bitstring[*i_pos]);
      ( to_change > 0 ) ? to_change -= 2: to_change += 2;

      // Fuzzy math at this point could create infinit loop
      if ( to_change < -1.0 or to_change > 1.0 ) break;
    }
    if ( not (to_change < -1.0 or to_change > 1.0) )
      throw std::runtime_error("Could not set the concentration of an object!!\n");
  }



  template <types::t_unsigned _D>
  void Concentration<_D> :: get( const Ising_CE::Structure &_str)
  {
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom_end = _str.atoms.end();
    for(; i_atom != i_atom_end; i_atom += _D )
      x += BitString::spin_up(i_atom->type) ?  1.0: -1.0;
    x /= (types::t_real) N; 
  }
  template <types::t_unsigned _D> template<class T_CONT>
  void Concentration<_D> :: get( const BitString::Object<T_CONT> &_obj )
  {
    if ( not single_c ) return;

    // computes concentrations first
    typedef typename BitString :: Object<T_CONT> :: const_iterator const_iterator;
    const_iterator i_bit = _obj.bitstring.begin();
    const_iterator i_bit_end = _obj.bitstring.end();
    types::t_unsigned conc = 0;
    for(; i_bit != i_bit_end; ++i_bit )
      conc += BitString::spin_up(*i_bit) ? 1: -1;

    // add frozen bits
    conc += Nfreeze;

    // finally normalizes
    x = (types::t_real) conc / (types::t_real) N;
  }

  template <types::t_unsigned _D>
  void Concentration<_D> :: setfrozen( const Ising_CE::Structure &_str )
  {
    N = _str.atoms.size();
    if( N % _D )
      throw std::runtime_error("Number of atoms in structure is not a multiple of _D\n");

    N /= _D;

    Ising_CE::Structure::t_Atoms::const_iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms::const_iterator i_atom_end = _str.atoms.end();
    Nfreeze = 0;
    for(; i_atom != i_atom_end; i_atom += _D )
      if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T )
        Nfreeze += BitString::spin_up(i_atom->type) ? 1 : -1; 
  }

  template <types::t_unsigned _D>
  std::string Concentration<_D> :: print() const 
  {
    if ( not single_c ) return "Concentration Range";
    std::ostringstream sstr;
    sstr << "Single Concentration, x0 = " << x0;
    return sstr.str();
  }

  template <types::t_unsigned _D>
  void Concentration<_D> :: LoadAttribute ( const TiXmlAttribute &_att )
  {
    std::string name = _att.Name();
    if ( name.compare("x0") ) return;
    
    double d;
    d = _att.DoubleValue();
    if( d < 0 or d > 1 ) goto errorout;
    single_c = true;
    x0 = 2.0 * (types::t_real) d - 1.0;
errorout:
    std::cerr << "Error while reading concentration input\n";
  }

  template<class T_INDIVIDUAL>
  inline void Evaluator<T_INDIVIDUAL> :: init( t_Individual &_indiv )
  {
    t_Base :: init( _indiv );
    // sets structure to this object 
    structure << *current_object;
  }


  template<class T_INDIVIDUAL>
  bool Evaluator<T_INDIVIDUAL> :: Load( const TiXmlElement &_node )
  {
    if ( not lattice.Load( _node ) )
    {
      std::cerr << " Could not load lattice type from input!! " << std::endl; 
      return false;
    }
    const TiXmlElement *parent = _node.FirstChildElement("Structure");
    if ( not parent )
      throw std::runtime_error("No Structure tag on input ?!\n");

    bool loadedstructure = false;
    Ising_CE::Structure::lattice = &lattice;
    if ( Load_Structure( *parent ) )  return true;

    std::cerr << "Found attributes for constructing a layered structure...\n" 
              << "But something went wrong.\n"
              << "Continuing with standard load.\n"; 
    
    return structure.Load( *parent );
  }

  template<class T_INDIVIDUAL>
  bool Evaluator<T_INDIVIDUAL> :: Load_Structure( const TiXmlElement &_node )
  {
    if(     ( not _node.Attribute("direction") )
        and ( not _node.Attribute("multiplicity") ) ) return false;
    if(     ( not _node.Attribute("direction") )
        or  ( not _node.Attribute("multiplicity") ) )
    {
      std::cerr << "Either cell direction or multiplicity is missing on input\n";
      return false;
    }
    atat::rVector3d cdir;
    types :: t_unsigned multiplicity; 
    atat::rMatrix3d &cell = structure.cell;
    
    // First, Load Attributes 
    bool couldload = true;
    std::istringstream sstr; sstr.str( _node.Attribute("direction") );
    sstr >> direction[0]; if ( sstr.fail() ) couldload = false;
    sstr >> direction[1]; if ( sstr.fail() ) couldload = false;
    sstr >> direction[2]; if ( sstr.fail() ) couldload = false;

    if ( atat::norm2( direction ) < types::tolerance ) couldload = false;

    if ( not couldload )
    {
      std::cerr << "Error while trying to read direction of cell\n";
      return false; 
    }
    
    int u;
    _node.Attribute( "multiplicity", &u );
    if ( u < 0 ) 
    {
      std::cerr << "Cell multiplicity is not found or negative "
                << "input\n";
      return false;
    }
    multiplicity = (types::t_unsigned) std::abs( u );

    // Load in scale
    if( _node.Attribute("scale") )
      _node.Attribute("scale", &structure.scale);
    else if ( structure.lattice ) 
      structure.scale = structure.lattice->scale;
    else
      structure.scale = 2;

    // Load in PI and Energy
    structure.Pi_name = 0;
    if( _node.Attribute("PI") )
    {
      _node.Attribute("PI", &u);
      structure.Pi_name = (types::t_int) u;
    }
    double d;
    structure.energy = 666.666;
    if( _node.Attribute("energy") )
    {
      _node.Attribute("energy", &d);
      structure.energy = (types::t_real) d;
    }


    // Then constructs unit cell
    types::t_real a = (types::t_real) multiplicity;
    cell = lattice.cell;
    cell.set_column(0, lattice.cell * direction * a ); 

    // Checks that cell is not singular
    if ( std::abs( det(cell) ) < types::tolerance )
      cell.set_column(1, lattice.cell.get_column( 0 ) );
    if ( std::abs( det(cell) ) < types::tolerance )
      cell.set_column(2, lattice.cell.get_column( 1 ) );
    if ( std::abs( det(cell) ) < types::tolerance )
    {
      std::cerr << "Could not construct unit-cell\n" << cell << std::endl;
      return false;
    }

    // Makes sure the triad is direct
    if ( det(cell) < 0 )
    {
      atat::rVector3d d = cell.get_column(2);
      cell.set_column(2, cell.get_column(1) );
      cell.set_column(1, d);
    }

    // Fills in a structure lattice points.
    Ising_CE::Structure copy_structure = structure;
    FillStructure( copy_structure );
    
    // We now sort the real space atoms according to layer
    std::sort( copy_structure.atoms.begin(), copy_structure.atoms.end(),
               Depth( cell ) );
    // The lattice sites are also sorted the same way
    std::sort( lattice.sites.begin(), lattice.sites.end(),
               Depth( cell ) );
    Ising_CE::Lattice::t_Sites::iterator i_site = lattice.sites.begin(); 
    Ising_CE::Lattice::t_Sites::iterator i_site_end = lattice.sites.end();
    for( types::t_unsigned n=0; i_site != i_site_end; ++i_site, ++n )
      i_site->site = n;


    // Finally, we copy the kvector positions as atoms, and the related sites
    Ising_CE::Structure::t_Atoms::const_iterator i_vec = copy_structure.atoms.begin();
    Ising_CE::Structure::t_Atoms::const_iterator i_vec_end = copy_structure.atoms.end();
    structure.atoms.clear();
    structure.atoms.reserve( lattice.sites.size() * copy_structure.atoms.size() );
    bool only_one_site = lattice.sites.size() == 1;
    atat::rVector3d origin = lattice.sites.front().pos;
    for(; i_vec != i_vec_end; ++i_vec )
    {
      Ising_CE::Structure::t_Atom atom;
      atom.site = -1; atom.pos = i_vec->pos;
      atom.type = Ising_CE::Structure::t_Atom::t_Type(0);
      atom.freeze = lattice.sites.front().freeze;
      structure.atoms.push_back(atom);

      if( only_one_site ) continue;

      i_site = lattice.sites.begin(); 
      i_site_end = lattice.sites.end();
      for(++i_site; i_site != i_site_end; ++i_site )
      {
        Ising_CE::Structure::t_Atom atom2(atom);
        atom2.pos += ( i_site->pos - origin );
        atom2.site = i_site->site;
        // vewy impowtant. Otherwise GA::KRandom and GA::Random produce
        // individuals which are to big.
        atom2.freeze = i_site->freeze | Ising_CE::Structure::t_Atom::FREEZE_T;
        structure.atoms.push_back(atom2);
      }
    }
   
    // More of the same... but for kvecs
    structure.find_k_vectors();

    concentration.setfrozen( structure );

    return true;

  }

  template< class T_INDIVIDUAL >
  GA::Taboo_Base<T_INDIVIDUAL>*
     Evaluator<T_INDIVIDUAL> :: LoadTaboo(const TiXmlElement &_el )
  {
    std::string name = _el.Value();
    GA::Taboo_Base<t_Individual> *taboo = NULL;

    if( name == "Layer" ) 
      taboo = new Taboo< t_Individual >( structure );
    else if (     name == "Concentration" 
              and ( not concentration.is_single_c() ) )
      taboo = new GA::xTaboo< t_Individual >( concentration );


    if ( taboo and taboo->Load( _el ) )  return taboo;
    if ( taboo ) delete taboo;
    return NULL;
  }

  template< class T_INDIVIDUAL >
  inline bool Evaluator<T_INDIVIDUAL> :: initialize( t_Individual &_indiv )
  {
    GA::Random< t_Individual > random( concentration, structure, _indiv );
    _indiv.invalidate(); return true;
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
    t_FourierRtoK( s.atoms.begin(),  s.atoms.end(),
                   s.k_vecs.begin(), s.k_vecs.end() );
    s.print_xml(_node);

    return true;
  }


  inline bool Depth::operator()( const atat::rVector3d &_1,
                                 const atat::rVector3d &_2 )
  {
    types::t_real a =   _1 * a0;
    types::t_real b =   _2 * a0;
    if ( not opt::Fuzzy<types::t_real> :: equal( a, b ) )
      return opt::Fuzzy<types::t_real> :: less( a, b );

    a =   _1 * a1;
    b =   _2 * a1;
    if ( not opt::Fuzzy<types::t_real> :: equal( a, b ) )
      return opt::Fuzzy<types::t_real> :: less( a, b );
    
    a =   _1 * a2;
    b =   _2 * a2;
    return opt::Fuzzy<types::t_real> :: less( a, b );
  }





  template< class T_INDIVIDUAL > 
  Taboo<T_INDIVIDUAL> :: Taboo   ( const Ising_CE::Structure &_str )
                                         : t_Base(), d0(0), d1(0)
  {
    Ising_CE::Structure::t_Atoms::const_iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms::const_iterator i_atom_end = _str.atoms.end();
    const types::t_unsigned is_frozen = Ising_CE::Structure::t_Atom::FREEZE_T;
    const types::t_unsigned _d = _str.lattice->sites.size();
    for(; i_atom != i_atom_end; i_atom += _d )
      if ( i_atom->freeze & is_frozen ) 
        sites.push_back( BitString::spin_up( i_atom->type ) ? 1: -1 );
      else sites.push_back(0); 
    Nmax = _str.atoms.size() * 2;
  }


  template< class T_INDIVIDUAL >
  bool Taboo<T_INDIVIDUAL> :: Load ( const TiXmlElement &_node )
  {
    if( not Ising_CE::Structure::lattice ) return false;
    std::string name = _node.Value();
    if( name.compare("Layer") != 0 ) return false;

    const Ising_CE::Lattice& lattice  = *Ising_CE::Structure::lattice;
    const std::vector< std::string > &strings = lattice.sites.front().type;
    std::vector< std::string > :: const_iterator i_string = strings.begin();

    if( _node.Attribute( i_string->c_str() ) )
    {
       types::t_int d = 0;
       _node.Attribute( i_string->c_str(), &d );
       if ( d > 0 and d < Nmax ) d0 = std::abs( d ); 
    }

    ++i_string;
    //! Only one type...
    if ( i_string == strings.end() ) return false;

      
    if( _node.Attribute( i_string->c_str() ) )
    {
       types::t_int d = 0;
       _node.Attribute( i_string->c_str(), &d );
       if ( d > 0 and d < Nmax ) d1 = std::abs( d ); 
    }

    // Returns false if neither d0 nor d1 impose a constraint
    if( d0 == 0 ) d0 = Nmax;
    if( d1 == 0 ) d1 = Nmax;
    return not ( d0 == Nmax and d1 == Nmax ); 
  }


  template< class T_INDIVIDUAL >
  bool Taboo<T_INDIVIDUAL> :: operator() ( const t_Individual &_indiv ) const
  {
    std::vector<types::t_int> :: const_iterator i_site = sites.begin();
    std::vector<types::t_int> :: const_iterator i_site_end = sites.end();
    typedef typename t_Individual :: t_IndivTraits :: t_Object t_Object;
    typedef typename t_Object :: t_Container :: const_iterator const_iterator;
    const_iterator i_var = _indiv.Object().Container().begin();
    const_iterator i_var_end = _indiv.Object().Container().end();

    typedef std::pair<bool, types::t_unsigned> t_layerpair;
    //! Keeps track of the width and type of the current layer
    t_layerpair  current( BitString::spin_up( *i_site == 0 ? *i_var: *i_site ), 0 );
    //! Keeps a record of the first layer in order to do periodicity check
    t_layerpair  first(0,0);
    for(; i_site != i_site_end and i_var != i_var_end;  ++i_site )
    {
      bool type;
      type = BitString::spin_up( *i_site == 0 ? *i_var: *i_site );

      if ( current.first != type )
      {
        if ( first.first == 0 ) first = current;
        current.first = type;
        current.second = 0;
      }
      ++current.second;
      
      if ( current.second > ( current.first ? d0: d1 ) ) return true;
      if ( *i_site == 0 )  ++i_var;
    }

    if ( i_site != i_site_end or i_var != i_var_end )
      throw std::runtime_error( "Layered::Taboo and individivual's container do not match !?\n" );

    //! periodicicity: adds first to current (eg last) layer if are the same type 
    if( current.first == first.first ) current.second += first.second;
    return current.second > ( current.first ? d0: d1 );
  }

} // namespace Layered


#endif
