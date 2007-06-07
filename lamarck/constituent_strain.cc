#include<numeric>
#include<complex>
#include<iomanip>
#include<math.h>
#include<limits.h>

#ifndef __PGI
  #include<ext/functional>
  using __gnu_cxx::compose1;
#else
  #include<functional>
  using std::compose1;
#endif


#include <opt/compose_functors.h>
#include <opt/ndim_iterator.h>
#include <analysis/analyze_code.h>

#include "constituent_strain.h"

namespace Ising_CE 
{
  std::vector<Harmonic> Constituent_Strain :: harmonics;
  const types::t_real Constituent_Strain::ZERO_TOLERANCE = 0.0000000001;

  void Constituent_Strain :: set_structure(const Ising_CE::Structure& str, const atat::rMatrix3d &lattice)
  {
     atat::rVector3d kvec;
     atat::rMatrix3d k_lat = !( ~(lattice) );
     std::cout << std::fixed << std::right
               << std::showpos << std::setprecision(6) << std::setw(14);

     r_cell = str.cell;
     k_cell = !( ~(r_cell) );
     r_vecs.clear(); k_vecs.clear();

     // gets r_vecs from str
     {
       std::vector<Atom> :: const_iterator i_atom = str.atoms.begin();
       std::vector<Atom> :: const_iterator i_end = str.atoms.end();
       for ( ; i_atom != i_end; i_atom++ )
         r_vecs.push_back( i_atom->pos );
     }

     // creates k vector list
     {
       // A is the basis used to determine "a" first brillouin zone
       atat::rMatrix3d A = (!k_lat) * k_cell;
       atat::iVector3d range, min;
       
       // computes range up to first periodic image
       find_range( A, range );
       
       // sets up the n-dimensional iterators
       opt::Ndim_Iterator< types::t_int, std::less_equal<types::t_int> > global_iterator;
       global_iterator.add( 0, range[0]);
       global_iterator.add( 0, range[1]);
       global_iterator.add( 0, range[2]);
       
       // the following loop creates all possible k-vectors,
       // it then refolds them and adds them to the k vector list
       // only one copy of each refolded vector is allowed
       types::t_real (*ptr_norm)(const atat::FixedVector<types::t_real, 3> &) = &atat::norm2;
       std::vector<atat::rVector3d> :: iterator i_begin = k_vecs.begin();
       std::vector<atat::rVector3d> :: iterator i_end = k_vecs.end();
       std::vector<atat::rVector3d> :: iterator i_which;
       do
       {
         // creates vector in A basis
         kvec[0] =  (types::t_real) global_iterator.access(0);
         kvec[1] =  (types::t_real) global_iterator.access(1);
         kvec[2] =  (types::t_real) global_iterator.access(2);
         kvec = A * kvec;
       
         kvec[0] -= rint(kvec[0]); 
         kvec[1] -= rint(kvec[1]); 
         kvec[2] -= rint(kvec[2]); 
         
         // switches to cartesian coordinates
         kvec = k_lat * kvec;
//        refold( kvec, k_lat);
         
         // find if vector is already in list
         i_which = std::find_if( i_begin, i_end, 
                        compose1( std::bind2nd(std::less<types::t_real>(), ZERO_TOLERANCE),
                        compose1( std::ptr_fun(ptr_norm),
                                  bind2nd(std::minus<atat::rVector3d>(), kvec) ) ) );
         // if it isn't, adds it
         if ( i_which == i_end  ) 
         {
           k_vecs.push_back( kvec );
           i_begin = k_vecs.begin();
           i_end = k_vecs.end();
         }
       
       } while( ++global_iterator );

       // refolds the vectors somewhat better
       i_begin = k_vecs.begin();
       i_end = k_vecs.end();
       for( ; i_begin != i_end; i_begin++ )
         refold(*i_begin, k_lat);

       // finally, puts vector 0,0,0 at front of list
       i_begin = k_vecs.begin();
       i_end = k_vecs.end();
       i_which = std::find_if( i_begin, i_end, 
                      compose1( std::bind2nd(std::less<types::t_real>(), ZERO_TOLERANCE),
                      compose1( std::ptr_fun(ptr_norm), std::_Identity<atat::rVector3d>() ) ) );
 
       if ( i_which != i_end  ) 
         std::iter_swap( i_which, k_vecs.begin() );
         
       // the refolding is not perfect, we now remove equivalent
       // vectors "by hand "
       remove_equivalents(k_lat);

       std::sort( k_vecs.begin(), k_vecs.end(),
                  opt::ref_compose2( std::less<types::t_real>(), std::ptr_fun(ptr_norm),
                                     std::ptr_fun(ptr_norm) ) );

//      std::cout << k_vecs.size() << std::endl;
//      i_begin = k_vecs.begin();
//      i_end = k_vecs.end();
//      for( ; i_begin != i_end; i_begin++ )
//      {
//        std::cout << " -> " << ( *i_begin )
//                  << std::endl;
//     }

     }

  }

  void Constituent_Strain :: cut_integer_part (atat::rVector3d &kvec)
  {
    kvec[0] -= rint(kvec[0]); 
    kvec[1] -= rint(kvec[1]); 
    kvec[2] -= rint(kvec[2]); 
  }


  void  Constituent_Strain :: find_range( const atat::rMatrix3d &A,
                                          atat::iVector3d &kvec )
  {
    atat::rVector3d a = A.get_column(0), b;
    types::t_int n = 1;
    b = a;
    while( not is_int(b) )
      { b += a; n++;  }
    kvec[0] = n;

    a = A.get_column(1);
    b = a; n = 1;
    while( not is_int(b) )
      { b += a; n++;  }
    kvec[1] = n;
    
    a = A.get_column(2);
    b = a; n = 1;
    while( not is_int(b) )
      { b += a; n++;  }
    kvec[2] = n;
  }


  // loads coefficients from input into static member variable x_coefs
  bool Constituent_Strain :: Load_Harmonics( const TiXmlElement &_element)
  {
    const TiXmlElement *child;
    Harmonic harmonic;
    types::t_real d=1.0;

    if ( _element.Attribute("attenuation", &d) )
      Harmonic::set_attenuation( d );
    
    child = _element.FirstChildElement( "Harmonic" );
    if ( !child )
      return false;
    for ( ; child; child=child->NextSiblingElement( "Harmonic" ) )
    {
      harmonic.clear(); // clear harmonic
      if ( not harmonic.Load( *child ) )
        return false;
      harmonics.push_back( harmonic );
    }

    return true;
  }


  types::t_real Constituent_Strain :: evaluate()
  {
    START_ANALYSIS("Constituent_Strain :: evaluate")
    #ifdef _DEBUG_LADA_
      if( variables->size() == 0 ) 
      {
        std::cerr << "variables have not been initialized in "
                  << "types::t_real Constituent_Strain :: evaluate()"
                  << std::endl;
        exit(0);
      }
    #endif // _DEBUG_LADA_
    types::t_real sum_harm;
    std::complex<types::t_real> sum_exp;
    const std::complex<types::t_real> imath(0, -2*3.1415926535897932384626433832795028841971693993751058208);
    types::t_real x = get_concentration();

    
    std::vector<atat::rVector3d> :: const_iterator i_k_vec = k_vecs.begin();
    std::vector<atat::rVector3d> :: const_iterator i_k_vec_end = k_vecs.end();
    std::vector<atat::rVector3d> :: const_iterator i_r_vec;
    std::vector<atat::rVector3d> :: const_iterator i_r_vec_begin = r_vecs.begin();
    std::vector<atat::rVector3d> :: const_iterator i_r_vec_end = r_vecs.end();
    std::vector<Harmonic> :: const_iterator i_harmonic;
    std::vector<Harmonic> :: const_iterator i_harmonic_begin = harmonics.begin();
    std::vector<Harmonic> :: const_iterator i_harmonic_end = harmonics.end();
    std::vector<types::t_real> :: const_iterator i_spin;
    std::vector<types::t_real> :: const_iterator i_spin_begin = variables->begin();
    types::t_real *interpolation = new types::t_real[ harmonics.size() ];
    types::t_real *i_inter;
    if (! interpolation )
    {
      std::cerr << "Memory allocation error in "
                << "Constituent_Strain :: evaluate(...) "
                << std::endl;
      exit(0);
    }
    i_harmonic = i_harmonic_begin;
    i_inter = interpolation;
    for ( ; i_harmonic != i_harmonic_end; ++i_harmonic, ++i_inter )
      *i_inter = i_harmonic->evaluate(x);

    types::t_real value = 0.0;
    for ( ; i_k_vec != i_k_vec_end; ++i_k_vec )
    {
      if ( norm2( *i_k_vec ) < ZERO_TOLERANCE ) // don't need to compute Gamma
        ++i_k_vec;

      i_harmonic = i_harmonic_begin;
      sum_harm = 0.0;
      i_inter = interpolation;
      for ( ; i_harmonic != i_harmonic_end; ++i_harmonic, ++i_inter )
        sum_harm += i_harmonic->evaluate(*i_k_vec) * (*i_inter);
      
      i_r_vec = i_r_vec_begin; i_spin = i_spin_begin;
      sum_exp = 0.0;
      for ( ; i_r_vec != i_r_vec_end; ++i_r_vec, ++i_spin )
      {

        sum_exp +=   exp( imath * ( (*i_r_vec)[0] * (*i_k_vec)[0] +
                                    (*i_r_vec)[1] * (*i_k_vec)[1] +
                                    (*i_r_vec)[2] * (*i_k_vec)[2] ) )
                   * (*i_spin);
      }
      value +=  (real(sum_exp * conj( sum_exp ))) * sum_harm;
    }
    value /= ( (types::t_real)  k_vecs.size() ) * (types::t_real) r_vecs.size();

    delete[] interpolation;
    END_ANALYSIS;
    return value;
  }

  types::t_real Constituent_Strain :: evaluate_one_gradient( types::t_unsigned _pos ) 
  {
    START_ANALYSIS("Constituent_Strain :: evaluate_gradient")
    #ifdef _DEBUG_LADA_
      if( variables->size() == 0 ) 
      {
        std::cerr << "variables have not been initialized in "
                  << "types::t_real Constituent_Strain :: evaluate_with_gradient( types::t_real *)"
                  << std::endl;
        exit(0);
      }
    #endif // _DEBUG_LADA_
    std::vector<atat::rVector3d> :: const_iterator i_k_vec = k_vecs.begin();
    std::vector<atat::rVector3d> :: const_iterator i_k_vec_end = k_vecs.end();
    std::vector<atat::rVector3d> :: const_iterator i_r_vec;
    std::vector<atat::rVector3d> :: const_iterator i_r_vec_begin = r_vecs.begin();
    std::vector<atat::rVector3d> :: const_iterator i_r_vec_end = r_vecs.end();
    std::vector<Harmonic> :: const_iterator i_harmonic;
    std::vector<Harmonic> :: const_iterator i_harmonic_begin = harmonics.begin();
    std::vector<Harmonic> :: const_iterator i_harmonic_end = harmonics.end();
    std::vector<types::t_real> :: const_iterator i_spin;
    std::vector<types::t_real> :: const_iterator i_spin_begin = variables->begin();
    std::vector<types::t_real> :: const_iterator i_spin_end = variables->end();
    types::t_real sum_harm, deriv_harm, result = 0.0;
    std::complex<types::t_real> sum_exp, c_h;
    const std::complex<types::t_real> imath(0, -2*3.1415926535897932384626433832795028841971693993751058208);
    types::t_real x = get_concentration();
    types::t_unsigned N= r_vecs.size();
    std::complex<types::t_real> exponential;
    types::t_real inv_N = 1/( (types::t_real)N  * (types::t_real) k_vecs.size()), dummy;

    // computes interpolation stuff and derivatives
    types::t_real *interpolation = new types::t_real[ 2*harmonics.size() ];
    types::t_real *i_inter, *ig_inter;
    if (! interpolation )
    {
      std::cerr << "Memory allocation error in "
                << "Constituent_Strain :: evaluate(...) "
                << std::endl;
      exit(0);
    }
    types::t_real *g_interpolation = interpolation + harmonics.size();
    
    i_harmonic = i_harmonic_begin;
    i_inter = interpolation;
    ig_inter = g_interpolation;
    for ( ; i_harmonic != i_harmonic_end; ++i_harmonic, ++i_inter, ++ig_inter )
    {
        *i_inter = i_harmonic->evaluate_with_gradient(x, *ig_inter);
        *ig_inter /= (types::t_real) N;
    }


    for ( ; i_k_vec != i_k_vec_end; ++i_k_vec )
    {
      if ( norm2( *i_k_vec ) < ZERO_TOLERANCE ) // don't need to compute Gamma
        ++i_k_vec; 

      i_harmonic = i_harmonic_begin;
      i_inter = interpolation;
      ig_inter = g_interpolation;
      sum_harm = 0.0; deriv_harm = 0.0;
      for ( ; i_harmonic != i_harmonic_end; ++i_harmonic, ++i_inter, ++ig_inter )
      {
        dummy = i_harmonic->evaluate(*i_k_vec);
        sum_harm += dummy * (*i_inter);
        deriv_harm += dummy * (*ig_inter);
      }
      
      i_r_vec = i_r_vec_begin;
      i_spin = i_spin_begin;
      sum_exp = 0.0;
      for ( types::t_unsigned n=0; i_r_vec != i_r_vec_end; ++i_r_vec,  ++i_spin, ++n )
      {
        std::complex<types::t_real> store;
        store = exp( imath * ( (*i_r_vec)[0] * (*i_k_vec)[0] +
                               (*i_r_vec)[1] * (*i_k_vec)[1] +
                               (*i_r_vec)[2] * (*i_k_vec)[2] ) );
        sum_exp +=  (*i_spin) * (store);
        if ( n == _pos )
          exponential = store;
      }

      result +=   real( sum_exp  * conj( sum_exp ) ) * deriv_harm 
                + real( exponential * conj( sum_exp ) ) * sum_harm * 2.0;

    }

    delete[] interpolation;
    END_ANALYSIS;
    return result * inv_N;
  }
  types::t_real Constituent_Strain :: evaluate_with_gradient( types::t_real* const gradient ) 
  {
    START_ANALYSIS("Constituent_Strain :: evaluate_with_gradient")
    #ifdef _DEBUG_LADA_
      if( variables->size() == 0 ) 
      {
        std::cerr << "variables have not been initialized in "
                  << "types::t_real Constituent_Strain :: evaluate_with_gradient( types::t_real *)"
                  << std::endl;
        exit(0);
      }
    #endif // _DEBUG_LADA_
    std::vector<atat::rVector3d> :: const_iterator i_k_vec = k_vecs.begin();
    std::vector<atat::rVector3d> :: const_iterator i_k_vec_end = k_vecs.end();
    std::vector<atat::rVector3d> :: const_iterator i_r_vec;
    std::vector<atat::rVector3d> :: const_iterator i_r_vec_begin = r_vecs.begin();
    std::vector<atat::rVector3d> :: const_iterator i_r_vec_end = r_vecs.end();
    std::vector<Harmonic> :: const_iterator i_harmonic;
    std::vector<Harmonic> :: const_iterator i_harmonic_begin = harmonics.begin();
    std::vector<Harmonic> :: const_iterator i_harmonic_end = harmonics.end();
    std::vector<types::t_real> :: const_iterator i_spin;
    std::vector<types::t_real> :: const_iterator i_spin_begin = variables->begin();
    std::vector<types::t_real> :: const_iterator i_spin_end = variables->end();
    types::t_real sum_harm, deriv_harm, *grad;
    std::complex<types::t_real> sum_exp, c_h;
    const std::complex<types::t_real> imath(0, -2*3.1415926535897932384626433832795028841971693993751058208);
    types::t_real x = get_concentration();
    types::t_unsigned N= r_vecs.size();
    std::complex<types::t_real> *exp_begin = new std::complex<types::t_real> [ N ];
    std::complex<types::t_real> *i_exp_end = exp_begin + r_vecs.size(), *i_exp;
    types::t_real inv_N = 1/( (types::t_real)N  * (types::t_real) k_vecs.size()), dummy;
    if ( not exp_begin )
    {
      std::cerr << "Memory allocation error in "
                << "Constituent_Strain::evaluate_with_gradient(..)"
                << std::endl;
      exit(0);
    }

    // computes interpolation stuff and derivatives
    types::t_real *interpolation = new types::t_real[ 2*harmonics.size() ];
    types::t_real *i_inter, *ig_inter;
    if (! interpolation )
    {
      std::cerr << "Memory allocation error in "
                << "Constituent_Strain :: evaluate(...) "
                << std::endl;
      exit(0);
    }
    types::t_real *g_interpolation = interpolation + harmonics.size();
    
    i_harmonic = i_harmonic_begin;
    i_inter = interpolation;
    ig_inter = g_interpolation;
    for ( ; i_harmonic != i_harmonic_end; ++i_harmonic, ++i_inter, ++ig_inter )
    {
        *i_inter = i_harmonic->evaluate_with_gradient(x, *ig_inter);
        *ig_inter /= (types::t_real) N;
    }

    grad = gradient;
    for (types::t_unsigned i=0; i< variables->size(); ++grad, ++i)
      *grad *= 0.0;

    types::t_real value = 0.0;
    for ( ; i_k_vec != i_k_vec_end; ++i_k_vec )
    {
      if ( norm2( *i_k_vec ) < ZERO_TOLERANCE ) // don't need to compute Gamma
        ++i_k_vec; 

      i_harmonic = i_harmonic_begin;
      i_inter = interpolation;
      ig_inter = g_interpolation;
      sum_harm = 0.0; deriv_harm = 0.0;
      for ( ; i_harmonic != i_harmonic_end; ++i_harmonic, ++i_inter, ++ig_inter )
      {
        dummy = i_harmonic->evaluate(*i_k_vec);
        sum_harm += dummy * (*i_inter);
        deriv_harm += dummy * (*ig_inter);
      }
      
      i_exp = exp_begin;
      i_r_vec = i_r_vec_begin;
      i_spin = i_spin_begin;
      sum_exp = 0.0;
      for ( ; i_r_vec != i_r_vec_end; ++i_r_vec, ++i_exp, ++i_spin )
      {
        *i_exp = exp( imath * ( (*i_r_vec)[0] * (*i_k_vec)[0] +
                                (*i_r_vec)[1] * (*i_k_vec)[1] +
                                (*i_r_vec)[2] * (*i_k_vec)[2] ) );
        sum_exp +=  (*i_spin) * (*i_exp);
      }

      grad = gradient;
      i_exp = exp_begin;
      for ( ; i_exp != i_exp_end; ++grad, ++i_exp)
        *grad +=   real( sum_exp  * conj( sum_exp ) ) * deriv_harm 
                 + real( (*i_exp) * conj( sum_exp ) ) * sum_harm * 2.0;

      value += real(sum_exp * conj( sum_exp )) * sum_harm;
    }
    value *= inv_N;
    grad = gradient;
    for (types::t_unsigned i=0; i< variables->size(); ++grad, ++i)
      *grad *= inv_N;

    delete[] exp_begin;
    delete[] interpolation;

    END_ANALYSIS;
    return value;
  }

  void Constituent_Strain :: remove_equivalents( const atat::rMatrix3d &lat)
  {
    std::vector<atat::rVector3d> :: iterator i_vec = k_vecs.begin();
    std::vector<atat::rVector3d> :: iterator i_end = k_vecs.end();
    std::vector<atat::rVector3d> :: iterator i_which;

    while( i_vec != i_end )
    {
      i_which = i_vec+1;
      for ( ; i_which != i_end; i_which++ )
        if ( are_equivalent( *i_which, *i_vec, lat ) )
          break;

      if ( i_which != i_end )
      {
        if ( norm2( *i_vec ) < norm2( *i_which ) )
          k_vecs.erase(i_which);
        else
          k_vecs.erase(i_vec);
        i_vec = k_vecs.begin();
        i_end = k_vecs.end();
      }
      else
        i_vec++;
    }

  }

  bool Constituent_Strain :: are_equivalent( const atat::rVector3d &vec_a,
                                             const atat::rVector3d &vec_b,
                                             const atat::rMatrix3d &lat) const
  {

    opt::Ndim_Iterator<types::t_int, std::less_equal<types::t_int> > i_cell;
    atat::rVector3d compute;

    i_cell.add(-1,1);
    i_cell.add(-1,1);
    i_cell.add(-1,1);

    do
    {
      compute(0) = (types::t_real) i_cell.access(0);
      compute(1) = (types::t_real) i_cell.access(1);
      compute(2) = (types::t_real) i_cell.access(2);

      compute = vec_a + lat*compute;
      if ( norm2( compute - vec_b ) <  ZERO_TOLERANCE ) 
        return true;

    } while ( (++i_cell) );

    return false;
  }

  // refold by one vector
  void Constituent_Strain :: refold( atat::rVector3d &vec, const atat::rMatrix3d &lat )
  {
    opt::Ndim_Iterator<types::t_int, std::less_equal<types::t_int> > i_cell;
    atat::rVector3d hold = vec;
    atat::rVector3d compute;
    atat::rVector3d current = vec;
    types::t_real norm_c = norm2(vec);

    i_cell.add(-2,2);
    i_cell.add(-2,2);
    i_cell.add(-2,2);

    do
    {
      compute(0) = (types::t_real) i_cell.access(0);
      compute(1) = (types::t_real) i_cell.access(1);
      compute(2) = (types::t_real) i_cell.access(2);

      vec = hold + lat*compute;
      if ( norm2( vec ) < norm_c ) 
      {
        current = vec;
        norm_c = norm2(vec);
      }

    } while ( (++i_cell) );

    vec = current;
  }


  // writes k_vecs into a Structure tag
  // first checks wether there already exists a Structure tag
  // if not, creates it
  // note that Atom tags should be written by some Structure object
  void Constituent_Strain :: print_xml( TiXmlElement &_node ) const
  {
    TiXmlElement *parent, *child;
    atat::rVector3d vec; 

    // checks if Structure exists
    parent = _node.FirstChildElement("Structure");
    if ( not parent )
    {
      parent = new TiXmlElement("Structure");
      _node.LinkEndChild(parent);
    }

    std::vector<atat::rVector3d> :: const_iterator i_vec = k_vecs.begin();
    std::vector<atat::rVector3d> :: const_iterator i_end = k_vecs.end();
    for( ; i_vec != i_end; ++i_vec )
    {
      child = new TiXmlElement( "kvec" );
      parent->LinkEndChild( child );
      child->SetDoubleAttribute( "x", vec(0) );
      child->SetDoubleAttribute( "y", vec(1) );
      child->SetDoubleAttribute( "z", vec(2) );
    }
  }

  // loads r_vecs and k_vecs from xml input
  bool Constituent_Strain :: Load( const TiXmlElement &_element )
  {
    const TiXmlElement *child, *str;
    Harmonic harmonic;
    atat::rVector3d vec; types::t_real d;

    // advance to Structure tag
    str = _element.FirstChildElement( "Structure" );
    if ( not str )
      return false;
    
    // read kvec tags
    child = str->FirstChildElement( "kvec" );
    if ( !child )
      return false;
    for ( ; child; child=child->NextSiblingElement( "kvec" ) )
    {
      d = 1.0; child->Attribute("x", &d); vec(0) = d;
      d = 1.0; child->Attribute("y", &d); vec(1) = d;
      d = 1.0; child->Attribute("z", &d); vec(2) = d;
      k_vecs.push_back( vec );
    }

    // read Atom tags
    child = str->FirstChildElement( "Atom" );
    if ( !child )
      return false;
    for ( ; child; child=child->NextSiblingElement( "Atom" ) )
    {
      d = 1.0; child->Attribute("x", &d); vec(0) = d;
      d = 1.0; child->Attribute("y", &d); vec(1) = d;
      d = 1.0; child->Attribute("z", &d); vec(2) = d;
      r_vecs.push_back( vec);
    }
    return true;
  }


  

  #ifdef _DEBUG_LADA_
    void Constituent_Strain :: check_derivative()
    {
      const types::t_int n_spin = 0;
      const types::t_int steps = 10000;
      types::t_int i;
      types::t_real integral;
      types::t_real *gradient = new types::t_real[variables->size()];
      const types::t_real start = 0.0, end = 0.5;
      const types::t_real inv_step = (end - start)/ ( (types::t_real) steps );

      integral = 0.0;
      for ( i = 0; i <= steps; ++i )
      {
        *(variables->begin() + n_spin ) = start + inv_step * (types::t_real) i;
        evaluate_with_gradient( gradient );
        integral += *(gradient + n_spin);
      }

      std::cout << " integral: " << integral * inv_step << std :: endl;

      *(variables->begin() + n_spin ) = end;
      integral = evaluate_with_gradient( gradient );
      *(variables->begin() + n_spin ) = start;
      integral -= evaluate_with_gradient( gradient );

      std::cout << " direct: " << integral << std :: endl;


      delete[] gradient;

    }
  #endif // _DEBUG_LADA_

} // namespace Ising_CE 




#ifdef _MPI

namespace mpi
{
  template<>
  bool BroadCast :: serialize<Ising_CE::Constituent_Strain>( Ising_CE::Constituent_Strain &_cs )
  {
    // first serializes cells
    if ( not serialize( _cs.r_cell.x[0], (_cs.r_cell.x[0]+3) ) ) return false;
    if ( not serialize( _cs.r_cell.x[1], (_cs.r_cell.x[1]+3) ) ) return false;
    if ( not serialize( _cs.r_cell.x[2], (_cs.r_cell.x[2]+3) ) ) return false;
    if ( not serialize( _cs.k_cell.x[0], (_cs.k_cell.x[0]+3) ) ) return false;
    if ( not serialize( _cs.k_cell.x[1], (_cs.k_cell.x[1]+3) ) ) return false;
    if ( not serialize( _cs.k_cell.x[2], (_cs.k_cell.x[2]+3) ) ) return false;
    
    // then serializes rvecs and kvecs
    types::t_int n = _cs.r_vecs.size();
    if( not serialize( n ) ) return false;
    if ( stage == COPYING_FROM_HERE )
      _cs.r_vecs.resize(n);
    n = _cs.k_vecs.size();
    if( not serialize( n ) ) return false;
    if ( stage == COPYING_FROM_HERE )
      _cs.k_vecs.resize(n);
    std::vector<atat::rVector3d> :: iterator i_vec = _cs.r_vecs.begin();
    std::vector<atat::rVector3d> :: iterator i_vec_end = _cs.r_vecs.end();
    for(; i_vec != i_vec_end; ++i_vec )
      if ( not serialize( i_vec->x, i_vec->x+3 ) ) return false;
    i_vec = _cs.k_vecs.begin();
    i_vec_end = _cs.k_vecs.end();
    for(; i_vec != i_vec_end; ++i_vec )
      if ( not serialize( i_vec->x, i_vec->x+3 ) ) return false;

    n = _cs.harmonics.size();
    if( not serialize( n ) ) return false;
    if ( stage == COPYING_FROM_HERE )
      _cs.harmonics.resize(n);
    Ising_CE::Constituent_Strain::t_Harmonics :: iterator i_harm = _cs.harmonics.begin();
    Ising_CE::Constituent_Strain::t_Harmonics :: iterator i_harm_end = _cs.harmonics.end();
    for(; i_harm != i_harm_end; ++i_harm )
      if ( not serialize( *i_harm ) ) return false;

    return true;
  }

}

#endif
