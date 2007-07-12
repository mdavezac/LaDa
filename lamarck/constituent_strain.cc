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


#include "analysis/analyze_code.h"

#include "constituent_strain.h"

namespace Ising_CE 
{
  std::vector<Harmonic> Constituent_Strain :: harmonics;
  const types::t_real Constituent_Strain::ZERO_TOLERANCE = 0.0000000001;

  Constituent_Strain :: Constituent_Strain   ( const Ising_CE::Structure& _str, t_Container *vars )
                                           : function::Base<t_Type>(vars)
  {
    if (    _str.atoms.size() < 1 
         or _str.k_vecs.size() < 1 )
      return;
    k_vecs.clear(); r_vecs.clear();
    k_vecs.reserve(_str.k_vecs.size()); r_vecs.reserve(_str.atoms.size());
    Ising_CE::Structure::t_kAtoms::const_iterator i_kvec = _str.k_vecs.begin();
    Ising_CE::Structure::t_kAtoms::const_iterator i_kvec_end = _str.k_vecs.end();
    for(; i_kvec != i_kvec_end; ++i_kvec )
      k_vecs.push_back( i_kvec->pos );
    Ising_CE::Structure::t_Atoms::const_iterator i_rvec = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms::const_iterator i_rvec_end = _str.atoms.end();
    for(; i_rvec != i_rvec_end; ++i_rvec )
      r_vecs.push_back( i_rvec->pos );
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
