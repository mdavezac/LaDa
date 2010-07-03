#include<numeric>
#include<complex>
#include<iomanip>
#include<math.h>
#include<limits.h>

#include <boost/lambda/lambda.hpp>

#include <boost/lambda/lambda.hpp>
#ifdef _MPI
#include <boost/mpi/collectives.hpp>
#endif

namespace LaDa
{
  namespace CE
  {
    namespace ConstituentStrain
    {
      template<class T_HARMONIC>
      typename Functional<T_HARMONIC>::t_Harmonics Functional<T_HARMONIC> :: harmonics;

      template<class T_HARMONIC>
      void Functional<T_HARMONIC> :: operator<<( const Crystal::Structure &_str )
      {
        __DOASSERT( _str.atoms.size() < 1,  "No atoms in structure.\n" 
                                            "Cannot create constituent strain"
                                            " for structure.\n" )
        __DOASSERT( _str.k_vecs.size() < 1, "No kvectors in structure.\n" 
                                            "Cannot create constituent strain"
                                            " for structure.\n" )
        k_vecs.clear(); r_vecs.clear();
        k_vecs.reserve(_str.k_vecs.size()); r_vecs.reserve(_str.atoms.size());
        Crystal::Structure::t_kAtoms::const_iterator i_kvec = _str.k_vecs.begin();
        Crystal::Structure::t_kAtoms::const_iterator i_kvec_end = _str.k_vecs.end();
        for(; i_kvec != i_kvec_end; ++i_kvec )
          k_vecs.push_back( i_kvec->pos );
        Crystal::Structure::t_Atoms::const_iterator i_rvec = _str.atoms.begin();
        Crystal::Structure::t_Atoms::const_iterator i_rvec_end = _str.atoms.end();
        for(; i_rvec != i_rvec_end; ++i_rvec )
          r_vecs.push_back( i_rvec->pos );
      }

      template<class T_HARMONIC>
      bool Functional<T_HARMONIC> :: Load_Harmonics( const TiXmlElement &_element)
      {
        try
        {
          const TiXmlElement *child;
          t_Harmonic harmonic;
          types::t_real d=1.0;
          
          if ( _element.Attribute("attenuation", &d) )
            t_Harmonic::set_attenuation( d );
          
          child = _element.FirstChildElement( "Harmonic" );
          __DOASSERT( not child, "Could not find <Harmonic> tag " 
                                 "when loading constituent strain" )
          harmonics.clear();
          for ( ; child; child=child->NextSiblingElement( "Harmonic" ) )
          {
            if( not harmonic.Load( *child) ) continue;
            harmonics.push_back( harmonic );
            harmonic.clear();
          }
        }
        __CATCHCODE(, "Error while parsing <Harmonic> tag for constituent strain.\n" )
        __DOASSERT( not harmonics.size(),
                    "Could find any harmonics in input.\n" )

        return true;
      }


      template<class T_HARMONIC>
      types::t_real Functional<T_HARMONIC> :: evaluate()
      {
        __ASSERT( variables->size() == 0, 
                  "variables have not been initialized.\n" )
        types::t_real sum_harm;
        std::complex<types::t_real> sum_exp;
        const std::complex<types::t_real> imath(0, -2*3.1415926535897932384626433832795028841971693993751058208);
        types::t_real x = get_concentration();

        
        std::vector<math::rVector3d> :: const_iterator i_k_vec = k_vecs.begin();
        std::vector<math::rVector3d> :: const_iterator i_k_vec_end __SERIALCODE( = k_vecs.end() );
        std::vector<math::rVector3d> :: const_iterator i_r_vec;
        std::vector<math::rVector3d> :: const_iterator i_r_vec_begin = r_vecs.begin();
        std::vector<math::rVector3d> :: const_iterator i_r_vec_end = r_vecs.end();
        typename t_Harmonics :: const_iterator i_harmonic;
        typename t_Harmonics :: const_iterator i_harmonic_begin = harmonics.begin();
        typename t_Harmonics :: const_iterator i_harmonic_end = harmonics.end();
        std::vector<types::t_real> :: const_iterator i_spin;
        std::vector<types::t_real> :: const_iterator i_spin_begin = variables->begin();
        types::t_real *interpolation = new types::t_real[ harmonics.size() ];
        types::t_real *i_inter;
        if (! interpolation )
        {
          std::cerr << "Memory allocation error in "
                    << "Functional :: evaluate(...) "
                    << std::endl;
          exit(0);
        }
        i_harmonic = i_harmonic_begin;
        i_inter = interpolation;
        for ( ; i_harmonic != i_harmonic_end; ++i_harmonic, ++i_inter )
          *i_inter = i_harmonic->evaluate(x);

        types::t_real value = 0.0;
        __MPICODE(
          __ASSERT( not comm, "Communicator not set.\n" )
          types :: t_unsigned nperproc = k_vecs.size() / comm->size(); 
          types :: t_unsigned remainder = k_vecs.size() % comm->size();
          i_k_vec +=  comm->rank() * nperproc + std::min( (types::t_int) remainder,
                                                          comm->rank() );
          i_k_vec_end = i_k_vec + nperproc;
          if( remainder and comm->rank() < remainder ) ++i_k_vec_end;
        )
        for ( ; i_k_vec != i_k_vec_end; ++i_k_vec )
        {
          if ( i_k_vec->squaredNorm() < types::tolerance ) // don't need to compute Gamma
            continue;

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
        __MPICODE( value = boost::mpi::all_reduce( *comm, value,
                                                   std::plus<types::t_real>() ); )
        value /= ( (types::t_real)  k_vecs.size() ) * (types::t_real) r_vecs.size();

        delete[] interpolation;
        return value;
      }

      template<class T_HARMONIC>
      types::t_real Functional<T_HARMONIC> :: evaluate_one_gradient( types::t_unsigned _pos ) 
      {
        __ASSERT( variables->size() == 0, 
                  "variables have not been initialized.\n" )
        std::vector<math::rVector3d> :: const_iterator i_k_vec = k_vecs.begin();
        std::vector<math::rVector3d> :: const_iterator i_k_vec_end __SERIALCODE( = k_vecs.end() );
        std::vector<math::rVector3d> :: const_iterator i_r_vec;
        std::vector<math::rVector3d> :: const_iterator i_r_vec_begin = r_vecs.begin();
        std::vector<math::rVector3d> :: const_iterator i_r_vec_end = r_vecs.end();
        typename t_Harmonics :: const_iterator i_harmonic;
        typename t_Harmonics :: const_iterator i_harmonic_begin = harmonics.begin();
        typename t_Harmonics :: const_iterator i_harmonic_end = harmonics.end();
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
                    << "Functional :: evaluate(...) "
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

        __MPICODE(
          types :: t_unsigned nperproc = k_vecs.size() / comm->size(); 
          types :: t_unsigned remainder = k_vecs.size() % comm->size();
          i_k_vec +=  comm->rank() * nperproc + std::min( (types::t_int) remainder,
                                                          comm->rank() );
          i_k_vec_end = i_k_vec + nperproc;
          if( remainder and comm->rank() < remainder ) ++i_k_vec_end;
        )
        for ( ; i_k_vec != i_k_vec_end; ++i_k_vec )
        {
          if ( i_k_vec->squaredNorm() < types::tolerance ) // don't need to compute Gamma
            continue;

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

        __MPICODE( result = boost::mpi::all_reduce( *comm, result, 
                                                    std::plus<types::t_real>() ); )
        delete[] interpolation;
        return result * inv_N;
      }
      template<class T_HARMONIC>
      types::t_real Functional<T_HARMONIC>
        :: evaluate_with_gradient( types::t_real* const gradient ) 
      {
        __ASSERT( variables->size() == 0, 
                  "variables have not been initialized.\n" )
        std::vector<math::rVector3d> :: const_iterator i_k_vec = k_vecs.begin();
        std::vector<math::rVector3d> :: const_iterator i_k_vec_end __SERIALCODE( = k_vecs.end() ); 
        std::vector<math::rVector3d> :: const_iterator i_r_vec;
        std::vector<math::rVector3d> :: const_iterator i_r_vec_begin = r_vecs.begin();
        std::vector<math::rVector3d> :: const_iterator i_r_vec_end = r_vecs.end();
        typename t_Harmonics :: const_iterator i_harmonic;
        typename t_Harmonics :: const_iterator i_harmonic_begin = harmonics.begin();
        typename t_Harmonics :: const_iterator i_harmonic_end = harmonics.end();
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
                    << "Functional::evaluate_with_gradient(..)"
                    << std::endl;
          exit(0);
        }

        // computes interpolation stuff and derivatives
        types::t_real *interpolation = new types::t_real[ 2*harmonics.size() ];
        types::t_real *i_inter, *ig_inter;
        if (! interpolation )
        {
          std::cerr << "Memory allocation error in "
                    << "Functional :: evaluate(...) "
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

        std::fill( gradient, gradient + variables->size(), t_Type(0) );

        types::t_real value = 0.0;
  #  ifdef __ADDHERE__
  #  error Please change __ADDHERE__ to something not already defined
  #  endif
  #  define __ADDHERE__ __MPISERIALCODE( array, gradient )
        __MPICODE(
          types :: t_unsigned nperproc = k_vecs.size() / comm->size(); 
          types :: t_unsigned remainder = k_vecs.size() % comm->size();
          i_k_vec +=  comm->rank() * nperproc + std::min( (types::t_int) remainder,
                                                          comm->rank() );
          i_k_vec_end = i_k_vec + nperproc;
          if( remainder and comm->rank() < remainder ) ++i_k_vec_end;
          types::t_real *__ADDHERE__ = new types::t_real[ variables->size() ];
          std::fill( __ADDHERE__, __ADDHERE__ + variables->size(), types::t_real(0) );
        )
        for ( ; i_k_vec != i_k_vec_end; ++i_k_vec )
        {
          if ( i_k_vec->squaredNorm() < types::tolerance ) // don't need to compute Gamma
            continue;

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

          grad = __ADDHERE__;
          i_exp = exp_begin;
          for ( ; i_exp != i_exp_end; ++grad, ++i_exp)
            *grad +=   real( sum_exp  * conj( sum_exp ) ) * deriv_harm 
                     + real( (*i_exp) * conj( sum_exp ) ) * sum_harm * 2.0;

          value += real(sum_exp * conj( sum_exp )) * sum_harm;
        }
        __MPICODE( 
          boost::mpi::all_reduce( *comm, grad, k_vecs.size(), grad,
                                  std::plus<types::t_real>() ); 
          std::transform( gradient, gradient + variables->size(), __ADDHERE__, gradient,
                          boost::lambda::_1 + boost::lambda::_2 * boost::lambda::constant(inv_N) ); 
          delete[] __ADDHERE__;
          value = boost::mpi::all_reduce( *comm, value,
                                          std::plus<types::t_real>() );
        )
        __SERIALCODE( 
          std::for_each( gradient, gradient + variables->size(),
                         boost::lambda::_1 *= boost::lambda::constant( inv_N ) );
        )
        value *= inv_N;

        delete[] exp_begin;
        delete[] interpolation;

        return value;
  #undef __ADDHERE__
      }

      // writes k_vecs into a Structure tag
      // first checks wether there already exists a Structure tag
      // if not, creates it
      // note that Atom tags should be written by some Structure object
      template<class T_HARMONIC>
      void Functional<T_HARMONIC> :: print_xml( TiXmlElement &_node ) const
      {
        TiXmlElement *parent, *child;
        math::rVector3d vec; 

        // checks if Structure exists
        parent = _node.FirstChildElement("Structure");
        if ( not parent )
        {
          parent = new TiXmlElement("Structure");
          _node.LinkEndChild(parent);
        }

        std::vector<math::rVector3d> :: const_iterator i_vec = k_vecs.begin();
        std::vector<math::rVector3d> :: const_iterator i_end = k_vecs.end();
        for( ; i_vec != i_end; ++i_vec )
        {
          child = new TiXmlElement( "kvec" );
          parent->LinkEndChild( child );
          child->SetDoubleAttribute( "x", (*i_vec)(0) );
          child->SetDoubleAttribute( "y", (*i_vec)(1) );
          child->SetDoubleAttribute( "z", (*i_vec)(2) );
        }
      }

      // loads r_vecs and k_vecs from xml input
      template<class T_HARMONIC>
      bool Functional<T_HARMONIC> :: Load( const TiXmlElement &_element )
      {
        const TiXmlElement *child, *str;
        t_Harmonic harmonic;
        math::rVector3d vec; types::t_real d;

        // advance to Structure tag
        str = _element.FirstChildElement( "Structure" );
        __DOASSERT( not str,    "Could not find <Structure> tag when loading"
                             << "constituent strain.\n" )
        
        // read kvec tags
        child = str->FirstChildElement( "kvec" );
        __DOASSERT( not child,    "Could not find <kvec> tag when loading"
                               << "constituent strain.\n" )
        for ( ; child; child=child->NextSiblingElement( "kvec" ) )
        {
          d = 1.0; child->Attribute("x", &d); vec(0) = d;
          d = 1.0; child->Attribute("y", &d); vec(1) = d;
          d = 1.0; child->Attribute("z", &d); vec(2) = d;
          k_vecs.push_back( vec );
        }

        // read Atom tags
        child = str->FirstChildElement( "Atom" );
        __DOASSERT( not child,    "Could not find <Atom> tag when loading"
                               << "constituent strain.\n" )
        for ( ; child; child=child->NextSiblingElement( "Atom" ) )
        {
          d = 1.0; child->Attribute("x", &d); vec(0) = d;
          d = 1.0; child->Attribute("y", &d); vec(1) = d;
          d = 1.0; child->Attribute("z", &d); vec(2) = d;
          r_vecs.push_back( vec);
        }
        return true;
      }


      

#     ifdef _LADADEBUG
        template<class T_HARMONIC>
        void Functional<T_HARMONIC> :: check_derivative()
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
#     endif // _LADADEBUG
    }
  } // namespace CE
} // namespace LaDa
