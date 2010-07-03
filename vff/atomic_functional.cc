#include "LaDaConfig.h"

#include <cstdlib>

#include <algorithm>
#include <functional>
#include <boost/lexical_cast.hpp>

#include <physics/physics.h>
#include <opt/debug.h>
#include <opt/tinyxml.h>

#include "atomic_functional.h"
#include "atomic_center.h"
  
namespace LaDa
{
  namespace Vff
  { 
    namespace details
    {
      //! Finds a bond node for species \a _A, \a _B.
      const TiXmlElement* find_bond_node( const TiXmlElement &_node, 
                                          const std::string &_A, 
                                          const std::string &_B );
      //! \brief Finds a angle node for species \a _A, \a _B, \a _C.
      //! \details \a _B is the center of the angle.
      const TiXmlElement* find_angle_node( const TiXmlElement &_node, 
                                           const std::string &_B, 
                                           const std::string &_A, 
                                           const std::string &_C );
    } // end of details namespace

    // constants obtained from bc -l with scale = 64
    const types::t_real // 2*sqrt(3)
      AtomicFunctional :: twos3 = 3.4641016151377545870548926830117447338856105076207612561116139588;
    const types::t_real // 1 / 16
      AtomicFunctional :: one16 = 0.0625;
    const types::t_real  // sqrt(3) / 8
      AtomicFunctional :: s3o160 = 0.2165063509461096616909307926882340458678506567262975785069758724;
    const types::t_real  // 1 / 640
      AtomicFunctional :: one640 = 0.0015625;
    const types::t_real  // 3 / 8
      AtomicFunctional :: three8 = 0.375;
    const  types::t_real // 3 * sqrt(3) / 8
      AtomicFunctional :: s33o8  = 0.6495190528383289850727923780647021376035519701788927355209276172;
    const  types::t_real // 3 * sqrt(3) / 16
      AtomicFunctional :: s33o16 = 0.3247595264191644925363961890323510688017759850894463677604638086;
    const  types::t_real // 3 / 16
      AtomicFunctional :: thre16 = 0.1875; 
    const  types::t_real // 3 / 32
      AtomicFunctional :: thre32 = 0.09375;
    const  types::t_real // 3/128 *sqrt(3) 
      AtomicFunctional :: s33128 = 0.0405949408023955615670495236290438836002219981361807959700579760;
    const  types::t_real // 3/256 *sqrt(3) 
      AtomicFunctional :: s33256 = 0.0202974704011977807835247618145219418001109990680903979850289880;
    const  types::t_real                   
      AtomicFunctional :: no1280 = 0.00703125;
    const  types::t_real                   
      AtomicFunctional :: no2560 = 0.003515625;

    types::t_real AtomicFunctional :: evaluate( const AtomicCenter &_center ) const
    {
      types :: t_real energy = 0;
      types :: t_real scale2 = structure->scale * structure->scale;

      AtomicCenter :: const_iterator i_bond_begin  = _center.begin();
      AtomicCenter :: const_iterator i_bond_end    = _center.end();
      AtomicCenter :: const_iterator i_bond ( i_bond_begin );
      for( ; i_bond != i_bond_end; ++i_bond )
      {
        // Bond stretch 
        // sets up parameters
        types::t_unsigned bond_kind = i_bond.kind();
        types::t_real bond_length = lengths[bond_kind];
        
        // adds energy only if center is site 0 
        // computes e0 for bond-angle now 
        types::t_real e0 = i_bond.norm2() * scale2 / bond_length - bond_length; 
        if ( _center.site_one() ) 
        {
          std::vector<types::t_real>::const_iterator i_alpha =   alphas.begin()
                                                               + 5 * bond_kind + 4;
          types::t_real dummy = e0 * one640 * (*i_alpha); --i_alpha;
          dummy =   e0 * ( (*i_alpha) * s3o160 + dummy ); --i_alpha;
          dummy =   e0 * ( (*i_alpha) * one16 + dummy ); --i_alpha;
          dummy =   e0 * ( (*i_alpha) / twos3 + dummy ); --i_alpha;
          energy +=   e0 * e0 * ( (*i_alpha)  + dummy ); 
        } 

        // Three body terms
        AtomicCenter :: const_iterator i_angle ( i_bond_begin );
        for (; i_angle != i_bond_end; ++i_angle )
        {
          if ( i_angle == i_bond ) continue;
          
          // sets up parameters
          types::t_unsigned end_kind = i_angle.kind();
          types::t_unsigned angle_kind = bond_kind + end_kind;
          types::t_real end_length = std::sqrt( bond_length * lengths[end_kind] );
          types::t_real gamma = gammas[angle_kind];
          
          // Bond bending
          types::t_real e1 =   i_bond.scalar_product( i_angle ) * scale2 / end_length 
                             - end_length * gamma;
          if ( i_bond - i_angle > 0 )
          {
            std::vector< types::t_real > :: const_iterator i_beta =   betas.begin()
                                                                    + 5 * angle_kind + 4;
            types::t_real dummy = e1 * one640 * (*i_beta); --i_beta;
            dummy =   e1 * ( (*i_beta) * s3o160 + dummy ); --i_beta;
            dummy =   e1 * ( (*i_beta) * one16 + dummy ); --i_beta;
            dummy =   e1 * ( (*i_beta) / twos3 + dummy ); --i_beta;
            energy +=   e1 * e1 * ( (*i_beta)  + dummy ); 
          }

          // Bond angle
          energy += e1 * e0 * sigmas[angle_kind];
        }

      } // for( ; i_bond != i_bond_end; ++i_bond )

      return energy * three8;
    }

    types::t_real AtomicFunctional
       :: evaluate_with_gradient( const AtomicCenter &_center,
                                  const math::rMatrix3d &_strain,
                                  math::rMatrix3d &_stress,
                                  const math::rMatrix3d &_K0 ) const
    {
      types :: t_real energy = 0;
      types :: t_real scale2 = structure->scale * structure->scale;

      AtomicCenter :: const_iterator i_bond_begin  = _center.begin();
      AtomicCenter :: const_iterator i_bond_end    = _center.end();
      AtomicCenter :: const_iterator i_bond ( i_bond_begin );
      for( ; i_bond != i_bond_end; ++i_bond )
      {
        // Bond stretch 
        // sets up parameters
        types::t_unsigned bond_kind = i_bond.kind();
        types::t_real bond_length = lengths[bond_kind];
        math::rVector3d d0; i_bond.vector( d0 );
        
        // adds energy only if center is site 0 
        // computes e0 for bond-angle now 
        types::t_real e0 = d0.squaredNorm() * scale2 / bond_length - bond_length;
        if ( _center.site_one() ) 
        {
          // energy 
          std::vector<types::t_real>::const_iterator i_alpha =   alphas.begin()
                                                               + 5 * bond_kind + 4;
          types::t_real dummy = e0 * one640 * (*i_alpha); --i_alpha;
          dummy = e0 * ( (*i_alpha) * s3o160 + dummy ); --i_alpha;
          dummy = e0 * ( (*i_alpha) * one16 + dummy ); --i_alpha;
          dummy = e0 * ( (*i_alpha) / twos3 + dummy ); --i_alpha;
          energy +=   e0 * e0 * ( (*i_alpha)  + dummy ); 

          // then gradient 
          i_alpha += 4; // alphas.begin() + 5 * bond_kind;
          dummy = e0 * no1280 * (*i_alpha); --i_alpha;
          dummy = e0 * ( (*i_alpha) * s33128 + dummy ); --i_alpha;
          dummy = e0 * ( (*i_alpha) * thre16 + dummy ); --i_alpha;
          dummy = e0 * ( (*i_alpha) * s33o8 + dummy ); --i_alpha;
          types::t_real e0grad =   types::t_real(2.0) * scale2 / bond_length 
                                 * e0 * ( types::t_real(1.5) * (*i_alpha) + dummy); 
          math::rVector3d hold = e0grad * ( _strain * d0 );
          _center.gradient -= hold; // with respect to atomic positions
          i_bond->gradient += hold;  

          // stress
          for( int i = 0; i < 3; ++i )
            for( int j = 0; j < 3; ++j )
              for( int k = 0; k < 3; ++k )
                _stress(i,j) += d0[i] * d0[k] * _K0(k,j) * e0grad * types::t_real(0.5);
        }

        // Three body terms
        AtomicCenter :: const_iterator i_angle ( i_bond_begin );
        for (types::t_int i=0; i_angle != i_bond_end; ++i, ++i_angle )
          if ( i_angle != i_bond )
          {
            // sets up parameters
            types::t_unsigned end_kind = i_angle.kind();
            types::t_unsigned angle_kind = bond_kind + end_kind;
            types::t_real mean_length = std::sqrt( bond_length * lengths[end_kind] );
            types::t_real gamma = gammas[angle_kind];
            types::t_real sigma = sigmas[angle_kind];
            math::rVector3d d1; i_angle.vector( d1 );
            
            // Bond bending
            types::t_real e1 =   d0.dot(d1) * scale2 / mean_length 
                               - mean_length * gamma;
            if ( i_bond - i_angle > 0 )
            {
              // energy
              std::vector< types::t_real > :: const_iterator i_beta =   betas.begin() 
                                                                      + 5 * angle_kind
                                                                      + 4;
              types::t_real dummy = e1 * one640 * (*i_beta); --i_beta;
              dummy = e1 * ( (*i_beta) * s3o160 + dummy ); --i_beta;
              dummy = e1 * ( (*i_beta) * one16 + dummy ); --i_beta;
              dummy = e1 * ( (*i_beta) / twos3 + dummy ); --i_beta;
              energy +=   e1 * e1 * ( (*i_beta)  + dummy ); 
              
              // then gradient 
              i_beta += 4; // goes back to beginning of array
              dummy = e1 * no2560 * (*i_beta); --i_beta;
              dummy = e1 * ( (*i_beta) * s33256 + dummy ); --i_beta;
              dummy = e1 * ( (*i_beta) * thre32 + dummy ); --i_beta;
              dummy = e1 * ( (*i_beta) * s33o16 + dummy ); --i_beta;
              types::t_real e1grad =   types::t_real(2.0) * scale2 / mean_length
                                     * e1 * ( *(  i_beta) * types::t_real(0.75) + dummy);
              math::rVector3d hold0 = e1grad * ( _strain * d0 );
              math::rVector3d hold1 = e1grad * ( _strain * d1 );
              _center.gradient -= ( hold0 + hold1); // with respect to atomic positions
              i_bond->gradient += hold1; 
              i_angle->gradient += hold0; 
              
              // stress
              for( int i = 0; i < 3; ++i )
                for( int j = 0; j < 3; ++j )
                  for( int k = 0; k < 3; ++k )
                    _stress(i,j) +=   (d1[i] * d0[k] + d0[i] * d1[k]) 
                                    * _K0(k,j) * e1grad * types::t_real(0.5);
            }

            // Bond angle energy
            energy += e1 * e0 * sigma;

            // Bond angle gradients
            { // position gradients
              math::rVector3d hold0 = types::t_real(1.5) * e1 * sigma / bond_length * scale2
                                      * ( _strain * d0 );
              math::rVector3d hold1 = types::t_real(0.75) * e0 * sigma / mean_length * scale2
                                      * ( _strain * d1 );
              math::rVector3d hold2 = types::t_real(0.75) * e0 * sigma / mean_length * scale2
                                      * ( _strain * d0 );
              _center.gradient -= (hold0 + hold1 + hold2);
              (*i_bond).gradient += (hold0 + hold1); //(hold0 + hold1);
              (*i_angle).gradient += hold2;
            }

            // stress
            for( int i = 0; i < 3; ++i )
              for( int j = 0; j < 3; ++j )
                for( int k = 0; k < 3; ++k )
                  _stress(i,j) +=   _K0(k,j) * types::t_real(0.375) * sigma * scale2 
                                  * (   types::t_real(2) * e1 / bond_length * d0[i] * d0[k] 
                                      +   e0 / mean_length 
                                        * (d0[i] * d1[k] + d1[i] * d0[k]) );
          } // angle loop

      } // for( ; i_bond != i_bond_end; ++i_bond )

      return energy * three8;
    }

    types::t_real AtomicFunctional :: MicroStrain( const AtomicCenter &_center, 
                                                   const Crystal::Structure &_str0 ) const 
    {
      if ( _center.size() != 4 )
      { 
        std::cerr << "Microscopic strain cannot be computed "
                  << "Because atom " << _center.get_index()
                  << " " << _center.Origin().pos
                  << " has only " << _center.size() << " bonds " << std::endl;
        return 0;
      }

      // computes the microstrain as the average strain over four paralellepipeds.
      // first computes four bond vectors, bond lengths, and angles.
      math::rVector3d bonds[4];
      types::t_real l0s[4]; 
      types::t_real angles[4][4]; 
      AtomicCenter :: const_iterator i_bond  = _center.begin();
      AtomicCenter :: const_iterator i_bond_end    = _center.end();
      for( size_t i(0); i_bond != i_bond_end; ++i_bond, ++i )
      {
        bonds[i] = structure->scale *  i_bond.vector( bonds[i] );
        const types::t_unsigned bond_kind( i_bond.kind() );
        l0s[i] = lengths[ bond_kind ];

        angles[i][i] = 0;
        AtomicCenter :: const_iterator i_angle( i_bond ); ++i_angle;
        for( size_t j(i+1); i_angle != i_bond_end; ++i_angle, ++j )
        {
          if( i == j ) {  continue; }
          angles[i][j] = gammas[ bond_kind + i_angle.kind() ];
          angles[j][i] = angles[i][j];
        }
      }

      // now performs loop over all parallelipipeds
      types::t_real strain(0);
      for( size_t i(0); i < 4; ++i )
      {
        const size_t j( i+1>3 ? i-3: i+1 );
        const size_t k( j+1>3 ? j-3: j+1 );
        const types::t_real vol( std::abs( bonds[i].dot(bonds[j] ^ bonds[k]) ));
        const types::t_real vol0
        ( 
            l0s[i] * l0s[j] * l0s[k]
          * std::sqrt
            (
              1e0 - angles[i][j] * angles[i][j] - angles[j][k] * angles[j][k]
                  - angles[k][i] * angles[k][i]
                  + 2e0 * angles[i][j] * angles[j][k] * angles[k][i] 
            )
         );
        strain += ( vol / vol0 - 1e0 );
      }
      strain *= 0.25e0;
      return strain;
    }

    void AtomicFunctional :: print_out( std::ostream &stream ) const
    {
      stream << "Site " << specie << " " << site << std::endl << "  ";
      std::vector<types::t_real> :: const_iterator i_var = lengths.begin();
      std::vector<types::t_real> :: const_iterator i_end = lengths.end();
      for(; i_var != i_end; ++i_var )
        stream << *i_var << "  ";
      stream << std::endl << "  ";
      i_var = alphas.begin();
      i_end = alphas.end();
      for(; i_var != i_end; ++i_var )
        stream << *i_var << "  ";
      stream << std::endl << "  ";
      i_var = betas.begin();
      i_end = betas.end();
      for(; i_var != i_end; ++i_var )
        stream << *i_var << "  ";
      stream << std::endl << "  ";
      i_var = gammas.begin();
      i_end = gammas.end();
      for(; i_var != i_end; ++i_var )
        stream << *i_var << "  ";
      stream << std::endl << "  ";
      i_var = sigmas.begin();
      i_end = sigmas.end();
      for(; i_var != i_end; ++i_var )
        stream << *i_var << "  ";
      stream << std::endl;
    }

    boost::tuples::tuple< const types::t_real&, const types::t_real&, const types::t_real&,
                          const types::t_real&, const types::t_real&, const types::t_real& >
      AtomicFunctional::get_bond( size_t _kind ) const
      {
        namespace bt = boost::tuples;
        __ASSERT( _kind >= lengths.size(), "Index out-of-range.\n" )
        __ASSERT( 5*_kind + 4 >= alphas.size(), "Index out-of-range.\n" )
        return boost::tuples::tie( lengths[_kind], alphas[5*_kind], alphas[5*_kind+1], 
                                   alphas[5*_kind+2], alphas[5*_kind+3], alphas[5*_kind+4] );
      }

    boost::tuples::tuple< const types::t_real&, const types::t_real&, 
                          const types::t_real&, const types::t_real&,
                          const types::t_real&, const types::t_real&, const types::t_real& >
      AtomicFunctional::get_angle( size_t _kind ) const
      {
          __ASSERT( _kind >= gammas.size(), "Index out-of-range.\n" )
          __ASSERT( _kind >= sigmas.size(), "Index out-of-range.\n" )
          __ASSERT( 5*_kind + 4 >= betas.size(), "Index out-of-range.\n" )
          return boost::tuples::tie( gammas[_kind], sigmas[_kind],
                                     betas[5*_kind], betas[5*_kind+1], 
                                     betas[5*_kind+2], betas[5*_kind+3], betas[5*_kind+4] );
      }

    bool AtomicFunctional :: load( const TiXmlElement& _node,
                                    const types::t_unsigned &_site_index,
                                    const types::t_unsigned &_type_index,
                                    const Crystal::Lattice::t_Site &_othersite,
                                    Crystal :: Structure &_structure )
    {
      __ASSERT( std::string( _node.Value() ).compare("Functional") != 0, "Incorrect node.\n" )
      __ASSERT( std::string( _node.Attribute( "type" ) ).compare("vff") != 0,
                "Incorrect node.\n" )

      const Crystal::Lattice &lattice = *Crystal::Structure::lattice;
      __ASSERT( _site_index >= lattice.sites.size(), "Index out of range.\n" )
      const Crystal::Lattice::t_Site &lsite = lattice.sites[ _site_index ];
      __ASSERT( _type_index >= lsite.type.size(), "Index out of range.\n" )

      structure = &_structure;
      specie = lsite.type[ _type_index ];
      site = _site_index;
      type = _type_index;
      const types::t_unsigned nbonds( _othersite.type.size() );
      const types::t_unsigned nangles( 2*nbonds - 1 );

      // Resize parameter arrays.
      lengths.clear(); lengths.resize( nbonds, 0e0 );     
      alphas.clear(); alphas.resize( 5 * nbonds, 0e0 );   
      sigmas.clear(); sigmas.resize( nangles, 0e0 );     
      gammas.clear(); gammas.resize( nangles, 0e0 ); 
      betas.clear(); betas.resize( 5 * nangles, 0e0 );  

      // Now reads parameters.
      size_t bond_index(0), angle_index(0);
      Crystal::Lattice::t_Site::t_Type::const_iterator i_atype = _othersite.type.begin();
      Crystal::Lattice::t_Site::t_Type::const_iterator i_atype_end = _othersite.type.end();
      for(; i_atype != i_atype_end; ++i_atype, ++bond_index )
      {
        // first the bond parameters.
        const TiXmlElement *child = details::find_bond_node( _node, specie, *i_atype );
        __DOASSERT( not child,
                       "Could not find parameters for bond: "
                    << specie << "-" << (*i_atype) << ".\n" )
        __DOASSERT( not child->Attribute("d0"), 
                       "Could not find parameters d0 of bond: "
                    << specie << "-" << (*i_atype) << ".\n" )
        __DOASSERT( not child->Attribute("alpha"), 
                       "Could not find parameters alpha of bond: "
                    << specie << "-" << (*i_atype) << ".\n" )

        lengths[ bond_index ] = boost::lexical_cast<types::t_real>( child->Attribute( "d0" ) );
        alphas[ 5 * bond_index ]
          = boost::lexical_cast<types::t_real>( child->Attribute( "alpha" ) );
        for( size_t i(3); i < 7; ++i )
        {
          const std::string str(   std::string("alpha") 
                                 + boost::lexical_cast<std::string>(i) );

          if( not child->Attribute(str) ) continue;

          alphas[ 5 * bond_index + i - 2 ]
            = boost::lexical_cast<types::t_real>( child->Attribute(str.c_str()) );
        }
        
        // then the angle and bond-angle parameters.
        Crystal::Lattice::t_Site::t_Type::const_iterator i_btype = i_atype;
        for(; i_btype != i_atype_end; ++i_btype, ++angle_index )
        {
          child = details::find_angle_node( _node, specie, *i_atype, *i_btype );
          
          __DOASSERT( not child,
                         "Could not find parameters for angle: "
                      << (*i_atype) << "-" << specie << "-" << (*i_btype) << ".\n" )
          __DOASSERT( not child->Attribute("gamma"), 
                         "Could not find parameters gamma of angle: "
                      << (*i_atype) << "-" << specie << "-" << (*i_btype) << ".\n" )
          __DOASSERT( not child->Attribute("beta"), 
                         "Could not find parameters beta of angle: "
                      << (*i_atype) << "-" << specie << "-" << (*i_btype) << ".\n" )

          // find gamma
          std::string str( child->Attribute("gamma") );
          if ( str.compare("tet") == 0 or str.compare("tetrahedral") == 0  )
            gammas[angle_index] = -1e0/3e0;
          else gammas[angle_index] = boost::lexical_cast<types::t_real>( str );
          __DOASSERT( std::abs(gammas[angle_index]) > 1,
                      " gamma must be comprised between 1 and -1.\n" )

          // find sigma
          if( child->Attribute( "sigma" ) )
            sigmas[angle_index]
              = boost::lexical_cast<types::t_real>( child->Attribute( "sigma" ) );
          // find betas
          betas[ 5 * angle_index ] 
             = boost::lexical_cast<types::t_real>( child->Attribute( "beta" ) );
          for( size_t i(3); i < 7; ++i )
          {
            const std::string str(   std::string("beta") 
                                   + boost::lexical_cast<std::string>(i) );
          
            if( not child->Attribute(str.c_str()) ) continue;
          
            betas[ 5 * angle_index + i - 2 ]
              = boost::lexical_cast< types::t_real>( child->Attribute(str.c_str()) );
          } // end of loop over betas.
        } // end of loop over angle kinds.
      } // end of loop over bond kinds.

      return true;
    }

    namespace details
    {
      const TiXmlElement* find_bond_node( const TiXmlElement &_node, 
                                          const std::string &_A, 
                                          const std::string &_B )
      {
        const TiXmlElement *result = NULL;
        const TiXmlElement *child( _node.FirstChildElement("Bond") );
        size_t partial_match(0);
        size_t exact_match(0);
        for(; child; child = child->NextSiblingElement( "Bond" ) )
        {
          if( not ( child->Attribute( "A" ) and child->Attribute( "B" ) ) ) 
          {
            std::cerr << "Found incomplete Bond tag in xml input.\nTag is ignored.\n";
            continue;
          }
          const std::string A = child->Attribute("A");
          const std::string B = child->Attribute("B");
          // exact match.
          if( A.compare( _A ) == 0 and B.compare( _B ) == 0 )
            { ++exact_match; result = child; }
          // partial match.
          if( B.compare( _A ) == 0 and A.compare( _B ) == 0 ) 
            { ++partial_match; if( not exact_match ) result = child; }
        }

        if( exact_match == 1 ) return result;
        if( exact_match > 1 ) 
        {
          std::cerr << "Found more than one match for Bond tag " 
                    << _A << "-" << _B <<".\nAborting.\n";
          return NULL;
        }

        if( partial_match < 2 ) return result;

        std::cerr << "Found more than one acceptable Bond tag for "
                  << _A << "-" << _B
                  << ".\n Aborting\n.";

        return NULL;
      }

      const TiXmlElement* find_angle_node( const TiXmlElement &_node, 
                                           const std::string &_B, 
                                           const std::string &_A, 
                                           const std::string &_C )
      {
        const TiXmlElement *result = NULL;
        const TiXmlElement *child( _node.FirstChildElement("Angle") );
        size_t match(0);
        for(; child; child = child->NextSiblingElement( "Angle" ) )
        {
          if( not (      child->Attribute( "A" )
                     and child->Attribute( "B" ) 
                     and child->Attribute( "C" ) ) ) 
          {
            std::cerr << "Found incomplete Angle tag in xml input.\nTag is ignored.\n";
            continue;
          }
          const std::string B = child->Attribute("B");
          // no match
          if( B.compare( _B ) != 0 ) continue;

          const std::string A = child->Attribute("A");
          const std::string C = child->Attribute("C");
          // match.
          if(     ( A.compare( _A ) != 0 or C.compare( _C ) != 0 )
              and ( C.compare( _A ) != 0 or A.compare( _C ) != 0 )  ) continue;
          ++match;
          result = child;
        }

        if( match < 2 ) return result;

        std::cerr << "Found more than one acceptable Angle tag for "
                  << _A << "-" << _B << "-" << _C
                  << ".\n Aborting\n.";

        return NULL;
      } 
    } // end of details namespace
  } // namespace vff
} // namespace LaDa
