//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <cstdlib>

#include <algorithm>
#include <functional>
#include <boost/filesystem/operations.hpp>

#include <physics/physics.h>
#include <opt/ndim_iterator.h>
#include <opt/atat.h>
#include <opt/debug.h>
#include <opt/atat.h>
#include <opt/tinyxml.h>

#include "atomic_functional.h"
#include "atomic_center.h"
  
namespace Vff
{ 
  namespace details
  {
    //! Finds a bond node for species \a _A, \a _B.
    const TiXmlElement* find_bond_node( const TiXmlElement &_node, 
                                        const std::string &_A, 
                                        const std::string &_B );
    //! \brief Finds a angle node for species \a _A, \a _B, \a _C.
    //! \details \a _C is the center of the angle.
    const TiXmlElement* find_bond_node( const TiXmlElement &_node, 
                                        const std::string &_C, 
                                        const std::string &_A, 
                                        const std::string &_B );
  } // end of details namespace

  // constants obtained from bc -l with scale = 64
  const types::t_real // 2*sqrt(3)
    Atomic_Functional :: twos3 = 3.4641016151377545870548926830117447338856105076207612561116139588;
  const types::t_real // 1 / 16
    Atomic_Functional :: one16 = 0.0625;
  const types::t_real  // sqrt(3) / 8
    Atomic_Functional :: s3o160 = 0.2165063509461096616909307926882340458678506567262975785069758724;
  const types::t_real  // 1 / 640
    Atomic_Functional :: one640 = 0.0015625;
  const types::t_real  // 3 / 8
    Atomic_Functional :: three8 = 0.375;
  const  types::t_real // 3 * sqrt(3) / 8
    Atomic_Functional :: s33o8  = 0.6495190528383289850727923780647021376035519701788927355209276172;
  const  types::t_real // 3 * sqrt(3) / 16
    Atomic_Functional :: s33o16 = 0.3247595264191644925363961890323510688017759850894463677604638086;
  const  types::t_real // 3 / 16
    Atomic_Functional :: thre16 = 0.1875; 
  const  types::t_real // 3 / 32
    Atomic_Functional :: thre32 = 0.09375;
  const  types::t_real // 3/128 *sqrt(3) 
    Atomic_Functional :: s33128 = 0.0405949408023955615670495236290438836002219981361807959700579760;
  const  types::t_real // 3/256 *sqrt(3) 
    Atomic_Functional :: s33256 = 0.0202974704011977807835247618145219418001109990680903979850289880;
  const  types::t_real                   
    Atomic_Functional :: no1280 = 0.00703125;
  const  types::t_real                   
    Atomic_Functional :: no2560 = 0.003515625;

  types::t_real Atomic_Functional :: evaluate( const Atomic_Center &_center ) const
  {
    types :: t_real energy = 0;
    types :: t_real scale2 = structure->scale * structure->scale;

    Atomic_Center :: const_iterator i_bond_begin  = _center.begin();
    Atomic_Center :: const_iterator i_bond_end    = _center.end();
    Atomic_Center :: const_iterator i_bond ( i_bond_begin );
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
      Atomic_Center :: const_iterator i_angle ( i_bond_begin );
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

  types::t_real Atomic_Functional
     :: evaluate_with_gradient( Atomic_Center &_center,
                                const atat::rMatrix3d &_strain,
                                atat::rMatrix3d &_stress,
                                const atat::rMatrix3d &_K0 ) const
  {
    types :: t_real energy = 0;
    types :: t_real scale2 = structure->scale * structure->scale;

    Atomic_Center :: const_iterator i_bond_begin  = _center.begin();
    Atomic_Center :: const_iterator i_bond_end    = _center.end();
    Atomic_Center :: const_iterator i_bond ( i_bond_begin );
    for( ; i_bond != i_bond_end; ++i_bond )
    {
      // Bond stretch 
      // sets up parameters
      types::t_unsigned bond_kind = i_bond.kind();
      types::t_real bond_length = lengths[bond_kind];
      atat::rVector3d d0; i_bond.vector( d0 );
      
      // adds energy only if center is site 0 
      // computes e0 for bond-angle now 
      types::t_real e0 = norm2(d0) * scale2 / bond_length - bond_length;
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
        types::t_real e0grad =   2.0 * scale2 / bond_length 
                               * e0 * ( 1.5e0 * (*i_alpha) + dummy); 
        atat::rVector3d hold = e0grad * _strain * d0;
        _center.get_gradient() -= hold; // with respect to atomic positions
        i_bond->get_gradient() += hold;  

        // stress
        for( int i = 0; i < 3; ++i )
          for( int j = 0; j < 3; ++j )
            for( int k = 0; k < 3; ++k )
              _stress(i,j) += d0[i] * d0[k] * _K0(k,j) * e0grad * 0.5;
      }

      // Three body terms
      Atomic_Center :: const_iterator i_angle ( i_bond_begin );
      for (types::t_int i=0; i_angle != i_bond_end; ++i, ++i_angle )
        if ( i_angle != i_bond )
        {
          // sets up parameters
          types::t_unsigned end_kind = i_angle.kind();
          types::t_unsigned angle_kind = bond_kind + end_kind;
          types::t_real mean_length = std::sqrt( bond_length * lengths[end_kind] );
          types::t_real gamma = gammas[angle_kind];
          types::t_real sigma = sigmas[angle_kind];
          atat::rVector3d d1; i_angle.vector( d1 );
          
          // Bond bending
          types::t_real e1 =   (d0*d1) * scale2 / mean_length 
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
            types::t_real e1grad =   2.0 * scale2 / mean_length
                                   * e1 * ( *(  i_beta) * 0.75e0 + dummy);
            atat::rVector3d hold0 = e1grad * _strain * d0;
            atat::rVector3d hold1 = e1grad * _strain * d1;
            _center.get_gradient() -= ( hold0 + hold1); // with respect to atomic positions
            i_bond->get_gradient() += hold1; 
            i_angle->get_gradient() += hold0; 
            
            // stress
            for( int i = 0; i < 3; ++i )
              for( int j = 0; j < 3; ++j )
                for( int k = 0; k < 3; ++k )
                  _stress(i,j) +=   (d1[i] * d0[k] + d0[i] * d1[k]) 
                                  * _K0(k,j) * e1grad * 0.5;
          }

          // Bond angle energy
          energy += e1 * e0 * sigma;

          // Bond angle gradients
          { // position gradients
            atat::rVector3d hold0 = 1.5 * e1 * sigma / bond_length * scale2
                                    * ( _strain * d0 );
            atat::rVector3d hold1 = 0.75 * e0 * sigma / mean_length * scale2
                                    * ( _strain * d1 );
            atat::rVector3d hold2 = 0.75 * e0 * sigma / mean_length * scale2
                                    * ( _strain * d0 );
            _center.get_gradient() -= (hold0 + hold1 + hold2);
            (*i_bond).get_gradient() += (hold0 + hold1); //(hold0 + hold1);
            (*i_angle).get_gradient() += hold2;
          }

          // stress
          for( int i = 0; i < 3; ++i )
            for( int j = 0; j < 3; ++j )
              for( int k = 0; k < 3; ++k )
                _stress(i,j) +=   _K0(k,j) * 0.375 * sigma * scale2 
                                * (   2.0 * e1 / bond_length * d0[i] * d0[k] 
                                    +   e0 / mean_length 
                                      * (d0[i] * d1[k] + d1[i] * d0[k]) );
        } // angle loop

    } // for( ; i_bond != i_bond_end; ++i_bond )

    return energy * three8;
  }

  types::t_real Atomic_Functional
    :: MicroStrain( const Atomic_Center &_center, 
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
    // constructs the two matrices
    // tetra is the original tetrahedron (from _str0), 
    // dtetra the deformed tetrahedron (from _center )
    atat::rMatrix3d tetra0, dtetra;
    atat::rVector3d R0, R1, dR0, dR1;

    Atomic_Center :: const_iterator i_bond  = _center.begin();

    // first vector
    i_bond.vector(dR0); 
    R0 =   _str0.atoms[i_bond->get_index()].pos
         - _str0.atoms[ _center.get_index() ].pos;
    i_bond.translate( R0, _str0.cell );
    types::t_real aeq = lengths[i_bond.kind()]; // equilibrium lattice constant
    types::t_real deq = std::sqrt(i_bond.norm2());
    types::t_real d0eq = std::sqrt(atat::norm2(R0));
    ++i_bond;
    for( types::t_unsigned i=0; i<3; ++i_bond, ++i )
    {
      aeq += lengths[i_bond.kind()];
      i_bond.vector( dR1 );
      R1 =   _str0.atoms[i_bond->get_index()].pos
           - _str0.atoms[ _center.get_index() ].pos;
      i_bond.translate( R1, _str0.cell ); 
      R0 -= R1; dR0 -= dR1;
      tetra0.set_row( i, R0 );
      dtetra.set_row( i, dR0 );
      R0 = R1; dR0 = dR1;
      deq += std::sqrt(i_bond.norm2());
      d0eq += std::sqrt(atat::norm2(R0));
    }

    return atat::trace( dtetra * (!tetra0) ) / aeq * d0eq * _str0.scale - 3.0;
  }

  void Atomic_Functional :: print_out( std::ostream &stream ) const
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

  void Atomic_Functional :: add_bond( const types::t_unsigned _typeB,
                                      const types::t_real _l,
                                      const std::vector<types::t_real> &_i )
  {
    if ( lengths.size() < _typeB + 1)
      lengths.resize( _typeB+1, types::t_real(0) );
    if ( alphas.size() < (_typeB+1) * 5  )
      alphas.resize( (_typeB+1)*5, types::t_real(0) );
    lengths[_typeB] = _l;
    std::copy( _i.begin(), _i.end(), alphas.begin() + _typeB * 5 );
  }
  void Atomic_Functional :: add_bond( const types::t_unsigned _typeB,
                                      const types::t_real _l,
                                      const types::t_real _i[5] )
  {
    if ( lengths.size() < _typeB + 1)
      lengths.resize( _typeB+1, types::t_real(0) );
    if ( alphas.size() < (_typeB+1) * 5  )
      alphas.resize( (_typeB+1)*5, types::t_real(0) );
    lengths[_typeB] = _l;
    const types::t_real *i_alpha = _i;
    std::copy( i_alpha, i_alpha+5, alphas.begin() + _typeB * 5 );
  }
  void Atomic_Functional :: add_angle( const types::t_unsigned _typeA,
                                       const types::t_unsigned _typeC,
                                       const types::t_real _gamma,
                                       const types::t_real _sigma, 
                                       const std::vector<types::t_real> &_i )
  {
    types::t_unsigned offset = _typeA+_typeC;
    if ( gammas.size() < offset + 1  )
      gammas.resize( offset + 1, types::t_real(0) );
    if ( sigmas.size() < offset + 1  )
      sigmas.resize( offset + 1, types::t_real(0) );
    if ( betas.size() < (offset + 1) * 5 )
      betas.resize( (offset+1)*5, types::t_real(0) );
    gammas[offset] = _gamma;
    sigmas[offset] = _sigma;
    std::copy( _i.begin(), _i.end(), betas.begin() + offset * 5 );
  }

  void Atomic_Functional :: add_angle( const types::t_unsigned _typeA,
                                       const types::t_unsigned _typeC,
                                       const types::t_real _gamma, const types::t_real _sigma, 
                                       const types::t_real _i[5] )
  {
    types::t_unsigned offset = _typeA+_typeC;
    if ( gammas.size() < offset + 1  )
      gammas.resize( offset + 1, types::t_real(0) );
    if ( sigmas.size() < offset + 1  )
      sigmas.resize( offset + 1, types::t_real(0) );
    if ( betas.size() < (offset + 1) * 5 )
      betas.resize( (offset+1)*5, types::t_real(0) );
    gammas[offset] = _gamma;
    sigmas[offset] = _sigma;
    const types::t_real *i_beta = _i;
    std::copy( i_beta, i_beta+5, betas.begin() + offset * 5 );
  }

  bool Atomic_Functional :: load( const TiXmlElement& _node,
                                  const types::t_unsigned &_site_index,
                                  const types::t_unsigned &_type_index,
                                  const Crystal::Lattice::t_Site &_othersite,
                                  Crystal :: Structure &_structure )
  {
    __ASSERT( std::string( _node.Value() ).compare("Functional") != 0;, "Incorrect node.\n" )
    __ASSERT( std::string( _node.Attribute( "type" ) ).compare("Vff") != 0;,
              "Incorrect node.\n" )

    const Crystal::Lattice &lattice = *Crystal::Structure::lattice;
    __ASSERT( _site_index >= lattice.sites.size(), "Index out of range.\n" )
    const Crystal::Lattice &lsite = lattice.sites[ _site_index ];
    __ASSERT( _type >= lsite.type.size(), "Index out of range.\n" )

    structure = &_structure;
    specie = lsite.type[ _type_index ];
    site = _site_index;
    type = _type_index;
    const types::t_unsigned nbonds( _othersite.type.size() );
    const types::t_unsigned nangles( 2*nbonds - 1 );

    // Resize parameter arrays.
    length.clear(); lengths.resize( nbonds, 0e0 );     
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
                  << type << "-" << (*i_atype) << ".\n" )
      __DOASSERT( not child->Attribute("d0"), 
                     "Could not find parameters d0 of bond: "
                  << type << "-" << (*i_atype) << ".\n" )
      __DOASSERT( not child->Attribute("alpha"), 
                     "Could not find parameters alpha of bond: "
                  << type << "-" << (*i_atype) << ".\n" )

      d0 = boost::lexical_cast< types::t_real >( child->Attribute( "d0" ) );
      alphas[ 5 * bond_index ]
        = boost::lexical_cast< types::t_real >( child->Attribute( "alpha" ) );
      for( size_t i(3); i < 7; ++i )
      {
        const std::string str(   std::string("alpha") 
                               + boost::lexical_cast<std::string>(i) );

        if( not child->Attribute(str) ) continue;

        alphas[ 5 * bond_index + i - 2 ]
          = boost::lexical_cast< types::t_real >( child->Attribute(str) );
      }
      
      // then the angle and bond-angle parameters.
      Crystal::Lattice::t_Site::t_Type::const_iterator i_btype = i_atype;
      for(; i_bsite != i_asite_end; ++i_bsite, ++angle_index )
      {
        child = find_angle_node( _node, specie, *i_atype, *i_btype );
        
        __DOASSERT( not child,
                       "Could not find parameters for bond: "
                    << type << "-" << (*i_atype) << "-" << (*i_btype) << ".\n" )
        __DOASSERT( not child->Attribute("gamma"), 
                       "Could not find parameters gamma of bond: "
                    << type << "-" << (*i_atype) << "-" << (*i_btype) << ".\n" )
        __DOASSERT( not child->Attribute("beta"), 
                       "Could not find parameters beta of bond: "
                    << type << "-" << (*i_atype) << "-" << (*i_btype) << ".\n" )

        // find gamma
        std::string str( child->Attribute("gamma") );
        if ( str.compare("tet") == 0 or str.compare("tetrahedral") == 0  )
          gammas[angle_index] = -0.33333333333333333333333333333333333333333;
        else gammas[angle_index] = boost::lexical_cast<types::t_real>( str );
        __DOASSERT( std::abs(gamma) > 1,
                    " gamma must be comprised between 1 and -1.\n" )

        // find sigma
        if( child->Attribute( "sigma" ) )
          sigmas[angle_index] = boost::lexical_cast<types::t_real>( child->Attribute( "sigma" ) );
        // find betas
        betas[ 5 * angle_index ] 
           = boost::lexical_cast<types::t_real>( child->Attribute( "betas" ) );
        for( size_t i(3); i < 7; ++i )
        {
          const std::string str(   std::string("betas") 
                                 + boost::lexical_cast<std::string>(i) );
        
          if( not child->Attribute(str) ) continue;
        
          betas[ 5 * angle_index + i - 2 ]
            = boost::lexical_cast< types::t_real >( child->Attribute(str) );
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
      const TiXmlElement *child( result );
      size_t partial_match(0);
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
        if( A.compare( _A ) == 0 and B.compare( _B ) == 0 ) return child;
        // partial match.
        if( B.compare( _A ) == 0 and A.compare( _B ) == 0 ) 
          { result = child; ++partial_match; }
      }

      if( partial_match < 2 ) return result;

      std::cerr << "Found more than one acceptable Bond tag for "
                << _A << "-" << _B
                << ".\n Aborting\n.";

      return NULL;
    }

    const TiXmlElement* find_bond_node( const TiXmlElement &_node, 
                                        const std::string &_C, 
                                        const std::string &_A, 
                                        const std::string &_B )
    {
      const TiXmlElement *result = NULL;
      const TiXmlElement *child( result );
      for(; child; child = child->NextSiblingElement( "Bond" ) )
      {
        if( not (      child->Attribute( "A" )
                   and child->Attribute( "B" ) 
                   and child->Attribute( "C" ) ) ) 
        {
          std::cerr << "Found incomplete Angle tag in xml input.\nTag is ignored.\n";
          continue;
        }
        const std::string A = child->Attribute("A");
        const std::string B = child->Attribute("B");
        const std::string C = child->Attribute("C");
        // no match
        if( C.compare( _C ) != 0 ) continue;
        // exact match.
        if( A.compare( _A ) == 0 and B.compare( _B ) == 0 ) return child;
        // partial match.
        if( B.compare( _A ) == 0 and A.compare( _B ) == 0 ) 
          { result = child; ++partial_match; }
      }

      if( partial_match < 2 ) return result;

      std::cerr << "Found more than one acceptable Angle tag for "
                << _A << "-" << _C << "-" << _B
                << ".\n Aborting\n.";

      return NULL;
    } 
  } // end of details namespace
} // namespace vff
