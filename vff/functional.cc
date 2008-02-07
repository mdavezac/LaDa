//
//  Version: $Id$
//
#include <cstdlib>

#include <algorithm>
#include <functional>

#include <physics/physics.h>
#include <opt/ndim_iterator.h>
#include <opt/debug.h>
#include <opt/atat.h>

#include "functional.h"
  
namespace Vff
{ 
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

  bool Functional :: initialize_centers()
  {
    centers.clear();
    
    // Creates a list of centers
    t_Atoms :: iterator i_atom = structure.atoms.begin();
    t_Atoms :: iterator i_atom_end = structure.atoms.end();

    for(types::t_unsigned index=0; i_atom != i_atom_end; ++i_atom, ++index )
      centers.push_back( Atomic_Center( structure, *i_atom, index ) );

    // Creates a list of closest neighbors
    std::vector< atat::rVector3d > neighbors;
    Ising_CE::Lattice::t_Sites :: iterator i_site_begin = structure.lattice->sites.begin();
    Ising_CE::Lattice::t_Sites :: iterator i_site, i_site2;
    Ising_CE::Lattice::t_Sites :: iterator i_site_end = structure.lattice->sites.end();
    
    for(i_site = i_site_begin; i_site != i_site_end; ++i_site )
    {
      for(i_site2 = i_site_begin; i_site2 != i_site_end; ++i_site2 )
      {
        opt::Ndim_Iterator< types::t_int, std::less_equal<types::t_int> > period;
        period.add(-1,1);
        period.add(-1,1);
        period.add(-1,1);
        do // goes over periodic image
        {
          // constructs perdiodic image of atom *i_bond
          atat::rVector3d frac_image, image;
          frac_image[0] =  (types::t_real) period.access(0);
          frac_image[1] =  (types::t_real) period.access(1);
          frac_image[2] =  (types::t_real) period.access(2);
          image = i_site2->pos + structure.lattice->cell * frac_image;
          if( atat::norm2( image - i_site->pos ) > types::tolerance )
            neighbors.push_back( image - i_site->pos );
        } while ( ++period ); 
      }
    }

    // Sorts the neighbors according to distance from origin
    std::sort( neighbors.begin(), neighbors.end(), atat::norm_compare() );
    // And reduces to first neighbors only
    neighbors.resize(4*structure.lattice->sites.size());

    t_Centers :: iterator i_begin = centers.begin();
    t_Centers :: iterator i_end = centers.end();
    t_Centers :: iterator i_center, i_bond;
    atat::rVector3d image;
    atat::rVector3d frac_image;
    atat::rVector3d cut;
    cut = neighbors.front();
    types::t_real cutoff = 0.25 * atat::norm2( neighbors.front() );
    for( i_center = i_begin; i_center != i_end; ++i_center )
    {
      for( i_bond = i_begin; i_bond != i_end; ++i_bond)
        if( i_bond != i_center )
        {
          std::vector<atat::rVector3d> :: const_iterator i_neigh = neighbors.begin();
          std::vector<atat::rVector3d> :: const_iterator i_neigh_end = neighbors.end();
          for(; i_neigh != i_neigh_end; ++i_neigh )
          {
            image = i_center->origin->pos - *i_neigh - i_bond->origin->pos;
            frac_image = (!structure.cell) * image;
            cut[0] = frac_image[0] - rint( frac_image[0] );
            cut[1] = frac_image[1] - rint( frac_image[1] );
            cut[2] = frac_image[2] - rint( frac_image[2] );
            cut = structure.cell * cut;
            if( atat::norm2( cut ) < cutoff )
            {
              i_center->bonds.push_back( t_Center ::__make__iterator__(i_bond) );
              i_center->translations.push_back( frac_image );
              i_center->do_translates.push_back( atat::norm2(frac_image) > atat::zero_tolerance );
            }
          }
        }
    }

#ifdef _DEBUG
    check_tree();
#endif
    return true;
  } // Functional :: construct_bond_list

  bool Functional :: construct_centers()
  {
    centers.clear();
    
    t_Atoms :: iterator i_atom = structure.atoms.begin();
    t_Atoms :: iterator i_atom_end = structure.atoms.end();

    for(types::t_unsigned index=0; i_atom != i_atom_end; ++i_atom, ++index )
      centers.push_back( Atomic_Center( structure, *i_atom, index ) );

    t_Centers :: iterator i_begin = centers.begin();
    t_Centers :: iterator i_end = centers.end();
    t_Centers :: iterator i_center, i_bond;

    for( i_center = i_begin; i_center != i_end; ++i_center )
    {
      types::t_int nb_bonds;
      
      for( i_bond = i_begin; i_bond != i_end; ++i_bond )
        if ( i_bond != i_center )
          nb_bonds = i_center->add_bond ( t_Center::__make__iterator__(i_bond), bond_cutoff);
    } // loop over atomic centers

    
    // consistency check
    for( i_center = i_begin; i_center != i_end; ++i_center )
      if ( i_center->size() != 4 )
      {
        std::cerr << " Atomic center at " << (atat::rVector3d) i_center->Origin()
                  << " has " << i_center->size() 
                  << " bonds!!" << std::endl;
        return false;
      } 

#ifdef _DEBUG
    check_tree();
#endif
    return true;
  } // Functional :: construct_bond_list

  bool Functional :: Load ( const TiXmlElement &_element )
  {
    // some consistency checking
    if ( structure.lattice->get_nb_sites() != 2 )
    { 
      std::cerr << "Cannot do vff on this lattice" 
                << std::endl
                << "Need 2 and only 2 different sites per unit cell" 
                << std::endl;
      return false;
    }
    if (     structure.lattice->get_nb_types(0) + structure.lattice->get_nb_types(1) < 2
         and structure.lattice->get_nb_types(0) + structure.lattice->get_nb_types(1) > 4 )
    { 
      std::cerr << "Cannot do vff on this lattice" 
                << std::endl
                << "Need at two sites with at most two different atomic types" 
                << std::endl;
      return false;
    }

    const TiXmlElement* parent = find_node( _element );
    if( parent ) return Load_(*parent);
    std::cerr << "Could not find an <Functional type=\"vff\"> tag in input file" 
              << std::endl;
    return false;
  }

  bool Functional :: Load_( const TiXmlElement &_element )
  {
    const TiXmlElement *child;
    std::string str;


    // reads and initializes bond cutoff
    _element.Attribute( "cutoff", &bond_cutoff );
    if ( bond_cutoff == 0 )
      bond_cutoff = 1.25; 
    bond_cutoff *= std::sqrt(3.0) / 4.0; // axes bs same as CE
    bond_cutoff *= bond_cutoff; // squared for simplicity
    
    // creates an unitialized array of atomic functionals
    functionals.clear();
    functionals.push_back( Atomic_Functional( structure.lattice->get_atom_string(0,0),
                           structure, 0, 0) );
    if ( structure.lattice->get_nb_types(0) == 2 )
      functionals.push_back( Atomic_Functional( structure.lattice->get_atom_string(0,1),
                             structure, 0, 1) );
    functionals.push_back( Atomic_Functional( structure.lattice->get_atom_string(1,0),
                           structure, 1, 0) );
    if ( structure.lattice->get_nb_types(1) == 2 )
      functionals.push_back( Atomic_Functional( structure.lattice->get_atom_string(1,1),
                             structure, 1, 1) );

    // **************************************************************
    // first reads bond interactions
    // **************************************************************
    child = _element.FirstChildElement( "Bond" );
    for( types::t_unsigned i=0; child and i < functionals.size();
         child = child->NextSiblingElement( "Bond" ), ++i )
    {
      // reads input
      std::string A, B;
      types::t_real d0;
      std::vector<types::t_real> alphas(5, types::t_real(0));
      if (    not child->Attribute("A") 
           or not child->Attribute("B") 
           or not child->Attribute("d0") 
           or not child->Attribute("alpha") )
      {
        std::cerr << "Bond input is incomplete in input file"
                  << std::cerr;
        return false;
      }
      A = child->Attribute("A");
      B = child->Attribute("B");
      child->Attribute("d0", &d0);
      child->Attribute("alpha",  &(alphas[0]));
      child->Attribute("alpha3", &(alphas[1]));
      child->Attribute("alpha4", &(alphas[2]));
      child->Attribute("alpha5", &(alphas[3]));
      child->Attribute("alpha6", &(alphas[4]));

      // finds out where to put it
      types::t_int where[2];
      bond_indices( A, B, where );

      functionals[ where[0] ].add_bond( where[1], d0, alphas );
      functionals[ where[1]+structure.lattice->get_nb_types(0)].add_bond( where[0], d0, alphas );
    }

    // **************************************************************
    // then reads angle and bond-angle interactions 
    // **************************************************************
    child = _element.FirstChildElement( "Angle" );
    for(; child; child = child->NextSiblingElement( "Angle" ) )
    {
      // reads input
      std::string A, B, C;
      types::t_real sigma, gamma;
      std::vector<types::t_real> betas(5,types::t_real(0));
      if (    not child->Attribute("A") 
           or not child->Attribute("B") 
           or not child->Attribute("C") 
           or not child->Attribute("gamma") 
           or not child->Attribute("sigma") 
           or not child->Attribute("beta") )
      {
        std::cerr << "Angle input is incomplete in input file"
                  << std::cerr;
        return false;
      }
      A = child->Attribute("A");
      B = child->Attribute("B");
      C = child->Attribute("C");
      child->Attribute("sigma",  &sigma);
      child->Attribute("beta",  &(betas[0]));
      child->Attribute("beta3", &(betas[1]));
      child->Attribute("beta4", &(betas[2]));
      child->Attribute("beta5", &(betas[3]));
      child->Attribute("beta6", &(betas[4]));
      // gamma is somewhat more complicated...
      str = child->Attribute("gamma");
      if ( str.compare("tet") == 0 or str.compare("tetrahedral") == 0  )
        gamma = -0.33333333333333333333333333333333333333333;
      else
        child->Attribute("gamma", &gamma);
      if( std::abs(gamma) > 1 )
      {
        std::cerr << " gamma must be comprised between 1 and -1 " << std::endl;
        return false;
      }
      
      // finds out where to put it
      types::t_int where[3];
      angle_indices( A, B, C, where );
      functionals[ where[1] ].add_angle( where[0], where[2], gamma, sigma, betas );
    }

    return true;
  }  // Functional :: Load

  const TiXmlElement* Functional :: find_node ( const TiXmlElement &_element )
  {
    const TiXmlElement *parent;
    std::string str;

    // This whole section tries to find a <Functional type="vff"> tag
    // in _element or its child
    str = _element.Value();
    if ( str.compare("Functional" ) != 0 )
      parent = _element.FirstChildElement("Functional");
    else parent = &_element;
    
    while (parent)
    {
      str = "";
      if ( parent->Attribute( "type" )  )
        str = parent->Attribute("type");
      if ( str.compare("vff" ) == 0 )
        break;
      parent = parent->NextSiblingElement("Functional");
    }
    if ( parent ) return parent;
    
    std::cerr << "Could not find an <Functional type=\"vff\"> tag in input file" 
              << std::endl;
    return NULL;
  }  // Functional :: find_node

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
        std::vector<types::t_real>::const_iterator i_alpha = alphas.begin() + 5 * bond_kind;
        energy +=   e0 * e0 * ( (*i_alpha) 
                  + e0 * ( *(++i_alpha) / twos3
                  + e0 * ( *(++i_alpha) * one16
                  + e0 * ( *(++i_alpha) * s3o160
                  + e0 * ( *(++i_alpha) * one640 )))));
      } 

      // Three body terms
      Atomic_Center :: const_iterator i_angle ( i_bond_begin );
      for (; i_angle != i_bond_end; ++i_angle )
        if ( i_angle != i_bond )
        {
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
            std::vector< types::t_real > :: const_iterator i_beta = betas.begin() + 5 * angle_kind;
            energy +=   e1 * e1 * ( (*i_beta) 
                      + e1 * ( *(++i_beta) / twos3 
                      + e1 * ( *(++i_beta) * one16
                      + e1 * ( *(++i_beta) * s3o160
                      + e1 * ( *(++i_beta) * one640 )))));
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
        std::vector<types::t_real>::const_iterator i_alpha = alphas.begin() + 5 * bond_kind;
        energy +=  e0 * e0 * ( (*i_alpha) 
                   + e0 * ( *(++i_alpha) / twos3
                   + e0 * ( *(++i_alpha) * one16
                   + e0 * ( *(++i_alpha) * s3o160
                   + e0 * ( *(++i_alpha) * one640 )))));

        // then gradient 
        i_alpha -= 4; // alphas.begin() + 5 * bond_kind;
        types::t_real e0grad = 2.0 * scale2 / bond_length *
                                 e0 * ( *(  i_alpha) * 1.5e0
                               + e0 * ( *(++i_alpha) * s33o8
                               + e0 * ( *(++i_alpha) * thre16
                               + e0 * ( *(++i_alpha) * s33128
                               + e0 * ( *(++i_alpha) * no1280 )))));
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
            std::vector< types::t_real > :: const_iterator i_beta = betas.begin() + 5 * angle_kind;
            energy +=   e1 * e1 * ( (*i_beta) 
                      + e1 * ( *(++i_beta) / twos3
                      + e1 * ( *(++i_beta) * one16
                      + e1 * ( *(++i_beta) * s3o160
                      + e1 * ( *(++i_beta) * one640 )))));
            
            // then gradient 
            i_beta -= 4; // goes back to beginning of array
            types::t_real e1grad = 2.0 * scale2 / mean_length *
                                     e1 * ( *(  i_beta) * 0.75
                                   + e1 * ( *(++i_beta) * s33o16
                                   + e1 * ( *(++i_beta) * thre32
                                   + e1 * ( *(++i_beta) * s33256
                                   + e1 * ( *(++i_beta) * no2560 )))));
            atat::rVector3d hold0 = e1grad * _strain * d0;
            atat::rVector3d hold1 = e1grad * _strain * d1;
            _center.get_gradient() -= ( hold0 + hold1); // with respect to atomic positions
            i_bond->get_gradient() += hold1; 
            i_angle->get_gradient() += hold0; 
            
            // stress
            for( int i = 0; i < 3; ++i )
              for( int j = 0; j < 3; ++j )
                for( int k = 0; k < 3; ++k )
                  _stress(i,j) += (d1[i] * d0[k] + d0[i] * d1[k]) * _K0(k,j) * e1grad * 0.5;
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
                _stress(i,j) += _K0(k,j) * 0.375 * sigma * scale2 *
                                ( 2.0 * e1 / bond_length * d0[i] * d0[k] 
                                  + e0 / mean_length * (d0[i] * d1[k] + d1[i] * d0[k]) );
        } // angle loop

    } // for( ; i_bond != i_bond_end; ++i_bond )

    return energy * three8;
  }


  Atomic_Center :: Atomic_Center  ( Ising_CE::Structure &_str, t_Atom &_e,
                                    types::t_unsigned _i )
                                 : origin(&_e), structure(&_str)
  {
     is_site_one = ( structure->lattice->get_atom_site_index( _e ) == 0 );
     is_site_one_two_species = ( structure->lattice->get_nb_types( 0 ) == 2 );
     index = _i;
  }

  types::t_unsigned  Atomic_Center :: kind() const
  {
    if ( is_site_one )
      return structure->lattice->convert_real_to_type_index( 0, origin->type );
    else if ( is_site_one_two_species )
      return 2 + structure->lattice->convert_real_to_type_index( 1, origin->type );
    return 1 + structure->lattice->convert_real_to_type_index( 1, origin->type );
  }

  void Functional :: angle_indices( const std::string &_A, const std::string &_B, 
                                    const std::string &_C, types::t_int _indices[3] ) const
  {
    types::t_int siteA, siteB, siteC;
    siteA = structure.lattice->get_atom_site_index( _A );
    if ( siteA == -1 ) goto failure;
    _indices[0] = structure.lattice->get_atom_type_index( _A );
    if ( _indices[0] == -1 ) goto failure;
    siteB = structure.lattice->get_atom_site_index( _B );
    if ( siteB == -1 ) goto failure;
    _indices[1] = structure.lattice->get_atom_type_index( _B );
    if ( _indices[1] == -1 ) goto failure;
    siteC = structure.lattice->get_atom_site_index( _C );
    if ( siteC == -1 ) goto failure;
    _indices[2] = structure.lattice->get_atom_type_index( _C);
    if ( _indices[2] == -1 ) goto failure;

    if ( siteA == siteB or siteA != siteC ) goto failure;
    if ( siteB == 0 ) return;
    
    _indices[1] += structure.lattice->get_nb_types(0);

    return;

failure:
    __THROW_ERROR(   "Something wrong with your input\n" 
                   << "Did not expect angle type " << _A << "-"
                   << _B << "-" << _C << "\n" )
  }

  void Functional :: bond_indices( const std::string &_A, const std::string &_B,
                                   types::t_int _indices[2] ) const
  {
    types::t_int siteA, siteB, swap;
    // finds out where to put it
    siteA = structure.lattice->get_atom_site_index( _A );
    if ( siteA == -1 ) goto failure;
    _indices[0] = structure.lattice->get_atom_type_index( _A );
    if ( _indices[0] == -1 ) goto failure;
    siteB = structure.lattice->get_atom_site_index( _B );
    if ( siteB == -1 ) goto failure;
    _indices[1] = structure.lattice->get_atom_type_index( _B );
    if ( _indices[1] == -1 ) goto failure;

    if ( siteA == siteB ) goto failure;

    // reorders things around
    if( siteA != 1 )  return;
    
    swap = _indices[1];
    _indices[1] = _indices[0];
    _indices[0] = swap;

    return;
failure:
    __THROW_ERROR(   "Something wrong with your input\n" 
                   << "Did not expect bond type " << _A << "-"
                   << _B << "-" << "\n" )
  }

  types::t_unsigned  Atomic_Center :: bond_kind( const Atomic_Center &_bond ) const
  {
    if ( is_site_one ) 
      return  structure->lattice->convert_real_to_type_index( 1, _bond.origin->type );
    else if ( is_site_one_two_species )
      return  structure->lattice->convert_real_to_type_index( 0, _bond.origin->type );
    return 0; 
  }

  types::t_int Atomic_Center :: add_bond( t_BondRefd _bond, 
                                          const types::t_real _cutoff ) 
  {
    bool found_bond = false;
    opt::Ndim_Iterator< types::t_int, std::less_equal<types::t_int> > period;
    period.add(-1,1);
    period.add(-1,1);
    period.add(-1,1);

    do // goes over periodic image
    {
      // constructs perdiodic image of atom *i_bond
      atat::rVector3d frac_image, image;
      frac_image[0] =  (types::t_real) period.access(0);
      frac_image[1] =  (types::t_real) period.access(1);
      frac_image[2] =  (types::t_real) period.access(2);
      image = _bond->origin->pos + structure->cell * frac_image;

      // checks if within 
      if( atat::norm2( image - origin->pos ) < _cutoff  )
      {
        // adds bond
        types :: t_unsigned site = structure->lattice->get_atom_site_index( image );
        if ( (not site) == is_site_one )
          return -1;  // error bond between same site

        bonds.push_back( _bond );
        translations.push_back( frac_image );
        do_translates.push_back( atat::norm2( frac_image ) > atat::zero_tolerance );
        found_bond = true;
      }
    } while ( ++period ); 

    return found_bond ? (types::t_int) bonds.size() : -1; // not a bond
  }


  types::t_real Functional :: energy() const
  {
    types::t_real energy = 0;
    
    t_Centers :: const_iterator i_center = centers.begin();
    t_Centers :: const_iterator i_end = centers.end();
    for (; i_center != i_end; ++i_center)
      energy += functionals[i_center->kind()].evaluate( *i_center );

    return energy;
  }

  // Unpacks opt::Function_Base::variables into Vff::Functional format
  void Functional :: unpack_variables(atat::rMatrix3d& strain)
  {
    __ASSERT( variables->size() == 0, "Too few variables.\n" )
    const_iterator i_x = variables->begin();

    strain(0,0) = ( structure.freeze & Ising_CE::Structure::FREEZE_XX ) ?
                  1.0 : (*i_x++);
    strain(1,1) = ( structure.freeze & Ising_CE::Structure::FREEZE_YY ) ?
                  1.0 : (*i_x++);
    strain(2,2) = ( structure.freeze & Ising_CE::Structure::FREEZE_ZZ ) ?
                  1.0 : (*i_x++);
    strain(0,1) = strain (1,0) = (structure.freeze & Ising_CE::Structure::FREEZE_XY) ?
                                 0.0 : (*i_x++);
    strain(0,2) = strain (2,0) = (structure.freeze & Ising_CE::Structure::FREEZE_XZ) ?
                                 0.0 : (*i_x++);
    strain(2,1) = strain (1,2) = (structure.freeze & Ising_CE::Structure::FREEZE_YZ) ?
                                 0.0 : (*i_x++);

    // compute resulting cell vectors
    structure.cell = strain * structure0.cell;
    unpack_positions( strain, i_x );
  }

  void Functional :: unpack_positions(atat::rMatrix3d& strain,
                                      const_iterator &_i_x )
  {
    // then computes positions
    t_Atoms :: iterator i_atom = structure.atoms.begin();
    t_Atoms :: iterator i_atom_end = structure.atoms.end();
    for(; i_atom != i_atom_end; ++i_atom )
    {
      atat::rVector3d pos;
      if ( not (i_atom->freeze & t_Atom::FREEZE_X ) )
        { pos[0] = 2.0 * (*_i_x ); ++_i_x; }
      else pos[0] = i_atom->pos[0];
      if ( not (i_atom->freeze & t_Atom::FREEZE_Y ) )
        { pos[1] = 2.0 * (*_i_x ); ++_i_x; }
      else pos[1] = i_atom->pos[1];
      if ( not (i_atom->freeze & t_Atom::FREEZE_Z ) )
        { pos[2] = 2.0 * (*_i_x ); ++_i_x; }
      else pos[2] = i_atom->pos[2];

      i_atom->pos = strain * pos;
      __ASSERT( variables->end() - _i_x < 0, "Too few variables.\n" )
    }

    // Correct stray movements of the center of mass
    __ASSERT(    fixed_index[0] < 0
              or fixed_index[1] < 0
              or fixed_index[2] < 0,
              "fixed_index contains negative indices. Was init() called?\n" )
    __ASSERT(    fixed_index[0] >= structure0.atoms.size()
              or fixed_index[1] >= structure0.atoms.size()
              or fixed_index[2] >= structure0.atoms.size(),
              "fixed_index contains out-of-range indices.\n" )
    atat::rMatrix3d cell_inv = !structure.cell;

    types::t_real x =   structure0.atoms[fixed_index[0]].pos[0] 
                      - structure.atoms [fixed_index[0]].pos[0];
    types::t_real y =   structure0.atoms[fixed_index[1]].pos[1] 
                      - structure.atoms [fixed_index[1]].pos[1];
    types::t_real z =   structure0.atoms[fixed_index[2]].pos[2] 
                      - structure.atoms [fixed_index[2]].pos[2];

    if ( Fuzzy::eq(x, 0e0) and Fuzzy::eq(y, 0e0) and Fuzzy::eq(z, 0e0) ) return;
    for(i_atom = structure.atoms.begin(); i_atom != i_atom_end; ++i_atom )
      { i_atom->pos[0] += x; i_atom->pos[1] += y; i_atom->pos[2] += z; }
  }


  // initializes stuff before minimization
  bool Functional :: init()
  {
    // sets up structure0, needed for fractional vs cartesian shit
    structure0 = structure;

    // Now counts the degrees of freedom
    types::t_unsigned dof = 0;
    if ( not (structure0.freeze & Ising_CE::Structure::FREEZE_XX ) ) ++dof;
    if ( not (structure0.freeze & Ising_CE::Structure::FREEZE_XY ) ) ++dof;
    if ( not (structure0.freeze & Ising_CE::Structure::FREEZE_XZ ) ) ++dof;
    if ( not (structure0.freeze & Ising_CE::Structure::FREEZE_YY ) ) ++dof;
    if ( not (structure0.freeze & Ising_CE::Structure::FREEZE_YZ ) ) ++dof;
    if ( not (structure0.freeze & Ising_CE::Structure::FREEZE_ZZ ) ) ++dof;
    dof += posdofs();

    __DOASSERT( not dof, "No degrees of freedom.\n" )
   
    __TRYCODE( function::Base<> :: resize( dof );,
               "Could not resize function.\n" )
    __DOASSERT( not variables, "Could not resize function.\n" )

    strain.zero(); 
    strain(0,0) = 1.0;
    strain(1,1) = 1.0;
    strain(2,2) = 1.0;

    pack_variables(strain);
    unpack_variables(strain);
    
    return true;
  }

  types::t_unsigned Functional :: posdofs()
  {
    fixed_index[0] = -1; fixed_index[1] = -1; fixed_index[2] = -1; 
    types::t_unsigned dof = 0;
    atat::rMatrix3d cell_inv = !structure.cell;
    t_Atoms :: iterator i_atom =  structure0.atoms.begin();
    t_Atoms :: iterator i_atom_end =  structure0.atoms.end();
    for( types::t_unsigned n = 0; i_atom != i_atom_end; ++i_atom, ++n ) 
    {
      if ( not (i_atom->freeze & t_Atom::FREEZE_X ) ) ++dof;
      else if (fixed_index[0] == -1 ) fixed_index[0] = n;
      if ( not (i_atom->freeze & t_Atom::FREEZE_Y ) ) ++dof;
      else if (fixed_index[1] == -1 ) fixed_index[1] = n;
      if ( not (i_atom->freeze & t_Atom::FREEZE_Z ) ) ++dof;
      else if (fixed_index[2] == -1 ) fixed_index[2] = n;
    }
    
    if( fixed_index[0] == -1 ) fixed_index[0] = 0;
    if( fixed_index[1] == -1 ) fixed_index[1] = 0;
    if( fixed_index[2] == -1 ) fixed_index[2] = 0;

    return dof;
  }

  // variables is expected to be of sufficient size!!
  // call init() first
  void Functional :: pack_variables( const atat::rMatrix3d& _strain)
  {
    __ASSERT( variables->size() == 0, "Too few variables\n" )
    // finally, packs vff format into function::Base format
    iterator i_var = begin();
    if ( not (structure0.freeze & Ising_CE::Structure::FREEZE_XX ) )
      { *i_var = _strain(0,0); ++i_var; }
    if ( not (structure0.freeze & Ising_CE::Structure::FREEZE_YY ) )
      { *i_var = _strain(1,1); ++i_var; }
    if ( not (structure0.freeze & Ising_CE::Structure::FREEZE_ZZ ) )
      { *i_var = _strain(2,2); ++i_var; }
    if ( not (structure0.freeze & Ising_CE::Structure::FREEZE_XY ) )
      { *i_var = 0.5*(_strain(1,0) + _strain(0,1)); ++i_var; }
    if ( not (structure0.freeze & Ising_CE::Structure::FREEZE_XZ ) )
      { *i_var = 0.5*(_strain(2,0) + _strain(0,2)); ++i_var; }
    if ( not (structure0.freeze & Ising_CE::Structure::FREEZE_YZ ) )
      { *i_var = 0.5*(_strain(2,1) + _strain(1,2)); ++i_var; }

    pack_positions( i_var );
  }

  void Functional :: pack_positions( iterator &_i_var )
  {
     t_Atoms :: const_iterator i_atom =  structure0.atoms.begin();
     t_Atoms :: const_iterator i_atom_end =  structure0.atoms.end();
     for(; i_atom != i_atom_end; ++i_atom )
     {
       if ( not (i_atom->freeze & t_Atom::FREEZE_X ) )
         { *_i_var = i_atom->pos[0] * 0.5; ++_i_var; }
       if ( not (i_atom->freeze & t_Atom::FREEZE_Y ) )
         { *_i_var = i_atom->pos[1] * 0.5; ++_i_var; }
       if ( not (i_atom->freeze & t_Atom::FREEZE_Z ) )
         { *_i_var = i_atom->pos[2] * 0.5; ++_i_var; }
       __ASSERT( variables->end() - _i_var < 0, "Too few variables.\n" )
     }
  }

  void Functional :: print_escan_input( const std::string &_f ) const
  {
    std::ostringstream stream;
    types::t_unsigned nb_pseudos=0;

    // prints cell vectors in units of a0 and other
    // whatever nanopes other may be
    for( types::t_unsigned i = 0; i < 3; ++i )
      stream << std::fixed << std::setprecision(7) 
             << std::setw(12) << std::setprecision(8)
             << structure.cell(0,i) * structure0.scale / Physics::a0("A")
             << std::setw(12) << std::setprecision(8)
             << structure.cell(1,i) * structure0.scale / Physics::a0("A")
             << std::setw(12) << std::setprecision(8)
             << structure.cell(2,i) * structure0.scale / Physics::a0("A")
             << std::setw(18) << std::setprecision(8) << structure.cell(0,i) 
             << std::setw(12) << std::setprecision(8) << structure.cell(1,i) 
             << std::setw(12) << std::setprecision(8) << structure.cell(2,i) << "\n";

    // prints atomic position, strain, weight, and atomic position in
    // "other unit"
    t_Centers :: const_iterator i_center = centers.begin();
    t_Centers :: const_iterator i_end = centers.end();
    for(; i_center != i_end; ++i_center )
    {
      // first gets pseudo index
      Ising_CE::StrAtom stratom;
      __TRYDEBUGCODE( 
        structure.lattice->convert_Atom_to_StrAtom(
           structure0.atoms[i_center->get_index()], stratom );, 
           "Error while printing escan input for atom " 
        << i_center->get_index() << ": \n"
        << structure0.atoms[i_center->get_index()] << "\n" << structure0 )
        
      types::t_unsigned index = Physics::Atomic::Z( stratom.type );
      types::t_real msstrain = functionals[i_center->kind()]
                                          .MicroStrain( *i_center, structure0 );

      // finally goes over bonds and finds number of pseudos and their
      // weights
      Atomic_Center :: const_iterator i_bond = i_center->begin();
      Atomic_Center :: const_iterator i_bond_end = i_center->end();
      typedef std::pair<types::t_unsigned, types::t_unsigned > t_pseudo;
      typedef std::vector< t_pseudo > t_pseudos;
      t_pseudos pseudos;
      for(; i_bond != i_bond_end; ++i_bond )
      { 
        __TRYDEBUGCODE( 
          structure.lattice->convert_Atom_to_StrAtom( 
            structure0.atoms[i_bond->get_index()], stratom );, 
             "Error while printing escan input for atoms\n" 
          << structure0.atoms[i_bond->get_index()] << "\n" )
        types::t_unsigned Z = Physics::Atomic::Z( stratom.type );
        t_pseudos::iterator i_pseudo = pseudos.begin();
        t_pseudos::iterator i_pseudo_end = pseudos.end();
        for(; i_pseudo != i_pseudo_end; ++i_pseudo )
          if ( i_pseudo->first == Z ) break;

        if ( i_pseudo == i_pseudo_end ) pseudos.push_back( t_pseudo( Z, 1 ) ); 
        else  ++(i_pseudo->second); 
      }

      // now goes over found pseudos and creates output
      t_pseudos::const_iterator i_pseudo = pseudos.begin();
      t_pseudos::const_iterator i_pseudo_end = pseudos.end();
      for( ; i_pseudo != i_pseudo_end; ++i_pseudo )
      {
        atat::rVector3d pos = (!structure.cell) * i_center->Origin().pos;
        ++nb_pseudos;
        stream << std::fixed    << std::setprecision(7)
               << std::setw(6)  << index << '0' << i_pseudo->first  // pseudo index
               << std::setw(12) << pos[0] // pseudo position
               << std::setw(12) << pos[1] 
               << std::setw(12) << pos[2] 
               << std::setw(18) << msstrain << " " // microscopic strain
               << std::setw(6) << std::setprecision(2)
                               << types::t_real( i_pseudo->second ) * 0.25  // weight
               << std::setw(18) << std::setprecision(7) << pos[0] // pseudo position
               << std::setw(12) << pos[1] 
               << std::setw(12) << pos[2] << "\n";
      }

    }
    std::ofstream file( _f.c_str(), std::ios_base::out|std::ios_base::trunc ); 
    // prints number of atoms
    file << nb_pseudos << "\n";
    // print rest of file
    file << stream.str();
    file.flush();
    file.close();
  }

  types::t_real Atomic_Functional
    :: MicroStrain( const Atomic_Center &_center, 
                    const Ising_CE::Structure &_str0 ) const 
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
    stream << "Site " << str << " " << site << std::endl << "  ";
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

  void Functional :: print_out( std::ostream &stream ) const
  {
    t_AtomicFunctionals :: const_iterator i_func = functionals.begin();
    t_AtomicFunctionals :: const_iterator i_func_end = functionals.end();
    for(; i_func != i_func_end; ++i_func )
      i_func->print_out(stream);
  }

#ifdef _DEBUG
  void Functional :: check_tree() const
  {
    t_Centers :: const_iterator i_center = centers.begin();
    t_Centers :: const_iterator i_center_end = centers.end();
    for(; i_center != i_center_end; ++i_center )
    {
      __DOASSERT( not i_center->origin, 
                  "Origin of the center is invalid\n"; )
      __DOASSERT( not i_center->structure, 
                  "Invalid pointer to structure\n"; )
      __DOASSERT( i_center->bonds.size() != 4,
                     "Invalid number of bonds: "
                  << i_center->bonds.size() << "\n"; )
      __DOASSERT( i_center->translations.size() != 4,
                     "Invalid number of translations: "
                  << i_center->translations.size() << "\n"; )
    }
  }
#endif
} // namespace vff

#ifdef _MPI

namespace mpi
{
  //! Allows mpi::BroadCast "-ing" and mpi::AllGather "-ing" of  Vff::Atomic_Functional
  template<>
  bool BroadCast :: serialize<Vff::Atomic_Functional>( Vff::Atomic_Functional &_func )
  {
    if ( not serialize( _func.str ) ) return false;
    if ( not serialize( _func.site ) ) return false;
    if ( not serialize( _func.type ) ) return false;
    if ( not serialize( _func.lengths ) ) return false;
    if ( not serialize( _func.alphas ) ) return false;
    if ( not serialize( _func.betas ) ) return false;
    if ( not serialize( _func.gammas ) ) return false;
    
    return serialize( _func.sigmas );
  }

  //! Allows mpi::BroadCast "-ing" and mpi::AllGather "-ing" of  Vff::Functional
  template<>
  bool BroadCast :: serialize<Vff::Functional>( Vff::Functional &_vff )
  {
    // first serializes bond_cutoff and freeze_none 
    if ( not serialize( _vff.bond_cutoff ) ) return false;
    if ( not serialize( _vff.fixed_index[0] ) ) return false;
    if ( not serialize( _vff.fixed_index[1] ) ) return false;
    if ( not serialize( _vff.fixed_index[2] ) ) return false;

    // finally, serializes atomic functionals
    types::t_int n = _vff.functionals.size();
    if ( not serialize( n ) ) return false;
    if ( stage == COPYING_FROM_HERE )
      _vff.functionals.resize(n, Vff::Atomic_Functional(_vff.structure));

    typedef Vff::Functional::t_AtomicFunctionals :: iterator t_iterator;
    t_iterator i_func = _vff.functionals.begin();
    t_iterator i_func_end = _vff.functionals.end();
    for(; i_func != i_func_end; ++i_func )
      if ( not serialize( *i_func ) ) return false;

    return true;
  }

}

#endif
