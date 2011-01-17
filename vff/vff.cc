#include "LaDaConfig.h"

#include <cstdlib>
#include <algorithm>
#include <functional>
#include <iomanip>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#ifdef LADA_MPI
# include <boost/mpi/collectives/all_reduce.hpp>
#endif

#include <physics/physics.h>
#include <opt/debug.h>
#include <opt/tinyxml.h>
#include <opt/path.h>
#include <opt/mpi.h>
#include <opt/debug.h>


#include "vff.h"
  
namespace LaDa
{
  namespace vff
  { 
    types::t_real const twos3  = 2e0 * std::sqrt(3e0); 
    types::t_real const one16  = 1e0 / 16e0;
    types::t_real const s3o160 = std::sqrt(3e0) / 8e0;
    types::t_real const one640 = 1e0 / 640e0;
    types::t_real const three8 = 3e0 / 8e0;
    types::t_real const s33o8  = 3e0 * std::sqrt(3) / 8e0;
    types::t_real const s33o16 = 3e0 * std::sqrt(3) / 16e0;
    types::t_real const thre16 = 3e0 / 16e0;
    types::t_real const thre32 = 3e0 / 32e0;
    types::t_real const s33128 = 3e0 / 128e0 * std::sqrt(3);
    types::t_real const s33256 = 3e0 / 256 * std::sqrt(3);
    types::t_real const no1280 = 0.00703125;
    types::t_real const no2560 = 0.003515625;
//   // 2*sqrt(3)
//   const types::t_real twos3 = 3.4641016151377545870548926830117447338856105076207612561116139588;
//   // 1 / 16
//   const types::t_real one16 = 0.0625;
//   // sqrt(3) / 8
//   const types::t_real s3o160 = 0.2165063509461096616909307926882340458678506567262975785069758724;
//   // 1 / 640 
//   const types::t_real one640 = 0.0015625;
//   // 3 / 8
//   const types::t_real three8 = 0.375;
//   // 3 * sqrt(3) / 8
//   const  types::t_real s33o8  = 0.6495190528383289850727923780647021376035519701788927355209276172;
//   // 3 * sqrt(3) / 16
//   const  types::t_real s33o16 = 0.3247595264191644925363961890323510688017759850894463677604638086;
//   // 3 / 16
//   const  types::t_real thre16 = 0.1875; 
//   // 3 / 32
//   const  types::t_real thre32 = 0.09375;
//   // 3/128 *sqrt(3) 
//   const  types::t_real s33128 = 0.0405949408023955615670495236290438836002219981361807959700579760;
//   // 3/256 *sqrt(3) 
//   const  types::t_real  s33256 = 0.0202974704011977807835247618145219418001109990680903979850289880;
//   const  types::t_real  no1280 = 0.00703125;
//   const  types::t_real  no2560 = 0.003515625;

    types::t_real Vff :: energy() const
    {
      types::t_real energy = 0;
      
      LADA_MPI_SPLIT_LOOP( t_Centers :: const_iterator, center, centers_, comm )
      for (; i_center != i_center_end; ++i_center) energy += evaluate_center( *i_center );
      LADA_MPI_CODE
      ( 
        energy = boost::mpi::all_reduce( comm, energy, std::plus<types::t_real>() ); 
      )

      return energy;
    }

    types::t_real Vff :: evaluate_center(const AtomicCenter &_center) const
    {
      types :: t_real energy = 0;
      types :: t_real scale2 = structure.scale * structure.scale;

      AtomicCenter :: const_iterator i_bond_begin  = _center.begin();
      AtomicCenter :: const_iterator i_bond_end    = _center.end();
      AtomicCenter :: const_iterator i_bond ( i_bond_begin );
      for( ; i_bond != i_bond_end; ++i_bond )
      {
        // Bond parameters.
        BondData const &bond0_param = get_bond_data(i_bond);
        types::t_real const bond0_length = bond0_param.length;
        
        // computes e0 for bond-angle now 
        types::t_real e0 = i_bond.norm2() * scale2 / bond0_length - bond0_length; 
        // avoids double counting
        if(_center.index < i_bond->index)
        {
          types::t_real const * i_alpha = bond0_param.alphas + (max_vff_expansion - 1); 
          types::t_real dummy = e0 * one640 * (*i_alpha); --i_alpha;
          dummy =   e0 * ( (*i_alpha) * s3o160 + dummy ); --i_alpha;
          dummy =   e0 * ( (*i_alpha) * one16 + dummy ); --i_alpha;
          dummy =   e0 * ( (*i_alpha) / twos3 + dummy ); --i_alpha;
          energy += e0 * e0 * ( (*i_alpha)  + dummy ); 
        }

        // Three body terms
        AtomicCenter :: const_iterator i_angle ( i_bond_begin );
        for (; i_angle != i_bond_end; ++i_angle )
        {
          if ( i_angle == i_bond ) continue;
          
          // sets up parameters
          BondData const &bond1_param = get_bond_data(i_angle);
          AngleData const &angle_param = get_angle_data(i_bond, i_angle);
          types::t_real const bond1_length = bond1_param.length;
          types::t_real const end_length = std::sqrt( bond0_length * bond1_length );
          types::t_real const gamma = angle_param.gamma;
          
          // Bond bending
          types::t_real e1 =   i_bond.scalar_product( i_angle ) * scale2 / end_length 
                             - end_length * gamma;
          if ( i_bond - i_angle > 0 ) // avoids double counting.
          {
            types::t_real const *i_beta = angle_param.betas + (max_vff_expansion - 1);
            types::t_real dummy = e1 * one640 * (*i_beta); --i_beta;
            dummy =   e1 * ( (*i_beta) * s3o160 + dummy ); --i_beta;
            dummy =   e1 * ( (*i_beta) * one16 + dummy ); --i_beta;
            dummy =   e1 * ( (*i_beta) / twos3 + dummy ); --i_beta;
            energy +=   e1 * e1 * ( (*i_beta)  + dummy ); 
          }

          // Bond angle
          energy += e1 * e0 * angle_param.sigma;
        }

      } // for( ; i_bond != i_bond_end; ++i_bond )

      return energy * three8;
    }

    types::t_real Vff :: evaluate_center_with_gradient( const AtomicCenter &_center,
                                                        const math::rMatrix3d &_strain,
                                                        math::rMatrix3d &_stress,
                                                        const math::rMatrix3d &_K0 ) const
    {
      types :: t_real energy = 0;
      types :: t_real scale2 = structure.scale * structure.scale;

      AtomicCenter :: const_iterator i_bond_begin  = _center.begin();
      AtomicCenter :: const_iterator i_bond_end    = _center.end();
      AtomicCenter :: const_iterator i_bond ( i_bond_begin );
      for( ; i_bond != i_bond_end; ++i_bond )
      {
        // Bond parameters.
        BondData const &bond0_param = get_bond_data(i_bond);
        types::t_real const bond0_length = bond0_param.length;
        math::rVector3d d0; i_bond.vector( d0 );
        
        // computes e0 for bond-angle now 
        types::t_real e0 = d0.squaredNorm() * scale2 / bond0_length - bond0_length;
        // adds energy only if bond has not already been counted.
        if(_center.index < i_bond->index)
        {
          // energy 
          types::t_real const* i_alpha = bond0_param.alphas + (max_vff_expansion - 1); 
          types::t_real dummy = e0 * one640 * (*i_alpha); --i_alpha;
          dummy = e0 * ( (*i_alpha) * s3o160 + dummy ); --i_alpha;
          dummy = e0 * ( (*i_alpha) * one16 + dummy ); --i_alpha;
          dummy = e0 * ( (*i_alpha) / twos3 + dummy ); --i_alpha;
          energy += e0 * e0 * ( (*i_alpha)  + dummy ); 
          
          // then gradient 
          i_alpha += 4; // alphas.begin() + 5 * bond_kind;
          dummy = e0 * no1280 * (*i_alpha); --i_alpha;
          dummy = e0 * ( (*i_alpha) * s33128 + dummy ); --i_alpha;
          dummy = e0 * ( (*i_alpha) * thre16 + dummy ); --i_alpha;
          dummy = e0 * ( (*i_alpha) * s33o8 + dummy ); --i_alpha;
          types::t_real e0grad =   types::t_real(2.0) * scale2 / bond0_length 
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
            BondData const &bond1_param = get_bond_data(i_angle);
            AngleData const &angle_param = get_angle_data(i_bond, i_angle);
            types::t_real const bond1_length = bond1_param.length;
            types::t_real const mean_length = std::sqrt( bond0_length * bond1_length );
            types::t_real const gamma = angle_param.gamma;
            types::t_real const sigma = angle_param.sigma;
            math::rVector3d d1; i_angle.vector( d1 );
            
            // Bond bending
            types::t_real e1 = d0.dot(d1) * scale2 / mean_length - mean_length * gamma;
            if ( i_bond - i_angle > 0 )
            {
              // energy
              types::t_real const *i_beta = angle_param.betas + (max_vff_expansion - 1);
              types::t_real dummy = e1 * one640 * (*i_beta); --i_beta;
              dummy = e1 * ( (*i_beta) * s3o160 + dummy ); --i_beta;
              dummy = e1 * ( (*i_beta) * one16 + dummy ); --i_beta;
              dummy = e1 * ( (*i_beta) / twos3 + dummy ); --i_beta;
              energy +=   e1 * e1 * ( (*i_beta)  + dummy ); 
              
              // then gradient 
              i_beta += (max_vff_expansion - 1); // goes back to beginning of array
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
              math::rVector3d hold0 = types::t_real(1.5) * e1 * sigma / bond0_length * scale2
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
                                  * (   types::t_real(2) * e1 / bond0_length * d0[i] * d0[k] 
                                      +   e0 / mean_length 
                                        * (d0[i] * d1[k] + d1[i] * d0[k]) );
          } // angle loop

      } // for( ; i_bond != i_bond_end; ++i_bond )

      return energy * three8;
    }

    types::t_real Vff :: micro_strain(const AtomicCenter &_center) const
    {
      if ( _center.size() != 4 )
        BOOST_THROW_EXCEPTION(input() << error_string("Wrong number of bonds in structure."));

      // computes the microstrain as the average strain over four paralellepipeds.
      // first computes four bond vectors, bond lengths, and angles.
      math::rVector3d bonds[4];
      types::t_real l0s[4]; 
      types::t_real angles[4][4]; 
      AtomicCenter :: const_iterator i_bond  = _center.begin();
      AtomicCenter :: const_iterator i_bond_end    = _center.end();
      for( size_t i(0); i_bond != i_bond_end; ++i_bond, ++i )
      {
        bonds[i] = structure.scale *  i_bond.vector( bonds[i] );
        l0s[i]   = get_bond_data(i_bond).length;

        angles[i][i] = 0;
        AtomicCenter :: const_iterator i_angle( i_bond ); ++i_angle;
        for( size_t j(i+1); i_angle != i_bond_end; ++i_angle, ++j )
        {
          if( i == j ) {  continue; }
          angles[i][j] = get_angle_data(i_bond, i_angle).gamma;
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

    void Vff :: print_escan_input( const t_Path &_f ) const
    {
      namespace bfs = boost::filesystem;
      std::ostringstream stream;
      types::t_unsigned nb_pseudos=0;

      // prints cell vectors in units of a0 and other
      // whatever nanopes other may be
      for( types::t_unsigned i = 0; i < 3; ++i )
        stream << std::fixed << std::setprecision(7) 
               << std::setw(18) << std::setprecision(8)
               << structure.cell(0,i) * structure.scale / Physics::a0("A")
               << std::setw(18) << std::setprecision(8)
               << structure.cell(1,i) * structure.scale / Physics::a0("A")
               << std::setw(18) << std::setprecision(8)
               << structure.cell(2,i) * structure.scale / Physics::a0("A")
               << std::setw(18) << std::setprecision(8) << structure.cell(0,i) 
               << std::setw(18) << std::setprecision(8) << structure.cell(1,i) 
               << std::setw(18) << std::setprecision(8) << structure.cell(2,i) << "\n";

      // prints atomic position, strain, weight, and atomic position in
      // "other unit"
      t_Centers :: const_iterator i_center = centers_.begin();
      t_Centers :: const_iterator i_end = centers_.end();
      for(; i_center != i_end; ++i_center )
      {
        // first gets pseudo index
        size_t const str_index0 = i_center->get_index();
        LADA_BASSERT( str_index0 < structure.atoms.size(),
                      internal() << error_string("Atom index is incorrect in first-neighbor tree."));

        t_Atom const & atom0 = structure.atoms[str_index0];
        types::t_unsigned const index0 = Physics::Atomic::Z( atom0.type );
        types::t_real msstrain = micro_strain(*i_center);

        // finally goes over bonds and finds number of pseudos and their
        // weights
        AtomicCenter :: const_iterator i_bond = i_center->begin();
        AtomicCenter :: const_iterator i_bond_end = i_center->end();
        typedef std::pair<types::t_unsigned, types::t_unsigned > t_pseudo;
        typedef std::vector< t_pseudo > t_pseudos;
        t_pseudos pseudos;
        for(; i_bond != i_bond_end; ++i_bond )
        { 
          size_t const str_index1 = i_bond->get_index();
          LADA_BASSERT( str_index1 < structure.atoms.size(),
                        internal() << error_string("Atom index is incorrect in first-neighbor tree.") );
          t_Atom const & atom1 = structure.atoms[str_index1];
          types::t_unsigned const index1 = Physics::Atomic::Z( atom1.type );

          t_pseudos::iterator i_pseudo = pseudos.begin();
          t_pseudos::iterator i_pseudo_end = pseudos.end();
          for(; i_pseudo != i_pseudo_end; ++i_pseudo )
            if ( i_pseudo->first == index1 ) break;

          if ( i_pseudo == i_pseudo_end ) pseudos.push_back( t_pseudo( index1, 1 ) ); 
          else  ++(i_pseudo->second); 
        }

        // now goes over found pseudos and creates output
        t_pseudos::const_iterator i_pseudo = pseudos.begin();
        t_pseudos::const_iterator i_pseudo_end = pseudos.end();
        for( ; i_pseudo != i_pseudo_end; ++i_pseudo )
        {
          math::rVector3d pos = (!structure.cell) * i_center->atom().pos;
          for( size_t i(0); i < 3; ++i ) pos[i] -= std::floor( pos[i] + 0.000001 );
          ++nb_pseudos;
          stream << std::fixed    << std::setprecision(7)
                 << std::setw(6)  << index0 << '0'
                    << std::setw(2) << std::setfill('0') 
                    << std::right << i_pseudo->first  
                    << " " << std::setfill(' ') // pseudo index
                 << std::setw(12) << pos[0] << " "  // pseudo position
                 << std::setw(12) << pos[1] << " "  
                 << std::setw(12) << pos[2] << " "  
                 << std::setw(18) << msstrain << " " // microscopic strain
                 << std::setw(6) << std::setprecision(2)
                    << types::t_real( i_pseudo->second ) * 0.25  << " " // weight
                 << std::setw(18) << std::setprecision(7) << pos[0] << " " // pseudo position
                 << std::setw(12) << pos[1] << " "
                 << std::setw(12) << pos[2] << "\n";
        }

      }
      const t_Path directory( _f.parent_path() );
      LADA_TRY_BEGIN
        if( not ( directory.empty() or bfs::exists( directory ) ) )
          bfs::create_directory( directory );
        std::ofstream file( _f.string().c_str(), std::ios_base::out|std::ios_base::trunc ); 
        LADA_DO_NASSERT( file.bad(), "Could not open file " << _f << ".\n" ) 
        // prints number of atoms
        file << nb_pseudos << "\n";
        // print rest of file
        file << stream.str();
        file.flush();
        file.close();
        LADA_NASSERT( not bfs::exists( _f ), _f << " was not created.\n" )
      LADA_TRY_END(, "")
    }

#   ifdef LADA_DEBUG
      void Vff :: check_tree() const
      {
        t_Centers :: const_iterator i_center = centers_.begin();
        t_Centers :: const_iterator i_center_end = centers_.end();
        for(size_t index(0); i_center != i_center_end; ++i_center, ++index )
        {
          LADA_DO_NASSERT( not i_center->i_atom_, 
                      "Origin of the center is invalid\n"; )
          LADA_DO_NASSERT( not i_center->structure, 
                      "Invalid pointer to structure\n"; )
          LADA_DO_NASSERT( i_center->bonds.size() != 4,
                         "Invalid number of bonds: "
                      << i_center->bonds.size() << ", " << index << "\n"; )
          LADA_DO_NASSERT( i_center->translations.size() != 4,
                         "Invalid number of translations: "
                      << i_center->translations.size() << "\n"; )
        }
      }
#   endif

  } // namespace vff
} // namespace LaDa
