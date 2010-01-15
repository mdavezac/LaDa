//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "compare_sites.h"
#include "symmetry_operator.h"


namespace LaDa
{

  namespace Crystal 
  {

    void compose( SymmetryOperator const &_a, SymmetryOperator const &_b, SymmetryOperator &_out )
    {
      _out.op = _a.op * _b.op;
      _out.trans = _a.trans + _a.op * _b.trans;
    }
    
    bool SymmetryOperator::invariant(Eigen::Matrix3d const &_mat, types::t_real _tolerance) const
    {
      Eigen::Matrix3d const mat((!_mat) * op * _mat);
      for(size_t i(0); i < 3; ++i)
        for(size_t j(0); j < 3; ++j)
          if( std::abs( std::floor(0.01+mat(i,j)) - mat(i,j) ) > types::tolerance ) 
            return false;
      return true;
    }

    // returns true if matrix is the identity.
    bool is_identity( Eigen::Matrix3d const &_cell, types::t_real _tolerance )
    {
      if( std::abs(_cell(0,0)-1e0) > _tolerance ) return false;
      if( std::abs(_cell(1,1)-1e0) > _tolerance ) return false;
      if( std::abs(_cell(2,2)-1e0) > _tolerance ) return false;
      if( std::abs(_cell(0,1)) > _tolerance ) return false;
      if( std::abs(_cell(0,2)) > _tolerance ) return false;
      if( std::abs(_cell(1,0)) > _tolerance ) return false;
      if( std::abs(_cell(1,2)) > _tolerance ) return false;
      if( std::abs(_cell(2,0)) > _tolerance ) return false;
      if( std::abs(_cell(2,1)) > _tolerance ) return false;
      return true;
    }

    //! Returns point symmetries of a cell, except identity.
    boost::shared_ptr< std::vector<SymmetryOperator> >
      get_point_group_symmetries( Eigen::Matrix3d const &_cell, types::t_real _tolerance )
      {
        if( _tolerance <= 0e0 ) _tolerance = types::tolerance;
        boost::shared_ptr< std::vector<SymmetryOperator> >
           result( new std::vector<SymmetryOperator> ); 

        
        // Finds out how far to look.
        types::t_real const volume( std::abs(_cell.determinant()) );
        Eigen::Vector3d const a0( _cell.col(0) );
        Eigen::Vector3d const a1( _cell.col(1) );
        Eigen::Vector3d const a2( _cell.col(2) );
        types::t_real const max_norm = std::max( a0.norm(), std::max(a1.norm(), a2.norm()) );
        int const n0( std::ceil(max_norm*a1.cross(a2).norm()/volume) );
        int const n1( std::ceil(max_norm*a2.cross(a0).norm()/volume) );
        int const n2( std::ceil(max_norm*a0.cross(a1).norm()/volume) );
        types::t_real const length_a0( a0.squaredNorm() );
        types::t_real const length_a1( a1.squaredNorm() );
        types::t_real const length_a2( a2.squaredNorm() );

        // now creates a vector of all G-vectors in the sphere of radius max_norm. 
        std::vector< Eigen::Vector3d > gvectors0;
        std::vector< Eigen::Vector3d > gvectors1;
        std::vector< Eigen::Vector3d > gvectors2;
        for( int i0(-n0); i0 <= n0; ++i0 )
          for( int i1(-n1); i1 <= n1; ++i1 )
            for( int i2(-n2); i2 <= n2; ++i2 )
            {
              Eigen::Vector3d const g = _cell * Eigen::Vector3d(i0, i1, i2);
              types::t_real length( g.squaredNorm() );
              if( std::abs(length-length_a0) < _tolerance ) gvectors0.push_back(g); 
              if( std::abs(length-length_a1) < _tolerance ) gvectors1.push_back(g); 
              if( std::abs(length-length_a2) < _tolerance ) gvectors2.push_back(g); 
            }


        // Adds triplets which are rotations.
        typedef std::vector<Eigen::Vector3d> :: const_iterator t_cit;
        Eigen::Matrix3d const inv_cell( !_cell );
        foreach( Eigen::Vector3d const & rot_a0, gvectors0 )
          foreach( Eigen::Vector3d const & rot_a1, gvectors1 )
            foreach( Eigen::Vector3d const & rot_a2, gvectors2 )
            {
              // creates matrix.
              Eigen::Matrix3d rotation;
              rotation.set_column(0, rot_a0);
              rotation.set_column(1, rot_a1);
              rotation.set_column(2, rot_a2);

              // checks that this the rotation is not singular.
              if( std::abs(rotation.determinant()) < _tolerance ) continue;

              rotation = rotation * inv_cell;
              // avoids identity.
              if( is_identity(rotation, _tolerance) ) continue;
              // checks that the rotation is a rotation.
              if( not is_identity(rotation * (~rotation), _tolerance) ) continue;

              // adds to vector of symmetries.
              SymmetryOperator symop(rotation);
              if( result->end() == std::find( result->begin(), result->end(), symop) )
                result->push_back( symop );
            }
        return result;
      }


    boost::shared_ptr< std::vector<SymmetryOperator> >
      get_space_group_symmetries( Lattice const &_lattice, types::t_real _tolerance )
      {
        if( _tolerance <= 0e0 ) _tolerance = types::tolerance;
        // Checks that lattice has sites.
        if( _lattice.sites.size() == 0 )
        {
          std::cerr << "Lattice does not contain sites.\n"
                       "Will not compute symmetries for empty lattice.\n";
          return boost::shared_ptr< std::vector<SymmetryOperator> >();
        }
        { // Checks that lattice is primitive.
          Lattice lat(_lattice);
          if( not lat.make_primitive() )
          {
            std::cerr << "Lattice is not primitive.\nCannot compute symmetries.\n";
            return boost::shared_ptr< std::vector<SymmetryOperator> >();
          }
        }

        // Finds minimum translation.
        Lattice::t_Sites sites(_lattice.sites);
        Eigen::Vector3d translation(sites.front().pos);
        Eigen::Matrix3d const invcell(!_lattice.cell);
        foreach( Lattice::t_Site &site, sites )
          if( translation.squaredNorm() > site.pos.squaredNorm() )
            translation = site.pos;
        // Creates a list of sites centered in the cell.
        foreach( Lattice::t_Site &site, sites )
          site.pos = into_cell(site.pos-translation, _lattice.cell, invcell);

        // gets point group.
        boost::shared_ptr< std::vector<SymmetryOperator> >
          pg = get_point_group_symmetries(_lattice.cell);
        boost::shared_ptr< std::vector<SymmetryOperator> > 
          result( new std::vector<SymmetryOperator> );
        result->reserve(pg->size());
             
        // lists sites of same type as sites.front()
        std::vector<Eigen::Vector3d> translations;
        CompareSites compsites(sites.front(), _tolerance);
        foreach( Lattice::t_Site const &site, sites )
          if( compsites(site.type) ) translations.push_back(site.pos);
        

        // applies point group symmetries and finds out if they are part of the space-group.
        foreach(SymmetryOperator &op, *pg)
        {
          // loop over possible translations.
          std::vector<Eigen::Vector3d> :: const_iterator i_trial = translations.begin();
          std::vector<Eigen::Vector3d> :: const_iterator const i_trial_end = translations.end();
          for(; i_trial != i_trial_end; ++i_trial)
          {
            // possible translation.
            op.trans = *i_trial;
            Lattice::t_Sites::const_iterator i_site = sites.begin();
            Lattice::t_Sites::const_iterator const i_site_end = sites.end();
            for(; i_site != i_site_end; ++i_site)
            {
              CompareSites transformed(*i_site, _tolerance);
              transformed.pos = into_cell(op(i_site->pos), _lattice.cell, invcell);
              Lattice::t_Sites::const_iterator const
                i_found( std::find_if(sites.begin(), sites.end(), transformed) );
              if(i_found == sites.end()) break;
              if( not transformed(i_found->type) ) break;
            } // loop over all sites.

            if(i_site == i_site_end) break; // found a mapping
          } // loop over trial translations.

          if(i_trial != i_trial_end)
            result->push_back(SymmetryOperator(op.op, op.trans-op.op*translation+translation) );
        // else
        // {
        //   std::cout << op.op << "\n";
        //   std::vector<Eigen::Vector3d> :: const_iterator i_trial = translations.begin();
        //   std::vector<Eigen::Vector3d> :: const_iterator const i_trial_end = translations.end();
        //   for(; i_trial != i_trial_end; ++i_trial)
        //     foreach( Lattice::t_Site const & site, sites )
        //     {
        //       Eigen::Vector3d pos = into_cell(op.op * sites.front().pos, _lattice.cell, invcell);
        //       pos = into_cell( *i_trial - pos, _lattice.cell, invcell );
        //       pos =op.op * site.pos + pos;
        //       std::cout << "              " << into_cell(pos, _lattice.cell, invcell)<< "\n";
        //     }
        //   std::cout << "\n";
        // }
        } // loop over point group.

        return result;
      } 

  } // namespace Crystal

} // namespace LaDa
