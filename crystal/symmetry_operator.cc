#include "LaDaConfig.h"

#include "compare_sites.h"
#include "symmetry_operator.h"


namespace LaDa
{

  namespace Crystal 
  {
    //! Returns point symmetries of a cell, except identity.
    boost::shared_ptr<t_SpaceGroup>
      get_cell_symmetries( math::rMatrix3d const &_cell, types::t_real _tolerance )
      {
        if( _tolerance <= 0e0 ) _tolerance = types::tolerance;
        boost::shared_ptr<t_SpaceGroup> result(new t_SpaceGroup); 
        
        // Finds out how far to look.
        types::t_real const volume( std::abs(_cell.determinant()) );
        math::rVector3d const a0( _cell.col(0) );
        math::rVector3d const a1( _cell.col(1) );
        math::rVector3d const a2( _cell.col(2) );
        types::t_real const max_norm = std::max( a0.norm(), std::max(a1.norm(), a2.norm()) );
        int const n0( std::ceil(max_norm*(a1^a2).norm()/volume) );
        int const n1( std::ceil(max_norm*(a2^a0).norm()/volume) );
        int const n2( std::ceil(max_norm*(a0^a1).norm()/volume) );
        types::t_real const length_a0( a0.squaredNorm() );
        types::t_real const length_a1( a1.squaredNorm() );
        types::t_real const length_a2( a2.squaredNorm() );

        // now creates a vector of all G-vectors in the sphere of radius max_norm. 
        std::vector< math::rVector3d > gvectors0;
        std::vector< math::rVector3d > gvectors1;
        std::vector< math::rVector3d > gvectors2;
        for( int i0(-n0); i0 <= n0; ++i0 )
          for( int i1(-n1); i1 <= n1; ++i1 )
            for( int i2(-n2); i2 <= n2; ++i2 )
            {
              math::rVector3d const g = _cell * math::rVector3d(i0, i1, i2);
              types::t_real length( g.squaredNorm() );
              if( std::abs(length-length_a0) < _tolerance ) gvectors0.push_back(g); 
              if( std::abs(length-length_a1) < _tolerance ) gvectors1.push_back(g); 
              if( std::abs(length-length_a2) < _tolerance ) gvectors2.push_back(g); 
            }


        // Adds triplets which are rotations.
        typedef std::vector<math::rVector3d> :: const_iterator t_cit;
        math::rMatrix3d const inv_cell(_cell.inverse());
        foreach( math::rVector3d const & rot_a0, gvectors0 )
          foreach( math::rVector3d const & rot_a1, gvectors1 )
            foreach( math::rVector3d const & rot_a2, gvectors2 )
            {
              // creates matrix.
              math::rMatrix3d rotation;
              rotation.col(0) = rot_a0;
              rotation.col(1) = rot_a1;
              rotation.col(2) = rot_a2;

              // checks that this the rotation is not singular.
              if( math::is_null(rotation.determinant(), _tolerance) ) continue;

              rotation = rotation * inv_cell;
              // avoids identity.
              if( math::is_identity(rotation, _tolerance) ) continue;
              // checks that the rotation is a rotation.
              if( not math::is_identity(rotation * (~rotation), _tolerance) ) continue;

              // adds to vector of symmetries.
              SymmetryOperator symop; symop.linear() = rotation;
              if( result->end() == std::find( result->begin(), result->end(), symop) )
                result->push_back( symop );
            }
        return result;
      }



  } // namespace Crystal

} // namespace LaDa
