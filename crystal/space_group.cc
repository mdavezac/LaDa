#include "LaDaConfig.h"

#include <algorithm>

#include "space_group.h"


namespace LaDa
{

  namespace crystal 
  {
    struct Cmp
    {
      math::Affine3d const &aff;
      types::t_real tolerance;
      Cmp(math::Affine3d const &_aff, types::t_real _tol) : aff(_aff), tolerance(_tol) {}
      Cmp(Cmp const &_c) : aff(_c.aff), tolerance(_c.tolerance) {}
      bool operator()(math::Affine3d const &_a) const 
        { return math::eq(_a.linear(), aff.linear(), tolerance); }
    };
    //! Returns point symmetries of a cell, with identity as first element.
    boost::shared_ptr<t_SpaceGroup>
      cell_invariants( math::rMatrix3d const &_cell, types::t_real _tolerance )
      {
        if( _tolerance <= 0e0 ) _tolerance = types::tolerance;
        boost::shared_ptr<t_SpaceGroup> result(new t_SpaceGroup); 
        result->reserve(48);
        result->push_back(math::Affine3d(math::AngleAxis(0, math::rVector3d::UnitX())));
        
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
        typedef std::vector<math::rVector3d, Eigen::aligned_allocator<math::rVector3d> > t_vector;
        t_vector gvectors[3];
        for( int i0(-n0); i0 <= n0; ++i0 )
          for( int i1(-n1); i1 <= n1; ++i1 )
            for( int i2(-n2); i2 <= n2; ++i2 )
            {
              math::rVector3d const g = _cell * math::rVector3d(i0, i1, i2);
              types::t_real length( g.squaredNorm() );
              if( std::abs(length-length_a0) < _tolerance ) gvectors[0].push_back(g); 
              if( std::abs(length-length_a1) < _tolerance ) gvectors[1].push_back(g); 
              if( std::abs(length-length_a2) < _tolerance ) gvectors[2].push_back(g); 
            }


        // Adds triplets which are rotations.
        math::rMatrix3d const inv_cell(_cell.inverse());
        t_vector::const_iterator i_a0 = gvectors[0].begin();
        t_vector::const_iterator const i_a0_end = gvectors[0].end();
        t_vector::const_iterator const i_a1_end = gvectors[1].end();
        t_vector::const_iterator const i_a2_end = gvectors[2].end();
        for(; i_a0 != i_a0_end; ++i_a0)
        {
          t_vector::const_iterator i_a1 = gvectors[1].begin();
          for(; i_a1 != i_a1_end; ++i_a1)
          {
            t_vector::const_iterator i_a2 = gvectors[2].begin();
            for(; i_a2 != i_a2_end; ++i_a2)
            {
              // creates matrix.
              math::rMatrix3d rotation;
              rotation.col(0) = *i_a0;
              rotation.col(1) = *i_a1;
              rotation.col(2) = *i_a2;

              // checks that this the rotation is not singular.
              if( math::is_null(rotation.determinant(), _tolerance) ) continue;

              rotation = rotation * inv_cell;
              // avoids identity.
              if( math::is_identity(rotation, _tolerance) ) continue;
              // checks that the rotation is a rotation.
              if( not math::is_identity(rotation * (~rotation), _tolerance) ) continue;

              // adds to vector of symmetries.
              math::Affine3d symop; symop.linear() = rotation;
              if( result->end() == std::find_if(result->begin(), result->end(), Cmp(symop, _tolerance)) )
                result->push_back( symop );
            }
          }
        }
        return result;
      }



  } // namespace Crystal

} // namespace LaDa
