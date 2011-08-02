#include "LaDaConfig.h"

#include "compare_sites.h"
#include "symmetry_operator.h"


namespace LaDa
{

  namespace crystal 
  {

    void compose( SymmetryOperator const &_a, SymmetryOperator const &_b, SymmetryOperator &_out )
    {
      _out.op = _a.op * _b.op;
      _out.trans = _a.trans + _a.op * _b.trans;
    }
    
    bool SymmetryOperator::invariant(math::rMatrix3d const &_mat, types::t_real _tolerance) const
    {
      math::rMatrix3d const mat(_mat.inverse() * op * _mat);
      for(size_t i(0); i < 3; ++i)
        for(size_t j(0); j < 3; ++j)
          if( std::abs( std::floor(types::roundoff+mat(i,j)) - mat(i,j) ) > types::tolerance ) 
            return false;
      return true;
    }


  } // namespace Crystal
} // namespace LaDa
