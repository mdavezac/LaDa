//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/symmetry_operator.h>


namespace LaDa
{

  namespace Crystal 
  {

    void compose( SymmetryOperator const &_a, SymmetryOperator const &_b, SymmetryOperator &_out )
    {
      _out.op = _a.op * _b.op;
      _out.trans = _a.trans + _a.op * _b.trans;
    }
    
    template<class T_CONTAINER>
      void transform_impl_( atat::SpaceGroup const &_sg, std::vector<SymmetryOperator> &_symops )
      {
        for(size_t i(0); i < _sg.point_op.size(); ++i )
          _symops.push_back(SymmetryOperator(_sg.point_op[i], _sg.trans[i]));
      }
    void transform( atat::SpaceGroup const &_sg, std::vector<SymmetryOperator> &_symops )
      { transform_impl_(_sg, _symops); }

  } // namespace Crystal

} // namespace LaDa
