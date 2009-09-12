//
//  Version: $Id$
//
#ifndef LADA_SYMMETRY_OPERATOR_H_
#define LADA_SYMMETRY_OPERATOR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <opt/debug.h>
#include <opt/types.h>
#include <atat/vectmac.h>


namespace LaDa
{
  namespace Crystal 
  {
    //! \cond
    class Lattice;
    struct SymmetryOperator;
    //! \endcond
    

    //! A symmetry operator.
    struct SymmetryOperator
    {
      //! Constructor.
      SymmetryOperator   ( atat::rMatrix3d const &_op,
                           atat::rVector3d const& _t = atat::rVector3d(0,0,0) )
                       : op(_op), trans(_t) {}
      //! Constructor.
      SymmetryOperator   ( atat::rMatrix3d const &_op,
                           atat::rVector3d const& _t = atat::rVector3d(0,0,0) )
                       : op(_op), trans(_t) {}
      //! Copy Constructor.
      SymmetryOperator( SymmetryOperator const& _c ) : op(_c), trans(_t) {}

      //! Matrix operator.
      atat::rMatrix3d op;
      //! Vector translation
      atat::rVector3d trans;
      //! Applies operator.
      atat::rVector3d operator()( atat::rVector3d const &_a ) const
        { return op*_a + trans; }
    };

    //! \brief Composes two symmetry operations.
    //! \details \a _out.op = \a _a.op * _\a b.op, \a _out.trans = _a.trans + _a.op * _b.trans.
    void compose( SymmetryOperator const &_a,
                  SymmetryOperator const &_b, 
                  SymmetryOperator &_out );

    //! Transforms space group into array of SymmetryOperators.
    void transform( atat::SpaceGroup const &_sg, std::vector<SymmetryOperator> &_symops );
  }
}

#endif
