#ifndef LADA_CRYSTAL_SYMMETRY_OPERATOR_H_
#define LADA_CRYSTAL_SYMMETRY_OPERATOR_H_

#include "LaDaConfig.h"

#include <boost/serialization/serialization.hpp>


#include <opt/debug.h>
#include <opt/types.h>

#include "fuzzy.h"
#include "misc.h"
#include "set_cell.h"

namespace Eigen
{
}

namespace LaDa
{
  namespace math
  {
    //! True if invariant by this transform. Ignores translation.
    bool inline invariant(Affine3d const&_a, rMatrix3d const &_m) { return eq(_m, _a.linear() * _m); }
    //! True if invariant by this transform. Ignores translation.
    bool inline invariant(Affine3d const&_a, rMatrix3d const &_m, types::t_real const &_tol) 
      { return eq(_m, _a.linear() * _m, _tol); }
    //! True if invariant by this transform.
    bool inline invariant(Affine3d const&_a, rVector3d const &_m) { return eq(_m, _a * _m); }
    //! True if invariant by this transform.
    bool inline invariant(Affine3d const&_a, rVector3d const &_m, types::t_real const &_tol) 
      { return eq(_m, _a * _m, _tol); }
    //! Checks if this is a pure rotation. 
    bool is_rotation(Affine3d const &_a) { return is_null(_a.translation); }

    //! A symmetry operator.
    struct SymmetryOperator
    {
      //! Constructor.
      SymmetryOperator() : trans(math::rVector3d::Zero()), op(math::rMatrix3d::Identity()) {}
      //! Constructor.
      SymmetryOperator   ( math::rMatrix3d const &_op,
                           math::rVector3d const& _t = math::rVector3d(0,0,0) )
                       : op(_op), trans(_t) {}
      //! Constructor.
      SymmetryOperator   ( math::rVector3d const& _t)
                       : trans(_t), op(math::rMatrix3d::Identity()) {}
      //! Copy Constructor.
      SymmetryOperator( SymmetryOperator const& _c ) : op(_c.op), trans(_c.trans) {}

      //! Matrix operator.
      math::rMatrix3d op;
      //! Vector translation
      math::rVector3d trans;
      //! Applies operator.
      math::rVector3d operator()( math::rVector3d const &_a ) const
        { return op*_a + trans; }

      //! True if the matrix is invariant by this operator.
      bool invariant(math::rMatrix3d const &_mat, types::t_real _tolerance = types::tolerance) const;
      //! Comparison.
      bool operator==(SymmetryOperator const &_sym) const
      {
        return not(    math::neq(op(0,0), _sym.op(0,0))
                    or math::neq(op(1,0), _sym.op(1,0))
                    or math::neq(op(2,0), _sym.op(2,0))
                    or math::neq(op(0,1), _sym.op(0,1))
                    or math::neq(op(1,1), _sym.op(1,1)) 
                    or math::neq(op(2,1), _sym.op(2,1))
                    or math::neq(op(0,2), _sym.op(0,2)) 
                    or math::neq(op(1,2), _sym.op(1,2))
                    or math::neq(op(2,2), _sym.op(2,2)) 
                    or math::neq(trans(0), _sym.trans(0))
                    or math::neq(trans(1), _sym.trans(1)) 
                    or math::neq(trans(2), _sym.trans(2)) );
      }
      //! Serializes a symmetry operator.
      template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version)
        { _ar & op; _ar & trans; }

      //! Returns the inverse operator.
      SymmetryOperator inverse() const { return SymmetryOperator(~op, -(~op)*trans); }

      //! Initializer for rotation.
      math::details::SetCell< boost::mpl::int_<1> >
        set_rotation(types::t_real _x, types::t_real _y, types::t_real _z)
          { return math::details::SetCell< boost::mpl::int_<0> >(op)(_x, _y, _z); }
      //! Initializer for rotation.
      math::details::SetCell< boost::mpl::int_<1> >
        set_rotation(math::rVector3d _pos)
          { return math::details::SetCell< boost::mpl::int_<0> >(op)(_pos); }

      //! Initializer for the translation part of the operator.
      void set_translation(types::t_real _x, types::t_real _y, types::t_real _z)
        { trans(0) = _x; trans(1) = _y; trans(2) = _z; }
      //! Initializer for the translation part of the operator.
      void set_translation(math::rVector3d const &_trans) { trans = _trans; }

      //! Checks that this is indeed a symmetry operation.
      bool is_symmetry() const;

      //! Checks if this is a pure rotation. 
      bool is_rotation() const;
    };

    inline std::ostream& operator<<( std::ostream& _stream, SymmetryOperator const &_sym )
    {
      return _stream << "Trans: " << _sym.trans.transpose() << "\n" << _sym.op << "\n"; 
    }
    //! \brief Composes two symmetry operations.
    //! \details \a _out.op = \a _a.op * _\a b.op, \a _out.trans = _a.trans + _a.op * _b.trans.
    void compose( SymmetryOperator const &_a,
                  SymmetryOperator const &_b, 
                  SymmetryOperator &_out );

  }
}

#endif
