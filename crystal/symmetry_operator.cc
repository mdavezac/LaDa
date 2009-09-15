//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <crystal/lattice.h>
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
    
    boost::shared_ptr< std::vector<SymmetryOperator> > transform(atat::SpaceGroup const &_sg) 
    {
      boost::shared_ptr< std::vector<SymmetryOperator> >
        result( new std::vector<SymmetryOperator> );
      for(size_t i(0); i < _sg.point_op.get_size(); ++i )
      {
        atat::rMatrix3d const &op( _sg.point_op[i] );
        atat::rVector3d const &trans( _sg.trans[i] );
        if(    Fuzzy::neq(op.x[0][0], 1.0)
            or Fuzzy::neq(op.x[0][1], 0.0) 
            or Fuzzy::neq(op.x[0][2], 0.0)
            or Fuzzy::neq(op.x[1][0], 0.0) 
            or Fuzzy::neq(op.x[1][1], 1.0) 
            or Fuzzy::neq(op.x[1][2], 0.0)
            or Fuzzy::neq(op.x[2][0], 0.0) 
            or Fuzzy::neq(op.x[2][1], 0.0) 
            or Fuzzy::neq(op.x[2][2], 1.0)  
            or Fuzzy::neq(trans.x[0], 0.0) 
            or Fuzzy::neq(trans.x[1], 0.0) 
            or Fuzzy::neq(trans.x[2], 0.0) )
          result->push_back(op);
      }
      return result;
    }

    // boost::shared_ptr< std::vector<SymmetryOperator> > get_symmetries_multisite( Lattice const &_lat );

    boost::shared_ptr< std::vector<SymmetryOperator> > get_symmetries( Lattice const &_lat )
    {
      if( _lat.sites.size() == 0 ) return boost::shared_ptr< std::vector<SymmetryOperator> >();
      Lattice lattice;
      lattice.cell = _lat.cell;
      lattice.sites = _lat.sites;
      lattice.find_space_group();
      return transform( lattice.space_group );
    }

//   boost::shared_ptr< std::vector<SymmetryOperator> > get_symmetries_multisite( Lattice const &_lat )
//   {
//     // Creates a list of sites centered in the cell.
//     Lattice::t_Sites sites(_lat.sites);
//     atat::rMatrix3d center(0,0,0);
//     atat::rMatrix3d const invcell(!_lat.cell);
//     foreach( Lattice::t_Site &site, sites )
//     {
//       atat::rVector3d const inved( invcell * site.pos );
//       atat::rVector3d const centered
//       ( 
//         invcell * site.pos 
//       );
//     }
//
//     // Gets all possible rotations.
//     atat::Array<atat::rMatrix3d> array;
//     atat::find_pointgroup(array, _lat.cell);
//     for(size_t i(0); i < array.get_size(); ++i )
//     {
//       atat::rMatrix3d const &op( array[i] );
//       // avoid identity
//       if( not(    Fuzzy::neq(op.x[0][0], 1.0)
//                or Fuzzy::neq(op.x[0][1], 0.0)
//                or Fuzzy::neq(op.x[0][2], 0.0)
//                or Fuzzy::neq(op.x[1][0], 0.0)
//                or Fuzzy::neq(op.x[1][1], 1.0)
//                or Fuzzy::neq(op.x[1][2], 0.0)
//                or Fuzzy::neq(op.x[2][0], 0.0) 
//                or Fuzzy::neq(op.x[2][1], 0.0) 
//                or Fuzzy::neq(op.x[2][2], 1.0)  
//                or Fuzzy::neq(trans.x[0], 0.0) 
//                or Fuzzy::neq(trans.x[1], 0.0) 
//                or Fuzzy::neq(trans.x[2], 0.0) ) ) continue; 
//
//       // now finds possible translation.
//     }
//   }

  } // namespace Crystal

} // namespace LaDa
