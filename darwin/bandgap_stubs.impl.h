//
//  Version: $Id$
//
#ifndef _DARWIN_BANDGAP_STUBS_IMPL_H_
#define _DARWIN_BANDGAP_STUBS_IMPL_H_

namespace BandGap
{
  inline std::ostream& operator<<(std::ostream &_stream, const Keeper &_o)
  { 
    _stream << " CBM " << std::fixed << std::setw(12) << std::setprecision(6) << _o.cbm 
            << "  --  VBM " << std::fixed << std::setw(12) << std::setprecision(6) << _o.vbm; 
    return _stream; 
  } 

  inline void Darwin :: operator()( Keeper &_keeper )
  {
    Darwin::operator()();
    // copies band edges into object
    _keeper.vbm = bandgap.bands.vbm; 
    _keeper.cbm = bandgap.bands.cbm;
  }
  template <class T_BASE> 
  inline void Darwin::operator<<( const Vff::Darwin<T_BASE> &_vff )
  {
    // creates an mpi aware file name for atomic configurations
    std::ostringstream  sstr;
    sstr << "atom_config";
    __MPICODE( sstr << "." << mpi::main.rank(); )
    // prints atomic configurations
    _vff.print_escan_input(sstr.str());
    // tells bandgap where to find atomic configurations
    atomicconfig = sstr.str();
  }
} // namespace bandgap

#endif
