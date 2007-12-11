//
//  Version: $Id$
//
#ifndef _DARWIN_PESCAN_IMPL_H_
#define _DARWIN_PESCAN_IMPL_H_


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
#ifdef _MPI
    sstr << "." << mpi::main.rank();
#endif
    // prints atomic configurations
    _vff.print_escan_input(sstr.str());
    // tells bandgap where to find atomic configurations
    atomicconfig = sstr.str();
  }
} // namespace bandgap

namespace eMassSL
{
  inline std::ostream& operator<<(std::ostream &_stream, const Keeper &_o)
  { 
    _stream << " CBM " << std::fixed << std::setw(12) << std::setprecision(6) << _o.cbm 
            << "  --  eMass " 
            << std::fixed << std::setw(8) << std::setprecision(4) << emass(0,0) << " " 
            << std::fixed << std::setw(8) << std::setprecision(4) << emass(0,1) << " " 
            << std::fixed << std::setw(8) << std::setprecision(4) << emass(0,2) << " "
            << std::fixed << std::setw(8) << std::setprecision(4) << emass(1,1) << " " 
            << std::fixed << std::setw(8) << std::setprecision(4) << emass(1,2) << " "
            << std::fixed << std::setw(8) << std::setprecision(4) << emass(2,2);

    return _stream; 
  } 

  inline void Darwin :: operator()( Keeper &_keeper )
  {
    Darwin::operator()();
    // copies band edges into object
    _keeper.cbm = emass.eigenvalue.bacl();
    _keeper.emass = emass.tensor();
  }
  template <class T_BASE> 
  inline void Darwin::operator<<( const Vff::Darwin<T_BASE> &_vff )
  {
    // creates an mpi aware file name for atomic configurations
    std::ostringstream  sstr;
    sstr << "atom_config";
#ifdef _MPI
    sstr << "." << mpi::main.rank();
#endif
    // prints atomic configurations
    _vff.print_escan_input(sstr.str());
    // tells bandgap where to find atomic configurations
    atomicconfig = sstr.str();
  }
}

#endif
