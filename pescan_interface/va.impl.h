//
//  Version: $Id$
//

namespace LaDa
{
  namespace Pescan
  {
#   if defined( VAHEAD ) || defined(INVA) || defined(INVA2)
#     error "Macros with same names."
#   endif
#   define VAHEAD \
      VirtualAtom<T_VFF> 
#    define INVA( var ) \
       template< class T_VFF >  var VAHEAD
#    define INVA2( code1, code2 ) \
       template< class T_VFF >  code1, code2 VAHEAD
 
 
    INVA( typename VAHEAD::t_Type ) :: evaluate()
    { 
      namespace bfs = boost::filesystem;
      __DIAGA(
        vff.evaluate();
        boost::mpi::broadcast ( t_PescanBase::comm(), structure, 0 );
        const t_Path orig
        (
            LaDa::opt::InitialPath::path() / dirname
          / __DIAGASUFFIX( t_Path(vff.filename.filename()) )
        );
        vff.zero_order( orig );
        set_atom_input( orig );
      )
      __IIAGA( 
        __ROOTCODE( (*::mpi::main), vff.evaluate(); )
        vff.zero_order( vff.filename );
        set_atom_input( vff.filename );
      )
      t_PescanBase::escan.read_in.clear();
      t_PescanBase::escan.wavefunction_out = "zero_order";
      
      if( not t_PescanBase::operator()( structure ) ) return false; 
 
      t_PescanBase::escan.read_in.reserve( t_PescanBase::escan.nbstates );
      typedef std::vector<types::t_unsigned> :: iterator t_iterator;
      t_iterator i_r = t_PescanBase::escan.read_in.begin();
      t_iterator i_r_end = t_PescanBase::escan.read_in.end();
      for( types::t_unsigned u=1; i_r != i_r_end; ++i_r, ++u ) *i_r = u;
 
      return true;
    }
 
    INVA( typename VAHEAD::t_Type ) :: evaluate_with_gradient( t_Type* _grad )
    { 
      types :: t_real result = t_PescanBase::operator()( structure ); 
      evaluate_gradient( _grad );
      return result;
    }
 
    INVA( void ) :: evaluate_gradient( t_Type* _grad )
    {
      t_Container :: iterator i_var = t_VABase::va_vars.begin();
      t_Container :: iterator i_var_end = t_VABase::va_vars.end();
      t_Type* i_grad = _grad;
      for(types::t_unsigned n=0; i_var != i_var_end; ++i_var, ++i_grad, ++n )
        *i_grad += evaluate_one_gradient( n );
    }
 
    INVA( typename VAHEAD::t_Type ) :: apply_wfns()
    { 
      vff.evaluate();
      t_PescanBase::escan.wavefunction_in = "zero_order";
      t_PescanBase::escan.wavefunction_out = "first_order";
 
      types::t_int niter      = t_PescanBase::escan.niter;
      types::t_int nlines     = t_PescanBase::escan.nlines;
      types::t_real tolerance = t_PescanBase::escan.tolerance;
 
      t_PescanBase::escan.niter     = 0;
      t_PescanBase::escan.nlines    = 0;        
      t_PescanBase::escan.tolerance = 10000;    
      
      types::t_real result = t_PescanBase::operator()( structure ); 
 
      t_PescanBase::escan.niter     = niter;
      t_PescanBase::escan.nlines    = nlines;
      t_PescanBase::escan.tolerance = tolerance;
 
      return result;
    }
 
    INVA( bool ) :: init( bool _redocenters )
    {
      bool ret =     vff.init( _redocenters or is_vff_uninitialized ) 
                 and t_VABase::init(); 
      is_vff_uninitialized = false;
      return ret;
    }
 
    INVA( typename VAHEAD::t_Type ) :: evaluate_one_gradient( types::t_unsigned _pos )
    {
      //! finds index of atom to change
      types::t_unsigned index = 0;
      t_Atoms :: const_iterator i_atom = structure.atoms.begin();
      t_Atoms :: const_iterator i_atom_end = structure.atoms.end();
      for(++_pos; i_atom != i_atom_end and _pos; ++i_atom, ++index )
        if( not ( i_atom->freeze & t_Atom::FREEZE_T ) ) --_pos;
 
      types::t_real result = 0;
      if( do_gradients & CHEMICAL_GRADIENT )
      {
        vff.chemical( index );
        result = apply_wfns() / t_Vff::deriv_amplitude;
      }
      if( do_gradients & STRESS_GRADIENT )
      {
        vff.stress( index );
        result += apply_wfns() / t_Vff::deriv_amplitude;
      }
 
      return result;
    }
 
#  ifdef _MPI
    INVA( void ) :: set_mpi( boost::mpi::communicator *_comm, const std::string &_s ) 
    {
      __IIAGA( return; )
      __DIAGA(
        std::ostringstream sstr;
        sstr << vff.filename 
             __IIAGA( << "." << mpi::main.rank() );
        vff.filename = sstr.str();
        t_PescanBase::set_mpi( _comm ); 
      )
    }
#  endif
 
    template< class T_VFF >
      std::ostream& operator<<( std::ostream& _stream, const VAHEAD& _va )
      {
        Crystal::Fourier( _va.structure.atoms.begin(),  _va.structure.atoms.end(),
                           _va.structure.k_vecs.begin(), _va.structure.k_vecs.end() );
        return _stream << _va.structure;
      }
 
#   undef VAHEAD
#   undef INVA1
#   undef INVA2
  }
} // namespace LaDa
