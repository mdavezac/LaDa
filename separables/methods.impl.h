#include <boost/numeric/ublas/io.hpp>
namespace LaDa
{ 
  namespace CE
  {
    namespace Method
    {
      namespace Policy
      {   
        template< class T_SAVEDOBJECT > template< class T_COLLAPSE, class T_MINIMIZER >
          bool BestOf<T_SAVEDOBJECT> :: go( T_COLLAPSE &_col,
                                            T_MINIMIZER &_min, 
                                            types::t_int _verb )
          {
            LADA_NASSERT( restarts == 0, "No runs required.\n" )
            _col.randomize( howrandom ); 
            if( restarts == 1 ) return false;  
            LADA_NASSERT( which > 2, "Not sure what to check.\n" )
            
            opt::ErrorTuple intermed = _min( _col.coefficients(), _col );
            if( _verb >= 1 ) std::cout << "  trial " << nbrestarts
                                       << " " << intermed << "\n";
            if(    nbrestarts == 0 
                or ( which == 0 and math::gt( best.variance(), intermed.variance() )  )
                or ( which == 1 and math::gt( best.mean(), intermed.mean() )  )
                or ( which == 2 and math::gt( best.max(), intermed.max() )  ) )
            {
              best = intermed;
              state = _col;
            }
            ++nbrestarts;
            return nbrestarts < restarts;
          }
        template< class T_SAVEDOBJECT > template< class T_COLLAPSE, class T_MINIMIZER >
          opt::ErrorTuple BestOf<T_SAVEDOBJECT> :: end( T_COLLAPSE &_col,
                                                        T_MINIMIZER &_min,
                                                        types::t_int _verb ) const
          {
            if( restarts == 1 ) return _min( _col.coefficients(), _col );
            state.reset(_col);
            return best; 
          }
      } // end of Policy namespace.

      template< class T_COLLAPSE, class T_FIT, class T_MINIMIZER >
        opt::t_ErrorPair leave_one_out( T_COLLAPSE &_collapse,
                                        T_FIT &_fit,
                                        const T_MINIMIZER &_min,
                                        types::t_int _verbosity )
        {
          LADA_TRY_BEGIN
          opt::t_ErrorPair errors;
          for( _collapse.mapping().n = 0;
               _collapse.mapping().n < _collapse.mapping().size();
               ++_collapse.mapping().n )
          {
            opt::ErrorTuple intermediate;
            if( _verbosity >= 1 ) std::cout << " " << _collapse.mapping().n
                                            << ". Training Errors: ";
            if( _verbosity >= 2 ) std::cout << "\n";
            intermediate = _fit( _collapse, _min );
            if( _verbosity >= 1 ) std::cout << intermediate << "\n";
            errors.first += intermediate;

            if( _verbosity >= 1 ) std::cout << " " << _collapse.mapping().n
                                            << ". Prediction Errors: ";
            if( _verbosity >= 2 ) std::cout << "\n";
            const Crystal::Structure& structure = _fit.structures()[ _collapse.mapping().n];
            intermediate = check_one( _collapse, structure,
                                      _collapse.mapping().n, _verbosity >= 2 );
            if( _verbosity >= 1 ) std::cout << intermediate << "\n";
            errors.second += intermediate;
          }
          return errors;
          LADA_TRY_END(,"Error in CE::Methods::leave_one_out().\n" )
        }

      template< class T_COLLAPSE, class T_FIT, class T_MINIMIZER >
        opt::t_ErrorPair leave_many_out( Fitting::LeaveManyOut &_lmo,
                                         T_COLLAPSE &_collapse,
                                         T_FIT &_fit,
                                         const T_MINIMIZER &_min )
      {
        LADA_TRY_BEGIN
        opt::t_ErrorPair errors;
        if( not _lmo.do_perform ) return errors;
        if( not _lmo.sets.size() ) _lmo.create_sets( _collapse.mapping().size() );
   
        typedef std::vector< std::vector< types::t_unsigned > >
                                      :: const_iterator const_iterator;
        typedef typename T_COLLAPSE :: t_Traits t_CollapseTraits;
        typedef typename t_CollapseTraits :: t_Mapping t_Mapping;
        typedef typename t_Mapping :: t_Container t_Container;
        t_Container& excluded = _collapse.mapping().excluded;
        const_iterator i_set = _lmo.sets.begin();
        const_iterator i_set_end = _lmo.sets.end();
        for(size_t n(0); i_set != i_set_end; ++i_set, ++n )
        {
          { // Fitting
            opt::ErrorTuple intermediate;
            excluded.resize( i_set->size() );
            std::copy( i_set->begin(), i_set->end(), excluded.begin() );
            
            if( _lmo.verbosity >= 1 ) std::cout << " " << n
                                            << ". Training Errors: ";
            if( _lmo.verbosity >= 2 ) std::cout << "\n";
            intermediate = _fit( _collapse, _min );
            if( _lmo.verbosity >= 1 ) std::cout << intermediate << "\n";
            errors.first += intermediate;
          }

          { // Prediction
            opt::ErrorTuple intermediate;
            if( _lmo.verbosity >= 1 ) std::cout << " " << n
                                            << ". Prediction Errors: ";
            if( _lmo.verbosity >= 2 ) std::cout << "\n";
            for( size_t i(0); i < _collapse.mapping().size(); ++i )
            {
              if( not _collapse.mapping().do_skip( i ) ) continue;
              const Crystal::Structure& structure = _fit.structures()[ i ];
              intermediate += check_one( _collapse, structure, i, _lmo.verbosity >= 2 );
            }
            if( _lmo.verbosity >= 1 ) std::cout << intermediate << "\n";
            errors.second += intermediate;
          }
        }
   
        return errors;
        LADA_TRY_END(, "Error while performing leave-many-out.\n" )
      }

      template< class T_COLLAPSE >
        opt::ErrorTuple check_one( const T_COLLAPSE &_collapse,
                                   const Crystal::Structure &_structure,
                                   size_t _n, bool _verbose )
        {
          LADA_DEBUG_TRY_BEGIN
          namespace fs = boost::filesystem;
          namespace bblas = boost::numeric::ublas;
          const types::t_real weight = _structure.weight;
          const std::string name = fs::path( _structure.name ).leaf() ;
       
          typedef typename T_COLLAPSE::t_Matrix::value_type t_Type; 
          t_Type predic( _collapse.evaluate(_n) );
          opt::ErrorTuple error( _structure.energy - predic, _structure.weight );
          if( _verbose )
            std::cout << "  structure: " << std::setw(30) << name << "  "
                      << "  x=" << std::setw(5) << _structure.get_concentration() << "  "
                      << "Target: " << std::fixed << std::setw(8) 
                      << std::setprecision(2) << _structure.energy << " "
                      << "Separable: " << std::fixed << std::setw(8)
                      << std::setprecision(2) << predic << "   "
                      << "|Target-Separable| * weight: "
                      << std::fixed << std::setw(10) << std::setprecision(3) 
                      << error.mean()
                      << "\n";
          return error;
          LADA_DEBUG_TRY_END(, "Error in Fit::check_one().\n" )
        }

      template< class T_COLLAPSE, class T_STRUCTURES >
        opt::ErrorTuple check_all( const T_COLLAPSE &_collapse,
                                   const T_STRUCTURES &_strs,
                                   bool _verbose )
        {
          opt::ErrorTuple result;
          typename T_STRUCTURES :: const_iterator i_str = _strs.begin();
          typename T_STRUCTURES :: const_iterator i_str_end = _strs.end();
          for(size_t n(0); i_str != i_str_end; ++i_str, ++n )
          {
            if( _collapse.mapping().do_skip(n) ) continue;
            result += check_one( _collapse, *i_str, n, _verbose );
          }
          return result;
        }
    } // end of namespace Methods
  } // end of CE namespace.

  namespace Fitting
  {
    template< class T_SOLVER > template< class T_COLLAPSE >
      opt::ErrorTuple AlternatingLeastSquare<T_SOLVER> 
        :: operator()( typename T_COLLAPSE :: t_Matrix &_solution,
                       T_COLLAPSE &_collapse  ) const
        {
          LADA_TRY_BEGIN
          namespace bblas = boost::numeric::ublas;
          typedef typename T_COLLAPSE::t_Matrix t_Matrix;
          typedef typename bblas::matrix< typename t_Matrix::value_type> t_OMatrix;
          types::t_real convergence( 1e1 * tolerance );
          if( verbose ) std::cout << "Starting Alternating-least-square fit.\n";
          types::t_unsigned iter = 0;
          const size_t D( _solution.size2() );
          t_OMatrix A( _collapse.dof(), _collapse.dof() );
          t_Vector b( _collapse.dof() );
          _collapse.update_all();
          opt::ErrorTuple errors( _collapse.evaluate() );
          if( verbose ) std::cout << "Allsq start: " << errors << "\n";
          do
          {
            for(size_t dim(0); dim < D; ++dim )
            {
              typedef bblas::matrix_column<t_Matrix> t_Column;
              t_Column column( _solution, dim );
              _collapse( A, b, dim );
              if( A.size1() == column.size() ) 
                linear_solver( A, column, b );
              else
              {
                bblas::vector_range< t_Column >
                  range( column, bblas::range( 0, A.size1() ) );
                linear_solver( A, range, b );
              }
              _collapse.update( dim );
              opt::ErrorTuple e( _collapse.evaluate() );
            }
            ++iter;
         
            if( tolerance > 0e0  )
            {
              opt::ErrorTuple newerrors = _collapse.evaluate();
              convergence = newerrors.variance() - errors.variance();
              errors = newerrors;
              if( verbose )
              {
                if( iter == 1 )
                  std::cout << "\n  Allsq iter: " << iter << errors << "\n";
                else
                  std::cout << "\n  Allsq iter: " << iter << errors
                            << "  convergence: " << convergence << " \n";
              } // end of if( verbose )
              if( D == 1 ) return errors;
              if( iter > 1 and std::abs(convergence) < tolerance ) return errors; 
            } // end of if( tolerance > 0e0 )
          }
          while( iter < itermax or itermax == 0 );
         
          return errors;
          LADA_TRY_END(, "Error encountered in Alternating-least square fit.\n" )
        }

      template< class T_SOLVER > 
        void AlternatingLeastSquare<T_SOLVER> :: Load( const TiXmlElement &_node )
        {
          LADA_TRY_BEGIN
          std::string name = _node.Value();
          const TiXmlElement *parent = &_node;
          if( name.compare( "Allsq" ) ) 
           parent = _node.FirstChildElement( "Allsq" );
          LADA_DO_NASSERT( not parent, "Could not find Allsq tag in input.\n" )
          if( parent->Attribute( "tolerance" ) )
            parent->Attribute( "tolerance", &tolerance );
          if( parent->Attribute( "itermax" ) )
            parent->Attribute( "itermax", &itermax );
          LADA_TRY_END(,"Error in AlternatingLeastSquare::Load()\n" )
        }
  } // end of Fitting namespace
} // namespace LaDa
