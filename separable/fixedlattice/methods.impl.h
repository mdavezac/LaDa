//
//  Version: $Id$
//

namespace CE
{
  namespace Method
  {
    template< class T_SEPARABLE, class T_COLLAPSE,
              class T_MINIMIZER, class T_STRUCTURES >
      opt::ErrorTuple fit( T_SEPARABLE &_sep, 
                           T_COLLAPSE &_collapse,
                           const T_MINIMIZER &_min,
                           const T_STRUCTURES &_strs,
                           bool _verbose )
      {
        __TRYBEGIN
        _collapse.init( _sep );
        opt::ErrorTuple errors = _min( _sep.coefficients, _collapse );
        if( _verbose ) return check_all( _sep, _collapse, _strs, _verbose );
        return errors;
        __TRYEND(,"Error in CE::Methods::fit().\n" )
      }
    template< class T_SEPARABLE, class T_COLLAPSE,
              class T_MINIMIZER, class T_STRUCTURES >
      opt::t_ErrorPair leave_one_out( T_SEPARABLE &_sep, 
                                      T_COLLAPSE &_collapse,
                                      const T_STRUCTURES &_strs,
                                      const T_MINIMIZER &_min,
                                      types::t_int _verbosity )
      {
        __TRYBEGIN
        opt::t_ErrorPair errors;
        for( _collapse.mappin.n = 0;
             _collapse.mapping.n < _collapse.mapping.size();
             ++_collapse.mapping.n )
        {
          opt::ErrorTuple intermediate;
          if( _verbosity >= 1 ) std::cout << " " << _collapse.mapping.n
                                          << ". Training Errors: ";
          intermediate = fit( _sep, _collapse, _strs, _verbosity >= 2); 
          if( _verbosity ) std::cout << intermediate << "\n";
          errors.first += intermediate;

          if( _verbosity >= 1 ) std::cout << " " << _collapse.mapping.n
                                          << ". Prediction Errors: ";
          intermediate = check_one( _sep, _collapse, _strs,
                                    _collapse.mapping.n, _verbosity >= 2 );
          if( _verbosity ) std::cout << intermediate << "\n";
          errors.second += intermediate;
        }
        return errors;
        __TRYEND(,"Error in CE::Methods::leave_one_out().\n" )
      }

    template< class T_COLLAPSE, class T_SEPARABLES, class T_STRUCTURES >
      opt::ErrorTuple check_one( const T_SEPARABLES &_separables,
                                 const T_COLLAPSE &_collapse,
                                 const T_STRUCTURES &_strs,
                                 size_t _n, bool _verbose = false )
      {
        __DEBUGTRYBEGIN
        namespace fs = boost::filesystem;
        namespace bblas = boost::numeric::ublas;
        const Crystal :: Structure structure = _strs[_n];
        const types::t_real weight = structure.weight;
        const std::string name = fs::path( structure.name ).leaf() ;
     
        types::t_real predic(0);
        typedef bblas::matrix_range<typename T_SEPARABLES::t_Vector> t_MatrixRange;
        t_MatrixRange range( _collapse.configurations_,
                             bblas::range(0,_collapse.configurations_.size()),
                                          _collapse.mapping.range(_n) );
        for( size_t i(0); i < range.size2(); ++i )
        {
          typedef bblas::matrix_column< t_MatrixRange > t_Column;
          t_Column column( range, i );
          predic += _separables( column );
        }
  
        opt::ErrorTuple error( structure.energy - predic,  weight );
        if( _verbose )
          std::cout << "  structure: " << std::setw(30) << name << "  "
                    << "Target: " << std::fixed << std::setw(8) 
                    << std::setprecision(2) << structure.energy << " "
                    << "Separable: " << std::fixed << std::setw(8)
                    << std::setprecision(2) << predic << "   "
                    << "|Target-Separable| * weight: "
                    << std::fixed << std::setw(10) << std::setprecision(3) 
                    << error.mean()
                    << "\n";
        return error;
        __DEBUGTRYEND(, "Error in Fit::check_one().\n" )
      }

    template< class T_COLLAPSE, class T_SEPARABLES, class T_STRUCTURES >
      opt::ErrorTuple check_all( const T_SEPARABLES &_separables,
                                 const T_COLLAPSE &_collapse,
                                 const T_STRUCTURES &_strs,
                                 bool _verbose = false )
      {
        opt::ErrorTuple result;
        typename T_STRUCTURES :: const_iterator i_str = _strs.begin();
        typename T_STRUCTURES :: const_iterator i_str_end = _strs.end();
        for(size_t n(0); i_str != i_str_end; ++i_str, ++n )
        {
          if( _collapse.mapping.do_skip(n) ) continue;
          result += check_one( _separables, _collapse, *i_str, n, _verbose );
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
        __TRYBEGIN
        typedef typename T_COLLAPSE::t_Matrix t_Matrix;
        types::t_real convergence( 1e1 * tolerance );
        if( verbose ) std::cout << "Starting Alternating-least-square fit.\n";
        types::t_unsigned iter = 0;
        const size_t D( _solution.size2() );
        t_Matrix A;
        t_Vector b;
        opt::ErrorTuple errors;
        _collapse.update_all();
        if( verbose ) std::cout << "Allsq start: \n";
        do
        {
          for(size_t dim(0); dim < D; ++dim )
          {
            _collapse( b, A, dim );
            linear_solver( A, *i_sol, b );
            _collapse.update( dim );
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
            if( iter > 1 and std::abs(convergence) < tolerance ) return errors; 
          } // end of if( tolerance > 0e0 )
        }
        while( iter < itermax or itermax == 0 );
       
        return errors;
        __TRYEND(, "Error encountered in Alternating-least square fit.\n" )
      }

    template< class T_SOLVER > 
      void AlternatingLeastSquare<T_SOLVER> :: Load( const TiXmlElement &_node )
      {
        std::string name = _node.Value();
        const TiXmlElement *parent = &_node;
        if( name.compare( "Allsq" ) ) 
         parent = _node.FirstChildElement( "Allsq" );
        __DOASSERT( not parent, "Could not find Allsq tag in input.\n" )
        if( parent->Attribute( "tolerance" ) )
          parent->Attribute( "tolerance", &tolerance );
        if( parent->Attribute( "itermax" ) )
          parent->Attribute( "itermax", &itermax );
      }

} // end of Fitting namespace
