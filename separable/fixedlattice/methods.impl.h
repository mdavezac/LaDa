//
//  Version: $Id$
//

#include <boost/numeric/ublas/io.hpp>
namespace CE
{
  namespace Method
  {
    template< class T_COLLAPSE, class T_MINIMIZER, class T_STRUCTURES >
      opt::ErrorTuple fit( T_COLLAPSE &_collapse,
                           const T_MINIMIZER &_min,
                           const T_STRUCTURES &_strs,
                           bool _verbose )
      {
        __TRYBEGIN
        opt::ErrorTuple errors = _min( _collapse.separables().coefficients(), _collapse );
        if( _verbose ) return check_all( _collapse, _strs, _verbose );
        return errors;
        __TRYEND(,"Error in CE::Methods::fit().\n" )
      }
    template< class T_COLLAPSE, class T_MINIMIZER, class T_STRUCTURES >
      opt::t_ErrorPair leave_one_out( T_COLLAPSE &_collapse,
                                      const T_MINIMIZER &_min,
                                      const T_STRUCTURES &_strs,
                                      types::t_int _verbosity )
      {
        __TRYBEGIN
        opt::t_ErrorPair errors;
        for( _collapse.mapping().n = 0;
             _collapse.mapping().n < _collapse.mapping().size();
             ++_collapse.mapping().n )
        {
          opt::ErrorTuple intermediate;
          if( _verbosity >= 1 ) std::cout << " " << _collapse.mapping().n
                                          << ". Training Errors: ";
          intermediate = fit( _collapse, _min, _strs, _verbosity >= 2); 
          if( _verbosity >= 1 ) std::cout << intermediate << "\n";
          errors.first += intermediate;

          if( _verbosity >= 1 ) std::cout << " " << _collapse.mapping().n
                                          << ". Prediction Errors: ";
          intermediate = check_one( _collapse, _strs[ _collapse.mapping().n],
                                    _collapse.mapping().n, _verbosity >= 2 );
          if( _verbosity >= 1 ) std::cout << intermediate << "\n";
          errors.second += intermediate;
        }
        return errors;
        __TRYEND(,"Error in CE::Methods::leave_one_out().\n" )
      }

    template< class T_COLLAPSE >
      opt::ErrorTuple check_one( const T_COLLAPSE &_collapse,
                                 const Crystal::Structure &_structure,
                                 size_t _n, bool _verbose )
      {
        __DEBUGTRYBEGIN
        namespace fs = boost::filesystem;
        namespace bblas = boost::numeric::ublas;
        const types::t_real weight = _structure.weight;
        const std::string name = fs::path( _structure.name ).leaf() ;
     
        types::t_real predic(0);
        const bblas::range erange( _collapse.mapping().range(_n) );
        for( bblas::range::const_iterator i( erange.begin() ); i != erange.end(); ++i )
        {
          typedef bblas::matrix_column< const typename T_COLLAPSE :: t_Matrix > 
            t_Column;
          t_Column column( _collapse.configurations(), *i );
          predic +=   _collapse.separables()( column )
                    * _collapse.mapping().eweight( _n, *i - erange.start() );
        }
  
        opt::ErrorTuple error( _structure.energy - predic,  weight );
        if( _verbose )
          std::cout << "  structure: " << std::setw(30) << name << "  "
                    << "Target: " << std::fixed << std::setw(8) 
                    << std::setprecision(2) << _structure.energy << " "
                    << "Separable: " << std::fixed << std::setw(8)
                    << std::setprecision(2) << predic << "   "
                    << "|Target-Separable| * weight: "
                    << std::fixed << std::setw(10) << std::setprecision(3) 
                    << error.mean()
                    << "\n";
        return error;
        __DEBUGTRYEND(, "Error in Fit::check_one().\n" )
      }

    template< class T_COLLAPSE, class T_STRUCTURES >
      opt::ErrorTuple check_all( const T_COLLAPSE &_collapse,
                                 const T_STRUCTURES &_strs,
                                 bool _verbose = false )
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
        __TRYBEGIN
        namespace bblas = boost::numeric::ublas;
        typedef typename T_COLLAPSE::t_Matrix t_Matrix;
        types::t_real convergence( 1e1 * tolerance );
        if( verbose ) std::cout << "Starting Alternating-least-square fit.\n";
        types::t_unsigned iter = 0;
        const size_t D( _solution.size2() );
        t_Matrix A( _collapse.dof(), _collapse.dof() );
        t_Vector b( _collapse.dof() );
        opt::ErrorTuple errors( _collapse.evaluate() );
        _collapse.update_all();
        if( verbose ) std::cout << "Allsq start: " << errors << "\n";
        do
        {
          for(size_t dim(0); dim < D; ++dim )
          {
            bblas::matrix_column<t_Matrix> column( _solution, dim );
            _collapse( A, b, dim );
            linear_solver( A, column, b );
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
        __TRYBEGIN
        std::string name = _node.Value();
        const TiXmlElement *parent = &_node;
        if( name.compare( "Allsq" ) ) 
         parent = _node.FirstChildElement( "Allsq" );
        __DOASSERT( not parent, "Could not find Allsq tag in input.\n" )
        if( parent->Attribute( "tolerance" ) )
          parent->Attribute( "tolerance", &tolerance );
        if( parent->Attribute( "itermax" ) )
          parent->Attribute( "itermax", &itermax );
        __TRYEND(,"Error in AlternatingLeastSquare::Load()\n" )
      }

} // end of Fitting namespace
