//
//  Version: $Id: print_tree.h 1136 2009-05-22 00:34:43Z davezac $
//

#ifndef _LADA_LOADNSAVE_PRINT_TREE_H_
#define _LADA_LOADNSAVE_PRINT_TREE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/iterator_range.hpp>
#include <boost/fusion/include/begin.hpp>  
#include <boost/fusion/include/prior.hpp>  
#include <boost/fusion/include/end.hpp>  
#include <boost/fusion/include/is_view.hpp>  

#include <boost/mpl/size.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/static_assert.hpp>


namespace LaDa 
{
  namespace load_n_save
  {
    namespace fusion = boost::fusion;

    class Print 
    {
        template< class T >
          struct is_false : public boost::is_same< T, boost::mpl::bool_<false> > {};
      public:
        Print() : indent( "" ), tab( "  " ) {}
        Print(std::string &_i, const std::string &_tab) : indent( _i ), tab( _tab ) {}
        Print( const Print &_c) : indent( _c.indent ), tab( _c.tab ) {}

        template< class T >
          typename boost::disable_if< is_false<T> > :: type
            operator()( const T& _a ) const
            {
              namespace bm = boost::mpl;
              namespace bf = boost::fusion;
              typedef typename bf::result_of::as_vector< T > :: type const t_Vector;
              t_Vector vector( fusion::as_vector( _a ) );
              bf :: for_each( vector, section_(indent, tab) );
            }
        void operator()( const boost::mpl::bool_<false>& ) const
          { std::cout << "none\n"; }

      protected:
        struct subsections_;
        struct section_;
        struct options_;
        struct option_;
        struct values_;
        struct value_;
        struct action_;

        mutable std::string indent;
        std::string tab;
    };
    
    struct Print :: section_
    {
      section_( std::string &_i, const std::string &_t ) : indent(_i), tab(_t) {};
      section_( const section_ &_c ) : indent(_c.indent), tab(_c.tab) {};
      template< class T >
        typename boost::disable_if< is_false<T> > :: type
          operator()( const T& _expr ) const
          {
            namespace bf = boost::fusion;
            namespace bm = boost::mpl;
            typedef typename bf::result_of::as_vector< T > :: type const t_Vector;
            t_Vector expr( fusion::as_vector( _expr ) );
            
            BOOST_STATIC_ASSERT
            (( 
              bm::equal_to< typename bm::size<t_Vector>::type, bm::int_<4> > :: type :: value
            ));
            std::cout << indent << "section " << bf::at_c<0>( expr ) << "\n";
            indent += tab;
          
            std::cout << indent << "action: ";
            action_()( bf::at_c<3>(expr) );
            std::cout << "\n";
          
            options_( indent, tab )( bf::at_c<1>(expr) );
          
            subsections_(indent, tab)( bf::at_c<2>(expr) );
          
            indent = indent.substr( 0, indent.size() - tab.size() );
          }

      mutable std::string &indent;
      const std::string &tab;
    };

    struct Print :: options_ 
    {
      options_( std::string &_i, const std::string &_t ) : indent(_i), tab(_t) {};
      options_( const options_ &_c ) : indent(_c.indent), tab(_c.tab) {};
      template< class T >
        typename boost::disable_if< is_false<T> > :: type
          operator()( const T& _expr ) const
          {
            namespace bf = boost::fusion;
            namespace bm = boost::mpl;
            typedef typename bf::result_of::as_vector< T > :: type const t_Vector;
            t_Vector vector( fusion::as_vector( _expr ) );
            bf :: for_each( vector, option_(indent, tab) );
          }
      void operator()( const boost::mpl::bool_<false>& ) const
        { std::cout << indent << "options: none\n"; }

      mutable std::string &indent;
      const std::string &tab;
    };

    struct Print :: option_ 
    {
      option_( std::string &_i, const std::string &_t ) : indent(_i), tab(_t) {};
      option_( const option_ &_c ) : indent(_c.indent), tab(_c.tab) {};
      template< class T >
        typename boost::disable_if< is_false<T> > :: type
          operator()( const T& _expr ) const
          {
            namespace bf = boost::fusion;
            namespace bm = boost::mpl;
            typedef typename bf::result_of::as_vector< T > :: type const t_Vector;
            t_Vector vector( fusion::as_vector( _expr ) );
            std::cout << indent << "option[";
            action_()( bf::at_c<2>(vector) );
            std::cout << "] " << bf::at_c<0>( vector ) << ": ";
            values_()( bf::at_c<1>(vector) );
            std::cout << "\n";
          }
      void operator()( const boost::mpl::bool_<false>& ) const
        { std::cout << indent << "options: none\n"; }

      mutable std::string &indent;
      const std::string &tab;
    };

    struct Print :: values_ 
    {
      template< class T >
        typename boost::disable_if< is_false<T> > :: type
          operator()( const T& _a ) const
          {
            namespace bm = boost::mpl;
            namespace bf = boost::fusion;
            typedef typename bf::result_of::as_vector< T > :: type const t_Vector;
            typedef typename bf::result_of::begin< t_Vector > :: type t_First;
            typedef typename bf::result_of::end< t_Vector > :: type t_End;
            typedef typename bf::result_of::prior< t_End > :: type t_Prior;
            t_Vector vector( fusion::as_vector( _a ) );
            t_First const first(vector);
            t_Prior const prior(vector);
            bf::iterator_range< t_First, t_Prior > const v( first, prior );
            bf :: for_each( v, value_() );
            value_("")( bf::deref( prior ) );
          }
      void operator()( const boost::mpl::bool_<false>& ) const {}
    };

    struct Print :: value_ 
    {
      template< class T > struct isview : public boost::mpl::bool_<false> {};
      template< class T1, class T2 > 
        struct isview< boost::fusion::joint_view<T1, T2> >
          : public boost::mpl::bool_<true> {};

      value_( const std::string& _sep = ", " ) : sep( _sep ) {};
      value_( const value_& _c ) : sep( _c.sep ) {};
      template< class T >
        typename boost::disable_if< boost::mpl::or_< is_false<T>, isview<T> > > :: type
          operator()( const T& _a ) const 
            { std::cout << _a << sep; }

      template< class T1, class T2 >
        void operator()( boost::fusion::joint_view< T1, T2> const & _a ) const
        {
          namespace bf = boost::fusion;
          typedef bf::joint_view<T1, T2> T;
          typedef typename bf::result_of::as_vector< T > :: type const t_Vector;
          t_Vector vector( fusion::as_vector( _a ) );
          std::cout << bf::at_c<0>( vector ) << " == " << bf::at_c<1>( vector ) << sep;
        }
      void operator()( const boost::mpl::bool_<false>& ) const {}
      const std::string sep;
    };
    

    struct Print :: subsections_ 
    {
      subsections_( std::string &_i, const std::string &_t ) : indent(_i), tab(_t) {};
      subsections_( const subsections_ &_c ) : indent(_c.indent), tab(_c.tab) {};
      template< class T >
        typename boost::disable_if< is_false<T> > :: type
          operator()( const T& _a ) const
          {
            namespace bf = boost::fusion;
            typedef typename bf::result_of::as_vector< T > :: type const t_Vector;
            t_Vector vector( fusion::as_vector( _a ) );
            bf :: for_each( vector, section_(indent, tab) );
          }
      void operator()( const boost::mpl::bool_<false>& ) const
        { std::cout << indent << "subsections: none\n"; }

      mutable std::string &indent;
      const std::string &tab;
    };

    struct Print :: action_ 
    {
      template< class T >
        typename boost::disable_if< is_false<T> > :: type
          operator()( const T& _expr ) const
            { std::cout << _expr; }
      void operator()( const boost::mpl::bool_<false>& ) const
        { std::cout << "none"; }
    };
 
  } // namespace load_n_save

} // namespace LaDa


#endif
