#include <boost/lambda/lambda.hpp>

namespace LaDa
{
  namespace crystal
  {
    template< class T_TYPE >
      std::ostream& operator<<(std::ostream &_stream, TemplateStructure<T_TYPE> const &_str)
        { return _stream << *_str.impl_; }

    template<class TYPE> template<class T_ARCHIVE>
      bool TemplateStructure :: lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version) 
      {
        if( _ar.is_loading() )
        {
          boost::shared_ptr< StructureData<TYPE> > dummy(impl_);
          impl_.reset(new StructureData<TYPE>());
        }
        return _ar & *impl_;
      }
    
    //! \cond
    void  find_range( const math::rMatrix3d &A, math::iVector3d &kvec );

    template <class CONTAINER>
    void remove_equivalents( CONTAINER &_cont, const math::rMatrix3d &_cell)
    {
      typename CONTAINER :: iterator i_vec = _cont.begin();
      typename CONTAINER :: iterator i_end = _cont.end();
      typename CONTAINER :: iterator i_which;

      while( i_vec != i_end )
      {
        i_which = i_vec+1;
        for ( ; i_which != i_end; i_which++ )
          if ( are_equivalent( *i_which, *i_vec, _cell ) )
            break;

        if ( i_which == i_end )
        { 
          ++i_vec;
          continue;
        }
        
        ( i_vec->pos.squaredNorm() < i_which->pos.squaredNorm() ) ?  
              _cont.erase(i_which): _cont.erase(i_vec);
        i_vec = _cont.begin();
        i_end = _cont.end();
      }

    }
    //! \endcond


    template< class T_TYPE >
      bool TStructure<T_TYPE> :: set_site_indices()
      {
        if ( not lattice ) return false;
      
        bool result = true;
        typename TStructure<T_TYPE> :: t_Atoms :: iterator i_atom = atoms.begin();
        typename TStructure<T_TYPE> :: t_Atoms :: iterator i_atom_end = atoms.end();
        for(; i_atom != i_atom_end; ++i_atom )
        {
          i_atom->site = lattice->get_atom_site_index( i_atom->pos );
          (i_atom->site == -1) ?
            result = false:
            i_atom->freeze |= lattice->sites[ i_atom->site ].freeze;
        }
        return result;
      }

    template< class T_TYPE >
     bool TStructure<T_TYPE> :: operator== (const TStructure<T_TYPE> &_str ) const
     {
       return     math::eq( cell(0,0), _str.cell(0,0) )
              and math::eq( cell(1,0), _str.cell(1,0) )
              and math::eq( cell(2,0), _str.cell(2,0) )
              and math::eq( cell(0,1), _str.cell(0,1) )
              and math::eq( cell(1,1), _str.cell(1,1) )
              and math::eq( cell(2,1), _str.cell(2,1) )
              and math::eq( cell(0,2), _str.cell(0,2) )
              and math::eq( cell(1,2), _str.cell(1,2) )
              and math::eq( cell(2,2), _str.cell(2,2) )
              and atoms == _str.atoms;
     }
      
    template<class t_container >
      void Structure :: set_kvectors( const t_container &_container )
      {
        typename t_container :: const_iterator i_kvec =  _container.begin();
        typename t_container :: const_iterator i_end =  _container.end();
        CAtom kvec;
        k_vecs.clear();
        k_vecs.reserve( _container.size() );
        for( ; i_kvec != i_end; ++i_kvec ) 
        {
          kvec = (*i_kvec);
          k_vecs.push_back( kvec );
        }
      }


  }
} // namespace LaDa
