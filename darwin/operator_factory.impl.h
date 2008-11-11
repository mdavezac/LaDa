//
//  Version: $Id$
//

namespace LaDa
{
  namespace GA
  {
    namespace Factory
    {

        template< class T_INDIVIDUAL, class T_FACTORY >
          eoGenOp<T_INDIVIDUAL>* read_operator( T_FACTORY& _factory,
                                                const TiXmlElement &_node )
          {
            const std::string key( child->Value() );
            if( not _factory.exists( key ) ) return NULL;
            eoGenOpt<t_Individual>* _factory( name, _node, _factory );
            
            const TiXmlAttribute *att = _parent.FirstAttribute();
            for(; att; att = att->Next() )
            {
              const std::string attkey( att->Name() );
              const std::string attvalue( att->Value() );
              if( _factory.exists_att( attkey )  )
                _result = _factory( attname, attvalue, _result, _factory );
            }

            return result;
          }
    }
  }
}

