#ifndef LADA_RW_MISC_H
#define LADA_RW_MISC_H

#include <iostream>

namespace LaDa
{
  namespace rw
  {
    //! Output stream.
    class Ouput : protected tree::Node
    {
      public:
        //! Constructor.
        Ouput(t_string const &_name = "") : tree::Node(_name) {};
        //! Destructor.
        virtual ~Output() {};

        //! Adds a new node.
        template<class T> Ouput& operator<<(Node<T> const &_object); 
        //! Adds a new attribute.
        template<class T> Ouput& operator<<(Attribute<T> const &_object); 
        //! Creates a string from an object.
        template<class T> t_string to_string(T const &_object);
    };

    template<class T> 
      Output& Output::operator<<(Node<T> const &_object)
      {
        Ouput node(_object.name);
        _object.rw(node);
        branches->push_back(node);
        return *this;
      }

    template<class T> 
      Output& Output::operator<<(Attribute<T> const &_object)
      {
        attributes->push_back( tree::Node::t_Attribute(_object.name, to_stream(_object.value)) );
        return *this;
      }

    template<class T> 
      t_string Output::to_string(T const &_object)
      {
        std::ostringstream sstr;
        sstr << _object;
        return sstr.str();
      }
  }
}
#endif
