#ifndef LADA_RW_MISC_H
#define LADA_RW_MISC_H

#include <vector>
#include <string>

#include <boost/shared_ptr.hpp>

namespace LaDa
{
  namespace rw
  {
    //! String type across the whole library.
    typedef std::string t_string;
    //! Node tag.
    //! Opens a constant new node in the i/o stream.
    template<class T> struct Node
    {
      //! Type of the reference.
      typedef T & reference_type;
      //! Name of the node.
      t_string name;
      //! Object to reference.
      reference_type value;
    };

    //! Attribute tag.
    //! Opens a new attribute in the i/o stream.
    template<class T> struct Attribute
    {
      //! Type of the reference.
      typedef T & reference_type;
      //! Name of the attribute.
      t_string name;
      //! Object to reference.
      reference_type value;
    };

    namespace tree
    {
      //! Contains all details of a node.
      struct Node
      {
        //! Type of the inner nodes.
        typedef std::vector<Nodes> t_Branches;
        //! Type of an attribute.
        typedef std::pair<t_string, t_string> t_Attribute;
        //! Type of a list of attributes.
        typedef std::vector<t_Attribute> t_Attributes;

        //! Name of this node.
        t_string name;
        //! Inner nodes.
        boost::shared_ptr<t_Branches> branches;
        //! Attributes in this node.
        boost::shared_ptr<t_Attributes> attributes

        //! Constructs a node.
        Node   (t_string const &_name)
             : name(_name), branches(new t_Branches), attributes(new t_Attributes) {};
        //! Copies a node.
        Node   (Node const &_c)
             : name(_c.name), branches(_c.branches), attributes(_c.attributes) {};
        //! Destructor.
        ~Node() {}
      };
    }
  }
}
#endif
