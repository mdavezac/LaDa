#ifndef _CHECKPOINT_H_
#define _CHECKPOINT_H_

namespace LaDa 
{
  template<class Object>
  class Monitor : public eoUpdater
  {
    private:
      Object* object;
    public:
      Monitor( Object* _object) { object =_object; }
      virtual void lastCall()
        { object->print_xml(); }
      virtual void operator()(void)
        { object->print_xmgrace(); }
      
      virtual std::string className(void) const { return "Monitor"; } 
  };

} // namespace LaDa


#endif //  _CHECKPOINT_H_
