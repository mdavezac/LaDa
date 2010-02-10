from _params import Standard, NoPrintStandard, Iniwave, Algo
class Test(object):

  standard = Standard("KEY", 1.5, validity = lambda x: x > 1 and x < 5)
  noprint = NoPrintStandard("KEY", 1.5, validity = lambda x: x > 1 and x < 5)
  iniwave = Iniwave()
  algo = Algo()

  def __getstate__(self): return self.standard, self.noprint, self.iniwave, self.algo
  def __setstate__(self, arg): 
    self.standard, self.noprint, self.iniwave, self.algo = arg

  def __init__(self): pass

  def test_serial(self):
    def check(self, name, other):
      from pickle import dumps, loads
      
      string = dumps(self)
      e = loads(string)
  
      assert getattr(self, name).incar_string(None) == getattr(e, name).incar_string(None),\
             "%s == %s" % ( getattr(self, name).incar_string(None), \
                            getattr(e, name).incar_string(None) )
            
      getattr(e, name).value = other
      assert getattr(self, name).incar_string(None) != getattr(e, name).incar_string(None),\
             "%s != %s" % ( getattr(self, name).incar_string(None), \
                            getattr(e, name).incar_string(None) )
  
    this = [ ("standard", 2.0), ("noprint", 2.0), ("iniwave", "jellium"), 
             ("algo", "fAst") ]
    for name, other in this:
      yield check, self, name, other

  def test_serial_invalidity(self):

    def check(self, name, invalids):
      from pickle import dumps, loads
      
      string = dumps(self)
      e = loads(string)
  
      for i in invalids:
        try: getattr(e, name).value = i
        except ValueError: pass
        else: assert False

    this = [ ("standard", (1.0, 5.0)), ("noprint", (1.0,5.0)),
             ("iniwave", "whatever"), ("algo", "nor mal") ]
    for name, invalids in this:
      yield check, self, name, invalids

  def test_noprint(self):
    from copy import deepcopy
    e = deepcopy(self)
    e.noprint.value = 1.5
    assert e.noprint.incar_string(None)[0] == '#', "%s" % (e.noprint.incar_string(None))
    e.noprint.value = 2.5
    assert e.noprint.incar_string(None)[0] != '#', "%s" % (e.noprint.incar_string(None))

  def test_values(self):
    from copy import deepcopy

    def check(name, val, eq):
      e = deepcopy(self)
      getattr(e, name).value = val
      assert getattr(e, name).value == eq, "%s == %s" % (getattr(e, name).value, eq)
    
    jobs = [ ("iniwave", "random", 1), ("iniwave", "jellium", 0),
             ("algo", " verY-FasT", "Very_Fast"), ("algo", "nORmal", "Normal") ]
    for name, val, eq in jobs:
      yield check, name, val, eq

