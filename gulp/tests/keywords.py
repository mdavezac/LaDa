def test_optikeywords():
  from pylada.tools.input.block import AttrBlock
  from pylada.gulp.keywords import OptiKeyword
  
  # fake class to try out opti keywords
  class A(AttrBlock):
    def __init__(self):
      super(A, self).__init__()
      self.opti = False
      self.conv = OptiKeyword()

  a = A()
  assert a.opti == False
  assert getattr(a._input['conv'], 'keyword', 'false') == 'conv'
  assert a.conv is None
  assert a._input['conv'].output_map(gulp=a) is None
  
  a.opti = True
  assert a._input['conv'].output_map(gulp=a) is None

  a.conv = True
  assert a._input['conv'].output_map(gulp=a).get('conv', 'False') == 'True'
  a.opti = False
  assert a._input['conv'].output_map(gulp=a) is None

def test_twobody():
  from pylada.gulp.keywords import TwoBody

  a = TwoBody()
  a.keyword = 'morse'
  assert len(a) == 0
  a['Ti', 'O'] = range(5)
  assert len(a) == 1
  assert a.__iter__().next() == sorted(['Ti|core', 'O|core'])
  assert a['Ti', 'O'] == range(5)
  assert a['Ti', 'O'] is a['Ti|core', 'O']
  assert a['Ti', 'O'] is a['Ti', 'O|core']
  assert a['Ti', 'O'] is a['Ti|core', 'O|core']
 
  a['Ti|shell', 'O'] = range(3)
  assert a['Ti', 'O'] == range(5)
  assert a['Ti|shell', 'O'] != range(5)
  assert a['Ti|shell', 'O'] == range(3)
  assert a['Ti|shell', 'O'] is a['Ti|shell', 'O|core']
  assert a['Ti|shell', 'O'] is a['Ti shell', 'O|core']
  assert a['Ti|shell', 'O'] is a['Tishell', 'O core']

  assert a.output_map() is None
  a.enabled = True
  assert 'morse' in a.output_map()
  assert a.output_map()['morse'].split()                                       \
         == ['O', 'core', 'Ti', 'shell'] + [str(float(u)) for u in range(3)]   \
            + ['O', 'core', 'Ti', 'core'] + [str(float(u)) for u in range(5)]


if __name__ == '__main__': 
  test_optikeywords()
  test_twobody()
