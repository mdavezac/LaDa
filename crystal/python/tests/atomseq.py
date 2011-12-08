from numpy import all
from lada.crystal.cppwrappers.atom import AtomSequence as AtomSeq

# Try correct initialization. 
a = AtomSeq()
assert a.type == "" and a.site == -1 and a.freeze == 0 and all(a.pos == 0) and len(a.__dict__) == 0
a = AtomSeq('Au')
assert a.type == "Au" and a.site == -1 and a.freeze == 0 and all(a.pos == 0) and len(a.__dict__) == 0
a = AtomSeq(0.1, 0.1, 0.1, 'Au')
assert a.type == "Au" and a.site == -1 and a.freeze == 0 and all(a.pos == 0.1) and len(a.__dict__) == 0
a = AtomSeq([0.1, 0.1, 0.1], 'Au')
assert a.type == "Au" and a.site == -1 and a.freeze == 0 and all(a.pos == 0.1) and len(a.__dict__) == 0
a = AtomSeq([0.1, 0.1, 0.1], type='Au')
assert a.type == "Au" and a.site == -1 and a.freeze == 0 and all(a.pos == 0.1) and len(a.__dict__) == 0
a = AtomSeq(type='Au', pos=[0.1, 0.1, 0.1])
assert a.type == "Au" and a.site == -1 and a.freeze == 0 and all(a.pos == 0.1) and len(a.__dict__) == 0
a = AtomSeq('Au', pos=[0.1, 0.1, 0.1])
assert a.type == "Au" and a.site == -1 and a.freeze == 0 and all(a.pos == 0.1) and len(a.__dict__) == 0
a = AtomSeq('Au', pos=[0.1, 0.1, 0.1], m=1.3)
assert a.type == "Au" and a.site == -1 and a.freeze == 0 and all(a.pos == 0.1)\
       and len(a.__dict__) == 1 and abs(a.__dict__.get('m', 0) - 1.3) < 1e-12
a = AtomSeq(['Au', 'Pd'])
assert a.type == ["Au", 'Pd'] and a.site == -1 and a.freeze == 0 and all(a.pos == 0) and len(a.__dict__) == 0
a = AtomSeq(['Au', 'Pd'], pos=[0.1, 0.1, 0.1])
assert a.type == ["Au", 'Pd'] and a.site == -1 and a.freeze == 0 \
       and all(abs(a.pos - 0.1) < 1e-12) and len(a.__dict__) == 0
# Try incorrect initialization. Check for garbage collection.
try: a = AtomSeq([0.1, 0.1])
except TypeError: pass
else: raise RuntimeError("Should have failed.")
try: a = AtomSeq([0.1, 0.1, 0.1, 0.1])
except TypeError: pass
else: raise RuntimeError("Should have failed.")
try: a = AtomSeq(0.1, 0.1, 0.1, 0.1)
except TypeError: pass
else: raise RuntimeError("Should have failed.")
try: a = AtomSeq(0.1, 'Au', 0.1, 0.1)
except TypeError: pass
else: raise RuntimeError("Should have failed.")
try: a = AtomSeq('Au', type='Au')
except TypeError: pass
else: raise RuntimeError("Should have failed.")
try: a = AtomSeq(0.1, 0.1, 0.1, pos=[0.1, 0.1, 0.1])
except TypeError: pass
else: raise RuntimeError("Should have failed.")

assert repr(AtomSeq('Au', pos=[1, 1, 1], m=1))     == "AtomSequence(1, 1, 1, 'Au', m=1)"
assert str(AtomSeq('Au', pos=[1, 1, 1], m=1))      == "AtomSequence(1, 1, 1, 'Au', m=1)"
assert str(AtomSeq('Au', pos=[1, 1, 1], freeze=1)) == "AtomSequence(1, 1, 1, 'Au', freeze=1)"
assert str(AtomSeq('Au', pos=[1, 1, 1], site=1))   == "AtomSequence(1, 1, 1, 'Au', site=1)"

from copy import copy, deepcopy
b = AtomSeq(['Au', 'Pd'], m=0)
a = copy(b)
b.type = ['Au']
assert a.type == ["Au"] and a.site == -1 and a.freeze == 0 and all(a.pos == 0)
assert getattr(a, 'm', 1) == 0
b = AtomSeq(['Au', 'Pd'], m=0)
a = deepcopy(b)
b.type = ['Pd']
b.pos += 1
del b.m 
assert a.type == ["Au", 'Pd'] and a.site == -1 and a.freeze == 0 and all(a.pos == 0)
assert getattr(a, 'm', 1) == 0
assert b.type == ['Pd'] and b.site == -1 and b.freeze == 0 and all(b.pos == 1) and len(b.__dict__) == 0

a = a.to_dict()
assert a['type'] == ["Au", 'Pd'] and a['site']== -1 and a['freeze'] == 0\
       and all(a['pos']== 0) and a['m'] == 0

a = AtomSeq(0.1, 0.1, 0.1, ['Au', 'Pd'], m=0, site=0, freeze=1).cast()
assert a.__class__.__name__ == 'AtomStr' and a.type == 'Au, Pd' and a.site == 0 \
       and a.freeze == 1 and all(abs(a.pos - 0.1) < 1e-12) and len(a.__dict__) == 1 \
       and getattr(a, 'm', 1) == 0

b = AtomSeq(['Au', 'Pd'], m=0)
a = b.copy()
b.type = ['Pd']
b.pos += 1
assert a.type == ["Au", 'Pd'] and a.site == -1 and a.freeze == 0 and all(a.pos == 0)
assert getattr(a, 'm', 1) == 0

from pickle import loads, dumps
a = AtomSeq(pos=[0, 1, 2], type=["Au", 'Pd'], site=0, freeze=1, m=6)
b = loads(dumps(a))
assert all(abs(b.pos - [0, 1, 2]) < 1e-12) and b.type == ['Au', 'Pd'] and b.site == 0 and b.freeze == 1\
       and len(b.__dict__) == 1 and b.__dict__.get('m', None) == 6
b = loads(dumps((a, a)))
assert b[0] is b[1]
