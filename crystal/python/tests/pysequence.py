from lada.crystal.cppwrappers.atom import _new_sequence as Sequence
from lada import error

# test sequence to list transformation.
assert list(Sequence()) == []
assert list(Sequence(['Au', 'Pd'])) == ['Au', 'Pd'];

# test all possible comparison operators.
for a in [ [], ['Au'], ['Au', 'Pd'], ['Au', 'P'], ['Au', 'Pu'], ['Au', 'Pd', 'Cr'] ]:
  astr, aseq = "".join(a), Sequence(a)
  assert str(a) == str(aseq)
  for b in [ [], ['Au'], ['Au', 'Pd'], ['Au', 'P'], ['Au', 'Pu'], ['Au', 'Pd', 'Cr'] ]:
    bstr, bseq = "".join(b), Sequence(b)
    assert (astr == bstr) == (aseq == bseq)
    assert (astr == bstr) == (a    == bseq)
    assert (astr == bstr) == (aseq == b)

    assert (astr < bstr) == (aseq < bseq)
    assert (astr < bstr) == (a    < bseq)
    assert (astr < bstr) == (aseq < b)

    assert (astr <= bstr) == (aseq <= bseq)
    assert (astr <= bstr) == (a    <= bseq)
    assert (astr <= bstr) == (aseq <= b)

    assert (astr > bstr) == (aseq > bseq)
    assert (astr > bstr) == (a    > bseq)
    assert (astr > bstr) == (aseq > b)

    assert (astr >= bstr) == (aseq >= bseq)
    assert (astr >= bstr) == (a    >= bseq)
    assert (astr >= bstr) == (aseq >= b)

# test sequence methods
for a in [ [], ['Au'], ['Au', 'Pd'], ['Au', 'P'], ['Au', 'Pu'], ['Au', 'Pd', 'Cr'] ]:
  aseq = Sequence(a)
  assert len(a) == len(aseq)
  assert list(a) == a
  assert aseq.copy() == a
  assert a * 2 == aseq * 2
  aseq *= 2 
  assert a * 2 == aseq
  aseq = Sequence(a)
  for b in [ [], ['Au'], ['Au', 'Pd'], ['Au', 'P'], ['Au', 'Pu'], ['Au', 'Pd', 'Cr'] ]:
    bseq = Sequence(b)
    assert aseq + bseq == a + b
    assert aseq +    b == a + b
    aseq += b
    assert aseq == a + b
    aseq = Sequence(a)
    aseq.extend(bseq)
    assert aseq == a + b
    aseq = Sequence(a)
for a in [ [], ['Au'], ['Au', 'Pd'], ['Au', 'P'], ['Au', 'Pu'], ['Au', 'Pd', 'Cr'] ]:
  aseq = Sequence(a)
  for i, s in enumerate(a): 
    assert s in aseq
    assert aseq[i] == s
    assert aseq[-(len(a)-i)] == s
    aseq[i] = "V"
    assert "V" in aseq
    assert aseq[i] == "V"
    aseq[-(len(aseq)-i)] = s
    assert s in aseq
    assert aseq[i] == s
  assert "V" not in aseq
  for s0, s1 in zip(a, aseq): 
    assert s0 == s1

# test member functions.
a = Sequence()
a.append('Au')
assert a == ['Au']
a.extend(['Pd', 'V']);
assert a == ['Au', 'Pd', 'V']
assert a.index('Au') == 0
assert a.index('Pd') == 1
assert a.index('V') == 2
a.pop(1)
try: a.index('Pd')
except ValueError: pass
try: a.index('Pd')
except error.ValueError: pass
except ValueError: raise RuntimeError("SHIT")
try: a.pop(len(a) + 1)
except IndexError: pass
try: a.pop(len(a) + 1)
except error.IndexError: pass
except IndexError: raise RuntimeError("SHIT")
try: a.pop(-len(a) - 1)
except IndexError: pass
try: a.pop(-len(a) - 1)
except error.IndexError: pass
except IndexError: raise RuntimeError("SHIT")
a.insert(1, 'Pd')
assert a == ['Au', 'Pd', 'V']
assert a.index('Au') == 0
assert a.index('Pd') == 1
assert a.index('V') == 2
a.insert(len(a), 'Pd')
assert a == ['Au', 'Pd', 'V', 'Pd']
a.pop(-1)
assert a == ['Au', 'Pd', 'V']
a.clear()
assert len(a) == 0
b = ['Au', 'V', 'Pd', 'Cr']
a = Sequence(b)
a.reverse()
b.reverse()
assert a == b
a.sort()
b.sort()
assert a == b


# test some expected exceptions.
try: a[-len(a)-1]
except IndexError: pass
try: a[len(a)]
except error.IndexError: pass
except IndexError: raise RuntimeError("SHIT")
try: a[-len(a)-1] = '0'
except IndexError: pass
try: a[len(a)] = '0'
except error.IndexError: pass
except IndexError: raise RuntimeError("SHIT")
try: a[0] = 0
except TypeError: pass

# sequences do not accept other attributes. No __dict__. 
try: a.this = "0"
except AttributeError: pass
