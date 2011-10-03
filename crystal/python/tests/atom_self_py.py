from lada.crystal import @TYPE@
from atom_self_@TYPE@ import get_new_object, get_new_atom, get_reference,\
                              get_incorrect_reference, check_back_reference, get_pointer,\
                              change_pos

a = @TYPE@()
assert id(a) == id(a._self)
check_back_reference(a, a)
a.magnetic = 5
assert getattr(a._self, "magnetic", None) == 5

# a = get_new_object()
# assert id(a) == id(a._self)
# check_back_reference(a, a)
# a.magnetic = 5
# assert getattr(a._self, "magnetic", None) == 5

# a = get_new_atom()
# assert id(a) == id(a._self)
# check_back_reference(a, a)
# a.magnetic = 5
# assert getattr(a._self, "magnetic", None) == 5

# a = get_incorrect_reference()
# assert a._self == None
# a.pos[0] = 1
# b = get_incorrect_reference()
# print abs(b.pos[0] - 0e0)

a = get_reference()
assert id(a) == id(a._self)
check_back_reference(a, a)
b = get_reference()
assert id(a) == id(b)
check_back_reference(a, b)
check_back_reference(b, a)
a.magnetic = 5
assert getattr(a._self, "magnetic", None) == 5
assert getattr(b, "magnetic", None) == 5
# assert abs(a.pos[0]-0e0) < 1e-12
a.type.add("Au")
del a
print "WWW"
# assert abs(b.pos[0]-0e0) < 1e-12
print "WWW2"
# change_pos(a)
print "WWW3"
# assert abs(b.pos[0]-1e0) < 1e-12
print "WWW4"

# a = get_pointer()
# assert id(a) == id(a._self)
# check_back_reference(a, a)
# change_pos(a)
# assert abs(a.pos[0]-1e0) < 1e-12
