from numpy import array
from lada.math import check_extract, check_frompy_to_cpp_vec, check_frompy_to_cpp_mat,\
                      check_topy_from_cpp_mat, check_topy_from_cpp_vec

mat = array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])
vec = array([-1,0.5,0.5])

check_extract(mat)
print "extraction - OK"

check_frompy_to_cpp_mat(mat)
print "automatic matrix extraction - OK"
check_frompy_to_cpp_vec(vec)
print "automatic vector extraction - OK"

print check_topy_from_cpp_mat()
print "automatic matrix conversion - OK"
print check_topy_from_cpp_vec()
print "automatic vector conversion - OK"

