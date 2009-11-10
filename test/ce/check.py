#! /usr/bin/python

def convert(_classes):
  from lada import ce

  result = ce.ClusterClasses()
  for class_ in _classes:
    dummy = ce.Clusters()
    for mlcluster in class_:
      cluster = ce.Cluster()
      cluster.eci = class_.eci
      cluster.vectors.append(mlcluster.origin.pos)
      for i, vec in enumerate(mlcluster):
        cluster.vectors.append(vec.pos)
      dummy.append(cluster)
    result.append(dummy)
  return result

def main():
  from lada import ce, crystal
  from math import fabs
  import boost.mpi as mpi
  import numpy
  import pyublas

  functional = ce.Cubic()
  functional.set_mpi(mpi.world)
  functional.load( "input.xml" )
  structure = crystal.Structure()
  structure.fromXML("input.xml")

  mlclasses = ce.MLClusterClasses("input.xml", False)
# mlclasses[0]
# print mlclasses
  classes = convert(mlclasses)

# mlclusters = convert_clusters_to_mlclusters
  tests = [ ("000011100000", -74.377010),
            ("111000000101", -83.064845), 
            ("101000000011", -80.889828), 
            ("000111111000", -84.277151), 
            ("000111111001", -71.156430), 
            ("110000000111", -83.064845), 
            ("110000000101", -80.889828), 
            ("111000000111", -84.277151), 
            ("000111101000", -83.064845), 
            ("000110111000", -83.064845), 
            ("101000000111", -83.064845), 
            ("000110011000", -80.889828), 
            ("000011101000", -80.889828), 
            ("000111011000", -83.064845), 
            ("000011111000", -83.064845), 
            ("110000000011", -80.889828), 
            ("111000100111", -71.156430), 
            ("001000000000", -27.375282), 
            ("011000000111", -83.064845), 
            ("000010000000", -27.375282), 
            ("000000100000", -27.375282), 
            ("000000000001", -27.375282), 
            ("000100001000", -55.269849), 
            ("111010000111", -71.156430), 
            ("000000000100", -27.375282), 
            ("010000000000", -27.375282), 
            ("000111111100", -71.156430), 
            ("000001100000", -55.269849), 
            ("000011110000", -80.889828), 
            ("111001100111", -58.030174), 
            ("100000000010", -55.269849), 
            ("010000000100", -55.269849), 
            ("000001110000", -74.377010), 
            ("010000000001", -55.269849) ] 

  print len(mlclasses[0]), len(structure.atoms)
  for test in tests:
    for i, atom in enumerate(structure.atoms):
      if test[0][i] == "0": atom.type = -1e0
      else:               atom.type = 1e0
    c = functional(structure) 
    print "check %f, test %f, diff %f " \
          % ( c, test[1], fabs(c-test[1]) )

if __name__ == "__main__":
  main()
