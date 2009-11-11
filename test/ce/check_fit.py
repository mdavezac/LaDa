#
#  Version: $Id$
#


def main():

  from lada import crystal, ce

  lattice = crystal.Lattice("input.xml")
  mlclasses = ce.MLClusterClasses("input.xml", False)
  fit = ce.Fit(mlclasses)

  fit.read_directory("data")

if __name__ == "__main__":
  main()
