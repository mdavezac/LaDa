#
#  Version: $Id$
#


def main():

  from lada import crystal, ce

  lattice = crystal.Lattice("data/ce.xml")
  fit = ce.Fit()

  fit.read_directory("data")
  print fit._data 

if __name__ == "__main__":
  main()
