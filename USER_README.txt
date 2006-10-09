La(marck)Da(rwin) code combines va and ga to obtain the groundstates of a cluster expansion

it is based upon the following libraries

  _ liblamarck (NREL VA library )
  _ opt++ (external, for minization) 
  _ opt (NREL optimization library )
  _ EO, evolutionary object (external GA library )

make sure you have everything installed

then edit dependencie.pl and change include directories, library locations, and compilers to suit your need. 

run dependencies.pl and compile!

LaDa specific input follows, and should be included with liblamarck input

<?xml version="1.0" ?>
<LaDa>
  <Structure>
    <!-- structure on which to perform a LaDa search -->
    <!-- see liblamarck input description -->
  </Structure>

  <!-- GA specific stuff -->
  <!-- at this point, methods are "debug", "lamarck" for a GA evolution, and "multistart" -->
  <!-- maxgen is the maximum number of generations -->
  <GA method = "?" maxgen = "?" >

    <!-- the minimizer is either --> 
       <!-- wang: opt++ dependant minimizer with Wang, Berathan et al constraints -->
       <!-- linear: or own physical minimizer, see opt_minimize_linear.h for description -->
       <!-- sa: a T=0 simmulated annealing method -->
       <!-- maxeval limits the number of evaluation and gradients of the last two -->
    <Minimizer type="wang|linear|sa" maxeval="?">

    <!-- breeding operators --> 
      <!-- can be either Mutation and/or Crossover --> 
      <!-- and/or is specified in type -->
    <Operators type="or" prob="0.75" >
      <Mutation prob=1.0 />
      <Crossover prob=0.5 />
    </Operators>

    <!-- deterministic tournament size when selecting parents --> 
    <Selection size=2 />

    <!-- population size --> 
    <Population size="10"/>

    <!-- number of offsprings per generation = rate * population size --> 
    <Offsprings rate="0.25"/>

    <!-- if OnePointHull is present, finds lowest energy state, rather than a true convex hull -->
    <OnePointHull/>
  </GA>

  <!-- you need to specify Lattice, Clusters and CS in this file as well -->
  <!-- and put it between the LaDa tags -->
  <!-- see liblamarck for description -->

</LaDa>
