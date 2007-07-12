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
    <!-- can be applied recurrently --> 
    <!-- four types exists --> 
      <!-- Operators type="and"  applies all the operators it contains with their probability prob -->
      <!-- Operators type="or"  applies one of the operators it contains -->
         <!-- using a roulette based on their probability prob -->
      <!-- Mutation  applies point mutation with probability value -->
      <!-- Crossover  applies point crossovers between two parents, depending on probability value -->
      <!-- UtterRandom utterly randomizes a chromosome -->
      <!-- Minimizer minimizes the offspring -- see input description above -->
    <Operators type="and">
      <Operators type="or" prob="1.0" >
        <Mutation value=0.05 prob="0.25"/>
        <Crossover value=0.5 prob="0.75"/>
      </Operators>
      <Minimizer type="linear" maxeval=20 prob="1.0" />
    </Operators>

    <!-- deterministic tournament size when selecting parents --> 
    <Selection size=2 />

    <!-- population size --> 
    <Population size="10"/>

    <!-- number of offsprings per generation = rate * population size --> 
    <Offsprings rate="0.25"/>

    <!-- if OnePointHull is present, finds lowest energy state, rather than a true convex hull -->
    <OnePointHull/>

    <!-- if MinimizeBest is present, algorithm minimizes best individuals at end of each step -->
    <!-- rate specifies how many individuals to minimize -->
    <!-- every means minimization occurs every n generations -->
    <MinimizeBest rate="0.1" every="n" />
  </GA>

  <!-- you need to specify Lattice, Clusters and CS in this file as well -->
  <!-- and put it between the LaDa tags -->
  <!-- see liblamarck for description -->

</LaDa>
