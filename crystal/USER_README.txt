The lamarck libraries converts a specified generalized Ising model into a cell-shape specific Polynomial + CS. 
input is described below.

At this juncture, liblamarck does not run of itself. see LaDa for actual use.
to compile, run dependencies.pl 
you will need
   atat installed somewhere on your system (specify location in dependencies.pl)
   tinyxml at ./tinyxml/
as well as libopt++.a and newmat.a libraries, 
by default the above libraries should be compiled with -malign-double
to link properly with LaDa on a 32bit architecture
if not, remove this option from the  makefile

At this point input can be obtained using the following steps:
  (1) CE_to_Atat.pl 
  (2) convert_to_input.xml
Note that the attenuation (in <CS attenuation="?") should be given in
NREL format, not Atat format. There may or may not be a pb with the axis
matrix. When using something other than [2,0,0][0,2,0][0,0,2], please
check that CS is correct. You can use the Lamarck::plot_attenuation
subroutine for this purpose (uncomment it first, and put it
somewhere...).


Each liblamarck class expects a certain xml format.

Below, you will find the input for VA_CE::Functional_Builder
<?xml version="1.0" standalone=no>
<LaDa>
   <!--  specifie the lattice type here -->
   <Lattice>
     <!--  following three indicates cell matrix  -->
     <column x="-0.50000000" y="0.50000000" z="0.50000000"/>
     <column x="0.50000000" y="-0.50000000" z="0.50000000"/>
     <column x="0.50000000" y="0.50000000" z="-0.50000000"/>
     <!--  there should be only one substitutional site  -->
     <site x="0.00000000" y="0.00000000" z="0.00000000">
       <!--  and only two possibilities per site  -->
       <atom type="Cu"/>
       <atom type="Au"/>
     </site>
   </Lattice>
   <!--  clusters (atat parlance, figure in NREL parlance) >
   <Clusters> 
     <!--  interaction energies (eci) must be divided by multiplicity, atat format  -->
     <Cluster eci="-144.69354570">
       <!--  J0, zero body interaction  -->
     </Cluster>
     <Cluster eci="-0.20499451"> 
       <!-- four body interactions, with following list of spins  -->
       <spin x="0.00000000" y="0.00000000" z="0.00000000"/>
       <spin x="1.50000000" y="0.50000000" z="0.50000000"/>
       <spin x="0.50000000" y="0.50000000" z="0.50000000"/>
       <spin x="1.00000000" y="1.00000000" z="0.00000000"/>
     </Cluster>
     ...
   </Clusters>

   <!-- constituent strain: -->
   <!-- input for VA_CE::Constituent_Strain and VA_CE::Harmonic -->
   <!-- attenuation is given in NREL format -->
   <!-- attenuation = 4/ | norm of k|^2  -->
   <!-- If not fcc or bcc, please check that this statement is true  -->
   <CS attenuation="1">
     <!-- each harmonic warrants an xml group of its own -->
     <!-- present implementation goes up to cubic harmonics of rank 3 >
     <Harmonic rank="0"> // harmonic of rank 0
       <!-- coefficients are given such that CubicHarmonic(k) * Coef(x) is an energy -->
       <!-- this is the atat format, as well as the NREL format -->
       <Point x="0.0000" y="0.00000000"/> 
       <Point x="0.0250" y="52.29130627"/>
       ...
       <Point x="1.0000" y="0.00000000"/>
     </Harmonic>
     <Harmonic rank="1">
       ...
     </Harmonic>
     ...
   </CS>
</LaDa>

Input expected by Ising_CE::Structure :: Load (TiXmlElement * )
  
  <!-- always between two Structure tags -->
  <!-- Pi and enegy are optional tagnames for each structure, mostly for debugging -->
  <Structure PI="5687" energy="-131.33108900"> 
    <!-- unit cell is given as matrix -->
    <Cell> 
      <column x="1.50000000" y="0.50000000" z="0.50000000"/>
      <column x="-0.50000000" y="1.50000000" z="0.50000000"/>
      <column x="0.00000000" y="-1.00000000" z="2.00000000"/>
    </Cell> 
    <!-- each site is given in an Atom xmltag -->
    <!-- if needed, the type can be given. it can be any real number -->
    <Atom x="0.00000000"  y="0.00000000"  z="0.00000000" type="1"/>
    ...
  </Structure>


Input expected by VA_CE::Convex_Hull :: Load ( const TiXmlElement * ) 

  <!-- should alway be contained between ConvexHull tags -->
  <ConvexHull>
    <!-- each break point is given separately, with its energy and concentration -->
    <BreakPoint E="0.698904" x="-1.000000">
      <Structure> 
        <!-- see above -->
      </Structure> 
    </BreakPoint>
    ...
  <ConvexHull/>
