#!/usr/bin/perl

# MBCE v0.94
# Mayeul d'Avezac, NREL, 2005
# 
###############################################################3
# Transforms NREL-style CE (fort.11, comp-coef, etc...)  to  input files for
# ATAT's Monte-Carlo codes.
###############################################################3

print "\n  ", $0, " transforms Cluster Expansion to ATAT emc2 input\n\n";

use Cwd;

# sets home directory
my %params;
$params{'home directory'} = `cd ~; pwd`; chomp $params{'home directory'};
$params{'home directory'} .= "/";

###### tunable quantities come here

  # directory containing fort.11 and comp-coef.dat (produced by cefit). 
  # use getcwd() for current directory.
  $params{"input directory"} = getcwd();
  
  # where atat files go
  $params{"output directory"} = getcwd(); #$params{"home directory"} . "Monte_Carlo/epsilon-G/CuAu-PRB35";
 
  # interactions you've used in computing cluster expansion
  $params{"interaction types"} = "3B5.4B4.5B3.6B3";
 
  # PI file you've used to comput the  convex hull
  $params{"PI types"} = "PIs2thru8";

  # at this juncture, should be either bcc or fcc
  $params{"lattice type"} = "bcc"; 
  
  # directory of cluster expansion Jtypes and PI files
  $params{"figure directory"} =   $params{"home directory"} . "Cluster_Expansion/v0.94/lattice/"
                                . $params{"lattice type"} . "-based/figures";

  # Do you want constituent strain? need file comp-coef.dat in input directory
  # put a attenuation ks. empty, or whitespace, means no constituent strain;
  $params{"constituent strain attenuation"} = "1";

  # specify known ground state structures, eg [1,2,3]
  # [] means get structure from $params{"PI file"} . "gslplot.dat"
  @{$params{"structures"}} = (15,1,-15);


  # Hell if I know what this is. Here it will print the alloy constituents in the structure file gs_str.out
  @{$params{"code"}} = ("Au", "Pd");  # ? code, usually [1,2]


  #  below this line, don't touch but rather complain to appropriate authority
  #  below this line, don't touch but rather complain to appropriate authority
  #  below this line, don't touch but rather complain to appropriate authority
  $params{"jtypes file"} = "jtypes." . $params{"interaction types"};
  $params{"PI file"} = $params{"PI types"} . "." . $params{"interaction types"};
  @{$params{"coordinate system"}} = ( [2, 0, 0],    # Coordinate system (multiplication matrix to 
                                      [0, 2, 0],    # get integer atomic positions;
                                      [0, 0, 2] );  # this setting works for fcc and bcc, 
                                                    # change at your own risk)

  # unit cell of the cluster expansion lattice, only bcc and fcc
  if ( $params{'lattice type'} =~ /bcc/ ) 
  {
    @{$params{"lattice unit cell"}} = ( [ -0.5,  0.5,  0.5 ],    
                                        [  0.5, -0.5,  0.5 ],    
                                        [  0.5,  0.5, -0.5 ] );  
  }
  elsif ( $params{'lattice type'} =~ /fcc/ ) 
  {
    @{$params{"lattice unit cell"}} = ( [ 0.0, 0.5, 0.5 ],
                                        [ 0.5, 0.0, 0.5 ],
                                        [ 0.5, 0.5, 0.0 ] );
  }
  else { die "Never heard of lattice type $params{'lattice type'} \n"; }


###### end of tunable quantities


###### begin work #####################
my @pairs;



if (! -e $params{"output directory"}) 
{
  print $params{"output directory"}," does not exist! Creating ...\n";
  system "mkdir", $params{"output directory"};
}


print "\t_ writing $params{'output directory'}/gs_str.out\n";
write_gs_str();

print "\t_ writing $params{'output directory'}/lat.in\n";
write_lat();

print "\t_ writing $params{'output directory'}/eci.out and \n\t and $params{'output directory'}/clusters.out\n";
write_eci_and_clusters();

if ( $params{'constituent strain attenuation'} =~ /\S+/ ) 
{
  print "\t_ writing $params{'output directory'}/cs.in \n";
  write_cs();
}

print "\t_ writing $params{'output directory'}/wherefrom \n";
write_wherefrom();




sub write_gs_str(){

  $wherefromstruct = 1;
  if (     ( scalar( @{$params{'structures'}} ) == 0 ) 
       and ( -e "$params{'input directory'}/$params{'PI file'}.gslplot.data") ) 
  {
    $wherefromstruct = 0;
    open (STR, "$params{'input directory'}/$params{'PI file'}.gslplot.data")
       or die "failed to open $params{'input directory'}/$params{'PI file'}.gslplot.data";

    while( $_=<STR> ) {
      /\s+(\S+)+\s+/;
      if ( $1 != 0 ) { push @{$params{"structures"}}, $1 ; }
    }
    close(STR);
  }
  elsif ( ! ( -e "$params{'input directory'}/$params{'PI file'}.gslplot.data" ) )
  { 
    printf " WARNING -- no structures on input, no file %s/%s.gslplot.data\n",
                $params{'input directory'}, $params{'PI file'}; 
  };
   
  open(STR,">$params{'output directory'}/gs_str.out") or  die "Couldn't write to gs_str.out";

  # prints x=0 structure
  my $cell = \@{$params{'coordinate system'}};
  printf STR "\t%12.8f %12.8f %12.8f\n\t%12.8f %12.8f %12.8f\n\t%12.8f %12.8f %12.8f\n", 
             $cell->[0][0], $cell->[0][1], $cell->[0][2],
             $cell->[1][0], $cell->[1][1], $cell->[1][2], 
             $cell->[2][0], $cell->[2][1], $cell->[2][2]; 
  my $cell = \@{$params{'lattice unit cell'}};
  printf STR "%12.8f %12.8f %12.8f\n%12.8f %12.8f %12.8f\n%12.8f %12.8f %12.8f\n", 
             $cell->[0][0], $cell->[0][1], $cell->[0][2],
             $cell->[1][0], $cell->[1][1], $cell->[1][2], 
             $cell->[2][0], $cell->[2][1], $cell->[2][2]; 
  printf STR " %10.6f  %10.6f  %10.6f   %2s \n",0,0,0, $params{"code"}[0];
  print STR "end\n\n";



  open PI, "$params{'figure directory'}/$params{'PI file'}" 
       or die " no $params{'figure directory'}/$params{'PI file'} \n" ;

  my $foundstructure = 0;
  my $invertme = 0;

  # begin work
  for( $str=0; $str < scalar( @{$params{'structures'}} ); $str++) {

    $foundstructure = 0;

    # prints coordinate basis to file -- CE in orthogonal basis with a=2
    my $cell = \@{$params{'coordinate system'}};
    print STR "\t", $cell->[0][0], " ",  $cell->[0][1], " ",  $cell->[0][2], "\n";
    print STR "\t", $cell->[1][0], " ",  $cell->[1][1], " ",  $cell->[1][2], "\n";
    print STR "\t", $cell->[2][0], " ",  $cell->[2][1], " ",  $cell->[2][2], "\n";
    if ($params{'structures'}[$str] < 0 ) {
      $invertme = 1; $params{'structures'}[$str] *= -1 }
    else {
      $invertme = 0; }

    # prints unit cell basis to file, respect to basis above
    while(<PI>){
        $i = $params{'structures'}[$str];

        if(/NO\. *$i /){
            $foundstructure = 1; # true !
            ($dummy,$istruct,$minority,$Nsites,
             $binary,$dummy,@a1[1..3],@a2[1..3],@a3[1..3]) = split;
            for($i=1;$i<=3;$i++){ # Make lattice constant 1 (not 2)
                $a1[$i] /= 2;$a2[$i] /= 2;$a3[$i] /= 2;}
            printf STR "%12.8f %12.8f %12.8f\n", @a1[1..3];
            printf STR "%12.8f %12.8f %12.8f\n", @a2[1..3];
            printf STR "%12.8f %12.8f %12.8f\n", @a3[1..3];
            while(<PI>){if(/BASIS/){last;}} # Skip the PIs
            s/BASIS//; # Get the basis atom positions
            @basis = split;
            $iatom = 0;
            while(@basis){
                printf STR " %10.6f  %10.6f  %10.6f   %2s \n",
                ,$basis[0]/2,$basis[1]/2,$basis[2]/2, $params{"code"}[abs($invertme-($binary>>$iatom)%2)];
                shift @basis;shift @basis;shift @basis; # Get next vector
                $iatom++;
                if($iatom == $Nsites){last;} # Break out if all atoms printed
                if(!($iatom%8)) {$_=<PI>;@basis = split;}
            }

            if ($Nsites != $iatom){die "Error";}
            print STR "end\n\n";
            seek (PI, 0, 2);
        }
    }
    if ( $foundstructure==0 ) { die " couldn't find structure $params{'structures'}[$str] in PI file! ";}
    if ( $invertme==1) { $params{'structures'}[$str] *= -1; }
    seek (PI,1,0);
  }

  close PI;

  # prints x=1 structure
  my $cell = \@{$params{'coordinate system'}};
  printf STR "\t%12.8f %12.8f %12.8f\n\t%12.8f %12.8f %12.8f\n\t%12.8f %12.8f %12.8f\n", 
             $cell->[0][0], $cell->[0][1], $cell->[0][2],
             $cell->[1][0], $cell->[1][1], $cell->[1][2], 
             $cell->[2][0], $cell->[2][1], $cell->[2][2]; 
  my $cell = \@{$params{'lattice unit cell'}};
  printf STR "%12.8f %12.8f %12.8f\n%12.8f %12.8f %12.8f\n%12.8f %12.8f %12.8f\n", 
             $cell->[0][0], $cell->[0][1], $cell->[0][2],
             $cell->[1][0], $cell->[1][1], $cell->[1][2], 
             $cell->[2][0], $cell->[2][1], $cell->[2][2]; 
  printf STR " %10.6f  %10.6f  %10.6f   %2s \n",0,0,0, $params{"code"}[1];
  print STR "end\n\n";

  close STR;
}

sub write_lat(){

  open(STR, ">$params{'output directory'}/lat.in") or die;

  my $cell = \@{$params{'coordinate system'}};
  printf STR "\t%12.8f %12.8f %12.8f\n\t%12.8f %12.8f %12.8f\n\t%12.8f %12.8f %12.8f\n", 
             $cell->[0][0], $cell->[0][1], $cell->[0][2],
             $cell->[1][0], $cell->[1][1], $cell->[1][2], 
             $cell->[2][0], $cell->[2][1], $cell->[2][2]; 
  my $cell = \@{$params{'lattice unit cell'}};
  printf STR "%12.8f %12.8f %12.8f\n%12.8f %12.8f %12.8f\n%12.8f %12.8f %12.8f\n", 
             $cell->[0][0], $cell->[0][1], $cell->[0][2],
             $cell->[1][0], $cell->[1][1], $cell->[1][2], 
             $cell->[2][0], $cell->[2][1], $cell->[2][2]; 
  print STR "0.0 0.0 0.0 ", $params{"code"}[0], ", ", $params{"code"}[1], "\n"; 


  close (STR);
  
}


sub write_eci_and_clusters(){

  open (FORT, "$params{'input directory'}/fort.11") or die;
  
  # finds number of non-pairs 
  my $Nmb = <FORT>; chop $Nmb;

  # reads Labels of non-pairs
  my $i=0, $j,  $k, $l;
  my @Labels;
  while($i < $Nmb ){
    $_ = <FORT>;
    @Labels = (@Labels, split);
    $i = $#Labels + 1;
  }

  # reads ECIs of non-pairs
  $i=0; my @ECIs;
  while($i < $Nmb ){
    $_ = <FORT>; @ECIs = (@ECIs, split); $i = $#ECIs + 1;
  }

  # find number of pairs
  my $Nb_pairs = <FORT>; chop $Nb_pairs;

  # read ECIs of pairs
  while($_ = <FORT>) {
    @ECIs = (@ECIs, split); $i = $#ECIs;
  }

  close(FORT);

  # printing many body clusters
  open(JTYPE, "$params{'figure directory'}/$params{'jtypes file'}") 
      or die "Couldn't open $params{'figure directory'}/$params{'jtypes file'}\n" ;
  open (STR, ">$params{'output directory'}/clusters.out")
      or die "Couldn't open $params{'output directory'}/clusters.out \n";
  open (ECI, ">$params{'output directory'}/eci.out")
      or die "Couldn't open $params{'output directory'}/eci.out \n";
  for ($i=0; $i < $#Labels+1; $i++) {

    my $label = $Labels[$i];

    while( ($_=<JTYPE>) && ($_!~/$label/) ){}
    $_ = <JTYPE>;
    my $n, $mult, $dist;
    ($n, $mult, $dist) = (split);
    if ($mult == 0 ) { $mult = 1};
    print STR " $mult\n $dist\n $n\n";
    for ($j=0; $j < $n; $j++) {
      $_ = <JTYPE>;
      my @strng = (split);
      printf STR " %6.3f %6.3f %6.3f \n", $strng[0]/2.0, $strng[1]/2.0, $strng[2]/2.0;
    }
    print STR "\n";

    # print ECI divided by multiplicity
    $_=$ECIs[$i]; /(\S+)/;
    print ECI $1/$mult, "\n";
    
    
    seek(JTYPE, 0,0);
  }
  close(JTYPE);

  # pair order
  my $lattice = \@{$params{'lattice unit cell'}};
  get_pairs();

  # prints pairs
  for($i=0; $i < $Nb_pairs; $i++) {
    print STR " 0.0\n 0.0\n 2\n ";
    print STR "0 0 0 \n ";
    printf STR " %.5f %.5f %.5f \n\n", 
               $pairs[$i][0]*0.5, 
               $pairs[$i][1]*0.5, 
               $pairs[$i][2]*0.5;

    $_=$ECIs[$i+$Nmb]; /(\S+)/;
    print ECI $1/$pairs[$i][3], "\n";
  }

  close(STR);
  close(ECI);
}

sub write_cs(){

  if ( ! -e $params{'input directory'}."/comp\-coef.dat" ) {
    print "Couldn't find ".$params{'input directory'}."/comp\-coef.dat \n No constituent strain added";
    $params{'constituent strain'} = "none";
    return; 
  }
  open IN, $params{'input directory'}."/comp\-coef.dat"  
     or die "Couldn't open ".$params{'input directory'}."/comp\-coef.dat \n";

  my @strain;
  while ( ($_ = <IN>) ) 
  {
    if ( /\S+/ ) 
    {
      my @line, $i=0;
      do {
        /(\S+)/; $line[$i] = $1; $i++;
      } while ( ($_=<IN>) and (/\S+/) );
      push @strain, ( [ @line[1..$#line] ] );
    }
  }

  close IN;
  
  open OUT, ">$params{'output directory'}/cs.in" 
    or die "Couldn't write to $params{'output directory'}/cs.in \n";

  printf OUT " %f  %d  %d \n", $params{"constituent strain attenuation"}, 
         scalar( @{$strain[0]} ), scalar( @strain );
  for ( my $i=0; $i < scalar( @{$strain[0]} ); $i++ )
  {
    for ( my $j=0; $j < scalar( @strain ); $j++ )
    {
      printf OUT " %.8f \n", $strain[$j][$i];
    }
    printf OUT "\n";
  }

  close OUT;
}


sub write_wherefrom(){

  open (STR, ">$params{'output directory'}/wherefrom" ) or die;
  print STR "data in lat.in, gs_str.out, eci.out. clusters.out, cs.in obtained from \n";
  print STR "\t_ $params{'input directory'} files \n";
  print STR "\t_ $params{'figure directory'}/$params{'PI file'} \n";
  print STR "\t_ $params{'figure directory'}/$params{'interaction types'} \n";
  print STR "\t_ Structures: ";
  if (!$wherefromstruct) {
    print STR "in $params{'input directory'}/$params{'PI types'}.gslplot.data \n\t\t"; }
  else {
    print STR "given on input\n\t\t "; }
  my $i;
  for ($i=0; $i < scalar( @{$params{'structures'}} )-1; $i++) {
    print STR $params{'structures'}[$i], ", ";
  } 
  print STR $params{'structures'}[$i], "\n"; 
  print STR "\t_ code ";
  for ($i=0; $i < scalar( @{$params{"code"}} )-1 ; $i++) {
    print STR $params{"code"}[$i], ", ";
  } 
  print STR $params{"code"}[$i], "\n"; 
  print STR "\t_ Coordinate System: \n"; 
  my $Coor = \@{$params{'coordinate system'}};
  print STR "\t\t", $Coor->[0][0], " ",  $Coor->[0][1], " ",  $Coor->[0][2], "\n";
  print STR "\t\t", $Coor->[1][0], " ",  $Coor->[1][1], " ",  $Coor->[1][2], "\n";
  print STR "\t\t", $Coor->[2][0], " ",  $Coor->[2][1], " ",  $Coor->[2][2], "\n";
  print STR "\t_ Real space lattice vectors: \n"; 
  my $lattice = \@{$params{'lattice unit cell'}};
  print STR "\t\t", $lattice->[0][0], " ",  $lattice->[0][1], " ",  $lattice->[0][2], "\n";
  print STR "\t\t", $lattice->[1][0], " ",  $lattice->[1][1], " ",  $lattice->[1][2], "\n";
  print STR "\t\t", $lattice->[2][0], " ",  $lattice->[2][1], " ",  $lattice->[2][2], "\n";
  close (STR);

}


sub get_pairs()
{
  if ($params{'lattice type'} =~ /fcc/ ) 
  {
    push @pairs, ([1, 1, 0,  6]);
    push @pairs, ([2, 0, 0,  3]);
    push @pairs, ([2, 1, 1, 12]);
    push @pairs, ([2, 2, 0,  6]);
    push @pairs, ([3, 1, 0, 12]);
    push @pairs, ([2, 2, 2,  4]);
    push @pairs, ([3, 2, 1, 24]);
    push @pairs, ([4, 0, 0,  3]);
    push @pairs, ([3, 3, 0,  6]);
    push @pairs, ([4, 1, 1, 12]);
    push @pairs, ([4, 2, 0, 12]);
    push @pairs, ([3, 3, 2, 12]);
    push @pairs, ([4, 2, 2, 12]);
    push @pairs, ([4, 3, 1, 24]);
    push @pairs, ([5, 1, 0, 12]);
    push @pairs, ([5, 2, 1, 24]);
    push @pairs, ([4, 4, 0,  6]);
    push @pairs, ([4, 3, 3, 12]);
    push @pairs, ([5, 3, 0, 12]);
    push @pairs, ([4, 4, 2, 12]);
    push @pairs, ([6, 0, 0,  3]);
    push @pairs, ([5, 3, 2, 24]);
    push @pairs, ([6, 1, 1, 12]);
    push @pairs, ([6, 2, 0, 12]);
    push @pairs, ([5, 4, 1, 24]);
    push @pairs, ([6, 2, 2, 12]);
    push @pairs, ([6, 3, 1, 24]);
    push @pairs, ([4, 4, 4,  4]);
    push @pairs, ([5, 4, 3, 24]);
    push @pairs, ([5, 5, 0,  6]);
    push @pairs, ([7, 1, 0, 12]);
    push @pairs, ([6, 4, 0, 12]);
    push @pairs, ([5, 5, 2, 12]);
    push @pairs, ([6, 3, 3, 12]);
    push @pairs, ([7, 2, 1, 24]);
    push @pairs, ([6, 4, 2, 24]);
    push @pairs, ([7, 3, 0, 12]);
    push @pairs, ([6, 5, 1, 24]);
    push @pairs, ([7, 3, 2, 24]);
    push @pairs, ([8, 0, 0,  3]);
    push @pairs, ([5, 5, 4, 12]);
    push @pairs, ([7, 4, 1, 24]);
    push @pairs, ([8, 1, 1, 12]);
    push @pairs, ([6, 4, 4, 12]);
    push @pairs, ([8, 2, 0, 12]);
    push @pairs, ([6, 5, 3, 24]);
    push @pairs, ([6, 6, 0,  6]);
    push @pairs, ([8, 2, 2, 12]);
    push @pairs, ([7, 4, 3, 24]);
    push @pairs, ([7, 5, 0, 12]);
  } 
  elsif ($params{'lattice type'} =~ /bcc/ ) 
  {
    push @pairs, ([ 1, 1, 1,  4]);
    push @pairs, ([ 2, 0, 0,  3]);
    push @pairs, ([ 2, 2, 0,  6]);
    push @pairs, ([ 3, 1, 1, 12]);
    push @pairs, ([ 2, 2, 2,  4]);
    push @pairs, ([ 4, 0, 0,  3]);
    push @pairs, ([ 3, 3, 1, 12]);
    push @pairs, ([ 4, 2, 0, 12]);
    push @pairs, ([ 4, 2, 2, 12]);
    push @pairs, ([ 3, 3, 3,  4]);
    push @pairs, ([ 5, 1, 1, 12]);
    push @pairs, ([ 4, 4, 0,  6]);
    push @pairs, ([ 5, 3, 1, 24]);
    push @pairs, ([ 4, 4, 2, 12]);
    push @pairs, ([ 6, 0, 0,  3]);
    push @pairs, ([ 6, 2, 0, 12]);
    push @pairs, ([ 5, 3, 3, 12]);
    push @pairs, ([ 6, 2, 2, 12]);
    push @pairs, ([ 4, 4, 4,  4]);
    push @pairs, ([ 5, 5, 1, 12]);
    push @pairs, ([ 7, 1, 1, 12]);
    push @pairs, ([ 6, 4, 0, 12]);
    push @pairs, ([ 6, 4, 2, 24]);
    push @pairs, ([ 5, 5, 3, 12]);
    push @pairs, ([ 7, 3, 1, 24]);
    push @pairs, ([ 8, 0, 0,  3]);
    push @pairs, ([ 7, 3, 3, 12]);
    push @pairs, ([ 6, 4, 4, 12]);
    push @pairs, ([ 8, 2, 0, 12]);
    push @pairs, ([ 6, 6, 0,  6]);
    push @pairs, ([ 8, 2, 2, 12]);
    push @pairs, ([ 5, 5, 5,  4]);
    push @pairs, ([ 7, 5, 1, 24]);
    push @pairs, ([ 6, 6, 2, 12]);
    push @pairs, ([ 8, 4, 0, 12]);
    push @pairs, ([ 7, 5, 3, 24]);
    push @pairs, ([ 9, 1, 1, 12]);
    push @pairs, ([ 8, 4, 2, 24]);
    push @pairs, ([ 6, 6, 4, 12]);
    push @pairs, ([ 9, 3, 1, 24]);
    push @pairs, ([ 8, 4, 4, 12]);
    push @pairs, ([ 7, 5, 5, 12]);
    push @pairs, ([ 7, 7, 1, 12]);
    push @pairs, ([ 9, 3, 3, 12]);
    push @pairs, ([ 8, 6, 0, 12]);
    push @pairs, ([10, 0, 0,  3]);
    push @pairs, ([ 8, 6, 2, 24]);
    push @pairs, ([10, 2, 0, 12]);
    push @pairs, ([ 7, 7, 3, 12]);
    push @pairs, ([ 9, 5, 1, 24]);
  } 
}



sub host
  {

    $params{'home directory'} = `cd ~; pwd`; chomp $params{'home directory'};
    $params{'home directory'} .= "/";
    
    $_ = `hostname` ;
    if (/linux/) 
      {
        $params{'figure directory'} = $params{'home directory'} 
                                     ."Cluster_Expansion/v0.94/lattice/fcc-based/figures" ;
      }
    elsif (/master/) 
      {
        $params{'figure directory'} = $params{'home directory'} 
                                     ."Cluster_Expansion/v0.94/lattice/fcc-based/figures" ;
      }
    else
      {
        print `hostname`, " is an unknown host - please specify home directory in subroutine host!\n" ;
        die ;
      }

  } 

