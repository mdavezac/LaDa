#!/usr/bin/perl

my %params;
open IN, "COMMANDS";
my @COMMANDS;
push @COMMANDS, ( <IN> );
push @COMMANDS, ( <IN> );
push @COMMANDS, ( <IN> );
push @COMMANDS, ( <IN> );
push @COMMANDS, ( <IN> );
push @COMMANDS, ( <IN> );
close IN;

$params{"home"} = `cd; pwd`; chomp $params{"home"};
$params{'directory'}{'result'} =  $COMMANDS[0]; chomp $params{'directory'}{'result'};
$params{'directory'}{'work'}   =  $COMMANDS[1]; chomp $params{'directory'}{'work'};
$params{'GA style'}            =  $COMMANDS[2]; chomp $params{'GA style'};
$params{'GA style'}            =~ s/^\s+//;
$params{'GA style'}            =~ s/\s+$//;
$params{'minimizer'}           =  $COMMANDS[3]; chomp $params{'minimizer'};
$params{'minimizer'}            =~ s/\s+$//;
$params{'minimizer'}            =~ s/^\s+//;
$params{'CH'}                  =  $COMMANDS[4]; chomp $params{'CH'};
$_ = $COMMANDS[5]; /(\d+)/; $params{'max calls'} = $1;
$_ = $COMMANDS[6]; /generations:\s+(\d+)\s+replacement:\s+(\S+)/;
$params{"GA"}{"population"} = $1;
$params{"GA"}{"replace per generation"} = $2;
$params{'iaga call'}           =  "lada > out";

$params{'nb atoms'} = 20;
$params{"file"}{"Pi"} = "$params{'home'}/nanopse/cell_shapes/fcc_7-32";
$params{'max GA iters'} = 1000;



if ( $params{'CH'} =~ /one point/i )
{
  $params{'GA style'} = "one_point_$params{'GA style'}";
}
if ( $params{'GA style'} =~ /true/ )
{ 
  $params{'GA style'} =~ s/true//; 
  $params{'GA style'}            =~ s/^\s+//;
  $params{'GA style'}            =~ s/\s+$//;
  $params{'GA style'} = "true_$params{'GA style'}"; chomp $params{'GA style'};
}



$params{'agr'}{"filename"} = "$params{'GA style'}_$params{'minimizer'}";
$params{'xml'}{"filename"} = "$params{'GA style'}_$params{'minimizer'}";

if ( $params{'max calls'} != 0 )
{
  $params{'agr'}{"filename"} = "$params{'GA style'}_$params{'minimizer'}_n=$params{'max calls'}";
  $params{'xml'}{"filename"} = "$params{'GA style'}_$params{'minimizer'}_n=$params{'max calls'}";
}

if ( $params{"GA"}{"population"} != 100 )
{
  $params{'agr'}{"filename"} .= "_gen:$params{'GA'}{'population'}";
  $params{'xml'}{"filename"} .= "_gen:$params{'GA'}{'population'}";
}
if ( $params{"GA"}{"replace per generation"} != 0.1 )
{
  $params{'agr'}{"filename"} .= "_rep:$params{'GA'}{'replace per generation'}";
  $params{'xml'}{"filename"} .= "_rep:$params{'GA'}{'replace per generation'}";
}

$params{'agr'}{"filename"} = "$params{'agr'}{'filename'}.agr";
$params{'xml'}{"filename"} = "$params{'xml'}{'filename'}.xml";


# begin work
my %structure;

open (PI, "$params{'file'}{'Pi'}") 
  or die "Did not find PI file $params{'file'}{'Pi'}";

while(<PI>)
{
  if(/NT/)
  {
    read_structure(); # data passes from PIfile to %structure hash
	
    if ( scalar( @{$structure{'atomic positions'}} ) < $params{'nb atoms'}*3 ) 
      { next; }
    else 
      { launch_iaga(); };

  }
}

close PI;



exit;


sub launch_iaga()
{
  system "rm -f $params{'agr'}{'filename'}";
  for( my $count =0; $count < $params{'max GA iters'}; $count++)
  {
    write_lamarck_input();
    if ( $params{'GA style'} !~ /true/ )
      { system "rm -f convex_hull.xml"; }
    system "$params{'iaga call'}";
    system "echo '# new GA ' >> $params{'directory'}{'result'}/$params{'agr'}{'filename'}";
    system "cat convex_hull.agr >> $params{'directory'}{'result'}/$params{'agr'}{'filename'}";
    system "cp convex_hull.xml $params{'directory'}{'result'}/$params{'xml'}{'filename'}  ";
  }
exit;
}

sub read_structure()
{
  my $dummy;
  ($dummy,$structure{'nb sites'},$structure{'nb triads'},
   @{$structure{'cell'}[0]}[1..3],@{$structure{'cell'}[1]}[1..3],
   @{$structure{'cell'}[2]}[1..3]) = (split);
  
  for($i=1;$i<=3;$i++)
  { 
    $structure{'cell'}[0][$i] /= 2;
    $structure{'cell'}[1][$i] /= 2;
    $structure{'cell'}[2][$i] /= 2;
  }

  # Skip the PIs
  while( ($_=<PI>) and ($_ !~ /BASIS/) ){}

  # atoms
  s/BASIS//; # Get the basis atom positions
  @{$structure{'atomic positions'}} = (split);
}


sub write_lamarck_input()
{
  open IN, "<input_template.xml"
    or die "No input_template.xml file found";
  open OUT, ">input.xml"
    or die "Could not open input.xml";

  while( <IN> )
  {

    if( /\?structure\?/ )
    {
      printf OUT "  <Structure>\n";
      printf OUT "    <Cell>\n";
      printf OUT "      <column x=\"%.6f\" y=\"%.6f\" z=\"%.6f\" />\n",
                 @{$structure{'cell'}[0]}[1..3];
      printf OUT "      <column x=\"%.6f\" y=\"%.6f\" z=\"%.6f\" />\n",
                 @{$structure{'cell'}[1]}[1..3];
      printf OUT "      <column x=\"%.6f\" y=\"%.6f\" z=\"%.6f\" />\n",
                 @{$structure{'cell'}[2]}[1..3];
      printf OUT "    </Cell>\n";

      for( my $i = 0; $i < scalar( @{$structure{'atomic positions'}} ); 
           $i +=3 )
      {
        printf OUT "    <Atom x=\"%.6f\" y=\"%.6f\" z=\"%.6f\" type=\"%.6f\" />\n",
                 $structure{'atomic positions'}[$i]/2, 
                 $structure{'atomic positions'}[$i+1]/2, 
                 $structure{'atomic positions'}[$i+2]/2; 
      }
      printf OUT "  </Structure>\n";
    }
    elsif( /\?GA\?/ )
    {
      printf OUT "  <GA";
      if ( $params{'GA style'} =~ /lamarck_evolution/i )
        { printf OUT " method=\"lamarck evolve from start\""; }
      elsif ( $params{'GA style'} =~ /lamarck/i )
        { printf OUT " method=\"lamarck\""; }
      elsif ( $params{'GA style'} =~ /multistart/i )
        { printf OUT " method=\"lamarck multistart evolve from start\""; }
      elsif ( $params{'GA style'} =~ /darwin/i )
        { printf OUT " method=\"darwin\""; }
      printf OUT " maxgen=\"%i\">\n",
                 $params{'GA'}{'max generations'};

      if ( $params{'minimizer'} =~ /linear/i )
      {
        printf OUT "    <Minimizer type=\"linear\" maxeval=\"%i\" />\n",
                   $params{'max calls'};
      }
      elsif ( $params{'minimizer'} =~ /sa/i )
      {
        printf OUT "    <Minimizer type=\"SA\" maxeval=\"%i\" />\n",
                   $params{'max calls'};
      }
      elsif ( $params{'minimizer'} =~ /wang/i )
        { printf OUT "    <Minimizer type=\"wang\"/>\n"; }
      elsif ( $params{'minimizer'} =~ /physical/i )
        { printf OUT "    <Minimizer type=\"physical\"/>\n"; }

      printf OUT "    <Operators type=\"or\" prob=\"0.75\" >\n";
      printf OUT "      <Mutation prob=1.0 />\n";
      printf OUT "      <Crossover prob=0.5 />\n";
      printf OUT "    </Operators>\n";

      printf OUT "    <Selection size=2 />\n";

      printf OUT "    <Population size=\"%i\"/>\n",
                 $params{'GA'}{'population'};
      printf OUT "    <Offspings rate=\"%i\"/>\n",
                 $params{'GA'}{'replace per generation'};
      if( $params{'CH'} =~ /one point/i )
        { printf OUT "    <OnePointHull/>\n"; }
      printf OUT "  </GA>\n";
    }
    else
    {
      print OUT $_;
    }
  }
  

  close IN;
  close OUT;
}

