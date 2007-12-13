#! /usr/bin/perl
#
#  Version: $Id$
#

my %escan;
my %vff;
my %structure;

my $arguments;
for( my $i = 0; $i <= $#ARGV; $i++ )
  { $argument = sprintf "%s %s ", $argument, $ARGV[$i]; }

if( $argument =~ /(--help|-h)/ )
{
  print "help:\n --help, -h this message.\n";
  print "All options are of the form --option filename,\n";
  print "where filename is the file to convert to XML\n";
  print "--escan, -e converts nanopse escan input file to xml LaDa input\n";
  print "--genpot, -g converts nanopse genpot input to xml LaDa input\n";
  print "--vff, -v converts nanopse vff input to xml LaDa input\n";
  print "--position, -p converts nanopse vff position input to xml LaDa input\n";
}

if( $argument =~ /(--escan|-e)\s+(\S+)/ )
{
  my $inputfile = $2;
  read_escan($inputfile);
}
if( $argument =~ /(--genpot|-g)\s+(\S+)/ )
{
  my $inputfile = $2;
  read_potential($inputfile);
}
write_escan($inputfile)
  if(    $argument =~ /(--escan|-e)\s+(\S+)/ 
      or $argument =~ /(--genpot|-g)\s+(\S+)/ );

if( $argument =~ /(--vff|-v)\s+(\S+)/ )
{
  my $inputfile = $2;
  read_vff($inputfile);
  write_vff();
}

if( $argument =~ /(--position|-p)\s+(\S+)/ )
{
  my $inputfile = $2;
  read_structure($inputfile);
  write_structure();
}

sub read_potential(\$)
{
  my $filename = $_[0];
  open IN, "$filename" or die "Couldn't open file $filename\n";

  $_=<IN>; /(\S+)/; $escan{"atomic input file"} = $4;
  $_=<IN>; /(\d+)\s+(\d+)\s+(\d+)/; @{$escan{"pot mesh"}} = ($1, $2, $3);
  <IN>; <IN>; <IN>; <IN>;
  while( ($_=<IN>) and /(\S+)/ )
    { push @{$escan{"pseudos"}}, $1; }
}

sub read_escan(\$)
{
  my $filename = $_[0];
  open IN, "$filename" or die "Couldn't open file $filename\n";

  $_=<IN>; 
  /^(\s+|)(\d|)(\s+|)(\S+)/; $escan{"potential input file"} = $4;
  $_=<IN>;
  /^(\s+|)(\d|)(\s+|)(\S+)/; $escan{"wavefunction output file"} = $4;
  $_=<IN>;
  /^(\s+|)(\d|)(\s+|)(\d)/; $escan{"Solvation Method"} = $4;
  $_=<IN>;
  /^(\s+|)(\d|)(\s+|)(\S+)(,|)(\s+|)(\S+)(,|)(\s+|)(\S+)(,|)(\s+|)(\S+)/; 
  $escan{"Reference Energy"} = $4;
  $escan{"cutoff"} = $7;
  $escan{"Smooth"} = $10;
  $escan{"Kinetic scaling"} = $13;
  $_=<IN>;
  /^(\s+|)(\d|)(\s+|)(\d+)/; $escan{"nbstates"} = $4;
  $_=<IN>;
  /^(\s+|)(\d|)(\s+|)(\d+)(,|)(\s+|)(\d+)(,|)(\s+|)(\S+)/; 
  $escan{"niter"} = $4;
  $escan{"nline"} = $7;
  $escan{"tolerance"} = $10;
  <IN>; <IN>; $_=<IN>;
  /^(\s+|)(\d|)(\s+|)(\S+)/; $escan{"wavefunction input file"} = $4;
  <IN>; $_=<IN>;
  /^(\s+|)(\d\d|)(\s+|)(\d+)/;
  if( $4 == 0 )
  {
    $escan{"kpoint"}[0] = 0;
    $escan{"kpoint"}[1] = 0;
    $escan{"kpoint"}[2] = 0;
  }
  else
  {
    /^(\s+|)(\d\d|)(\s+|)(\d+)\s+(\S+)\s+(\S+)\s+(\s+)/;
    $escan{"kpoint"}[0] = $6;
    $escan{"kpoint"}[1] = $8;
    $escan{"kpoint"}[2] = $10;
  }
  $_=<IN>;
  /^(\s+|)(\d\d|)(\s+|)(\d)/; $escan{"Hamiltonian"} = $4;
  $_=<IN>; 
  /^(\s+|)(\d\d|)(\s+|)(\S+)/; $escan{"atomic input file"} = $4;
  $_=<IN>;
  /^(\s+|)(\d\d|)(\s+|)(\S+)/; $escan{"real space cutoff"} = $4;
  $_=<IN>;
  /^(\s+|)(\d\d|)(\s+|)(\d+)/; $escan{"number of atom types"} = $4;
  while( ($_=<IN>) )
  {
    my $line = $_;
    $line =~ s/\,//g;
    $line =~ /^(\s+|)(\d\d|)(\s+|)(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
    my $name = $4; my $type = $5;
    push @{$escan{"potentials"}{$name}}, ([ $type, $6, $7, $8, $9, $10 ]);
  }

  close IN;
}

sub write_escan ()
{
  printf "  <Functional type=\"escan\" method=\"%i\" >\n",
         $escan{"Solvation Method"};
  printf "    <Reference value=\"%8.4f\" />\n",
         $escan{"Reference Energy"};
  printf "    <Minimizer niter=\"%i\" nlines=\"%i\" tolerance=\"%2.1e\" />\n",
         $escan{"niter"}, $escan{"nline"}, $escan{"tolerance"};
  printf "    <Wavefunctions in=\"%s\" out=\"%s\" />\n",
         $escan{"wavefunction input file"}, $escan{"wavefunction output file"};
  printf "    <Hamiltonian nbstates=\"%i\" smooth=\"%8.4f\" kinscal=\"%6.3f\" ",
         $escan{"nbstates"}, $escan{"Smooth"}, $escan{"Kinetic scaling"};
  printf " potential=\"%i\" realcutoff=\"%8.4f\" launch=\"~/usr/bin/escanCNL\" >\n",
         $escan{"Hamiltonian"}, $escan{"real space cutoff"}, $escan{"launch"};
  foreach $name ( keys %{$escan{"potentials"}} )
  {
    for( my $i = 0; $i < scalar( @{$escan{"potentials"}{$name}} ); $i++ )
    {
      printf "      <SpinOrbit filename=\"%s\" izz=\"%5i\" ",
             $name, $escan{"potentials"}{$name}[$i][0];
      printf "s=\"%8f\" ", $escan{"potentials"}{$name}[$i][1]
        if( abs( $escan{"potentials"}{$name}[$i][1] ) > 0.00001 );
      printf "p=\"%8f\" ", $escan{"potentials"}{$name}[$i][2]
        if( abs( $escan{"potentials"}{$name}[$i][2] ) > 0.00001 );
      printf "d=\"%8f\" ", $escan{"potentials"}{$name}[$i][3]
        if( abs( $escan{"potentials"}{$name}[$i][3] ) > 0.00001 );
      printf "pnl=\"%8f\" ", $escan{"potentials"}{$name}[$i][4]
        if( abs( $escan{"potentials"}{$name}[$i][4] ) > 0.00001 );
      printf "dnl=\"%8f\" ", $escan{"potentials"}{$name}[$i][5]
        if( abs( $escan{"potentials"}{$name}[$i][5] ) > 0.00001 );
      print "/>\n";
    }
  }
  print "    </Hamiltonian>\n";
  printf "    <GenPot x=\"%i\" y=\"%i\" z=\"%i\" cutoff=\"%4.2f\" launch=\"~/usr/bin/getVLarg\" >\n",
         $escan{"pot mesh"}[0], $escan{"pot mesh"}[1], $escan{"pot mesh"}[2], $escan{"cutoff"};
  for( my $i=0; $i < scalar( @{$escan{"pseudos"}} ); $i++ )
    { printf "      <Pseudo filename=\"%s\" />\n", $escan{"pseudos"}[$i]; }
  print "    </GenPot>\n";
  print "  </Functional>\n";
}

sub read_vff(\$)
{
  my $input = $_[0];
  open IN, "$input" or die " Couldn't open $input\n";

  while( ($_=<IN>) )
  {
    if( /bond-length/ )
    {
      <IN>;
      while( ($_=<IN>) and ($_=~/(\S+)\s+(\S+)/ ) )
      {
        my $A = $1, $B = $2;
        my $line = $_;
        $line =~ s/$A//;
        $line =~ s/$B//;
        $line =~ s/^\s+//;
        @{$vff{"bonds"}{$A}{$B}} = (split /\s+/, $line );
      }
    }
    elsif( /bond-angle/ )
    {
      <IN>;
      while( ($_=<IN>) and ($_=~/(\S+)\s+(\S+)+\s+(\S+)/ ) )
      {
        my $A = $1, $B = $2, $C = $3;
        my $line = $_;
        $line =~ s/$A//;
        $line =~ s/$B//;
        $line =~ s/$C//;
        $line =~ s/^\s+//;
        @{$vff{"angles"}{$A}{$B}{$C}} = (split /\s+/, $line );
      }
    }
  } 
  close IN
}

sub write_vff()
{
  print "  <Functional type=\"vff\">\n";
  foreach $A ( keys %{$vff{"bonds"}} )
  {
    foreach $B ( keys %{$vff{"bonds"}{$A}} )
    {
      printf "    <Bond A=\"%s\" B=\"%s\" d0=\"%6.4f\" alpha=\"%7.4f\" ",
             $A, $B, $vff{"bonds"}{$A}{$B}[0],
             $vff{"bonds"}{$A}{$B}[1];
      for( my $i = 2; $i < scalar( @{$vff{"bonds"}{$A}{$B}} ); $i++ )
        { printf "alpha%i=\"%7.4f\" ", $i+1, $vff{"bonds"}{$A}{$B}[$i]; }
      print "/>\n";
    }
  }
  foreach $A ( keys %{$vff{"angles"}} )
  {
    foreach $B ( keys %{$vff{"angles"}{$A}} )
    {
      foreach $C ( keys %{$vff{"angles"}{$A}{$B}} )
      {
        printf "    <Angle A=\"%s\" B=\"%s\" C=\"%s\" ",
               $A, $B, $C;
        printf "gamma=\"tet\" " if ( $vff{"angles"}{$A}{$B}{$C}[0] =~ /tet/ );
        printf "gamma=\"%7.4f\" ", $vff{"angles"}{$A}{$B}{$C}[0]
          if ( $vff{"angles"}{$A}{$B}{$C}[0] !~ /tet/ );
        printf "sigma=\"%7.4f\" beta=\"%7.4f\" ",
               $vff{"angles"}{$A}{$B}{$C}[1],
               $vff{"angles"}{$A}{$B}{$C}[2];
        for( my $i = 3; $i < scalar( @{$vff{"angles"}{$A}{$B}{$C}} ); $i++ )
          { printf "beta%i=\"%7.4f\" ", $i, $vff{"angles"}{$A}{$B}{$C}[$i]; }
        print "/>\n";
      }
    }
  }
  print "  </Functional>\n";
}



sub read_structure(\$)
{
  my $input = $_[0];
  open IN, "$input" or die " Couldn't open $input\n";

  while( ($_=<IN>) )
  {
    if( /\>(\s+|)lattice\s+vectors/ )
    {
      $_=<IN>; /(\S+)\s+(\S+)\s+(\S+)/; @{$structure{"cell"}[0]} = ($1, $2, $3);
      $_=<IN>; /(\S+)\s+(\S+)\s+(\S+)/; @{$structure{"cell"}[1]} = ($1, $2, $3);
      $_=<IN>; /(\S+)\s+(\S+)\s+(\S+)/; @{$structure{"cell"}[2]} = ($1, $2, $3);
    }
    if( />(\s+|)atomic\s+positions/ )
    {
      while( ($_=<IN>) and /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/ )
      {
        my $type=$1, $x=$2, $y=$3, $z=$4;
        my $site = -1;
        $site = 0 if( /^(\s+|)(\S+)\s+(\d+)\.(0|5)/ ); # site 1
        $site = 1 if( /(\S+)\s+(\d+)\.(2|7)/ ); # site 2
        $structure{"site $site"}{$type} = 1;
        push @{$structure{"positions"}}, (["$type", $x, $y, $z, $site]);
      }
    }
  }

  close IN;


  return if( not exists $escan{"atomic input file"} );
  return if( not -e $escan{"atomic input file"} );

  open IN, "$escan{'atomic input file'}" or return;
  <IN>;
  $_=<IN>; /(\S+)/; my $a = $1; $b = $structure{"cell"}[0][0]; 
  if( $a == 0 ) { $_=<IN>; /(\S+)/; $a = $1; $b = $structure{"cell"}[1][0]; }
  if( $a == 0 ) { $_=<IN>; /(\S+)/; $a = $1; $b = $structure{"cell"}[2][0]; }
  return if( $a == 0 );
  $structure{"scale"} = $a / $b;

  close IN;
}

sub write_structure()
{
  printf "  <Lattice scale=\"%e\">\n", $structure{"scale"} 
    if ( abs( $structure{"scale"} ) > 0.001 );
  printf "  <Lattice>\n" if ( abs( $structure{"scale"} ) < 0.001 );
  print "    <column x=\"0\" y=\"0.5\" z=\"0.5\"/>\n";
  print "    <column x=\"0.5\" y=\"0\" z=\"0.5\"/>\n";
  print "    <column x=\"0.5\" y=\"0.5\" z=\"0\"/>\n";
  print "    <site x=\"0\" y=\"0\" z=\"0\">\n";
  foreach $type ( keys %{$structure{"site 0"}} )
    { printf "      <atom type=\"%s\"/>\n", $type; }
  print "    </site>\n";
  print "    <site x=\"0.25\" y=\"0.25\" z=\"0.25\">\n";
  foreach $type ( keys %{$structure{"site 1"}} )
    { printf "      <atom type=\"%s\"/>\n", $type; }
  print "    </site>\n";
  print "  </Lattice>\n";
  print "  <Structure>\n";
  print "    <Cell>\n";
  printf "      <column x=\"%6.4f\" y=\"%6.4f\" z=\"%6.4f\" />\n",
         $structure{"cell"}[0][0], $structure{"cell"}[1][0], $structure{"cell"}[2][0];
  printf "      <column x=\"%6.4f\" y=\"%6.4f\" z=\"%6.4f\" />\n",
         $structure{"cell"}[0][1], $structure{"cell"}[1][1], $structure{"cell"}[2][1];
  printf "      <column x=\"%6.4f\" y=\"%6.4f\" z=\"%6.4f\" />\n",
         $structure{"cell"}[0][2], $structure{"cell"}[1][2], $structure{"cell"}[2][2];
  print "    </Cell>\n";
  for( my $i = 0; $i < scalar( @{$structure{"positions"}} ); $i++ )
  {
    printf "    <Atom x=\"%6.4f\" y=\"%6.4f\" z=\"%6.4f\" type=\"%s\" site=\"%i\" />\n",
           $structure{"positions"}[$i][1], $structure{"positions"}[$i][2], 
           $structure{"positions"}[$i][3], $structure{"positions"}[$i][0], 
           $structure{"positions"}[$i][4];
  }

  print "  </Structure>\n";
}
