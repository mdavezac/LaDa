#! /usr/bin/perl
#
#  Version: $Id$
#

my %data;

read_lat_in();
read_clusters_out();
read_cs_in();
write_xml();


sub read_lat_in()
{
  open IN, "<lat.in" or die "could not open lat.in\n";
  $_ = <IN>; chomp;
  push @{$data{"Axes"}}, ([ split ]);
  $_ = <IN>; chomp;
  push @{$data{"Axes"}}, ([ split ]);
  $_ = <IN>; chomp;
  push @{$data{"Axes"}}, ([ split ]);

  $_ = <IN>; chomp;
  push @{$data{"Lattice"}{"cell"}}, ([ split ]);
  $_ = <IN>; chomp;
  push @{$data{"Lattice"}{"cell"}}, ([ split ]);
  $_ = <IN>; chomp;
  push @{$data{"Lattice"}{"cell"}}, ([ split ]);

  while ( $_ = <IN>  )
  {
    chomp; s/\,//g;
    push @{$data{"Lattice"}{"sites"}}, ([ split ]);
  }

  close IN;
}

sub read_clusters_out()
{
  open IN, "<clusters.out" or die "Could not open clusters.out\n";
  open ECI, "<eci.out" or die "Could not open clusters.out\n";
  my $cluster_nb = 0;

  while(<IN>)
  {
    <IN>; # ignores first and second line
    $_=<IN>; /(\d+)/; my $nb = $1;
    # eci
    $_=<ECI>; /(\S+)/;
    $data{"clusters"}[$cluster_nb]{"E"} = $1;
    # positions
    for (my $i=0; $i < $nb; $i++)
    {
      $_=<IN>; chomp;
      push @{$data{"clusters"}[$cluster_nb]{"positions"}[$i]}, (split);
    }

    $cluster_nb++;
    <IN>;
  }
  
  close IN;
  close ECI;
}

sub read_cs_in()
{
  open IN, "<cs.in" or die "no cs.in file\n";

  <IN>; # no need for first line
  my $i = -1;
  my $new_harmonic = 0;
  while ( $_ = <IN> )
  {
    if( /(\S+)/)
    {
      if ( $new_harmonic == 0 )
        { $i ++; }
      push @{$data{"harmonics"}[$i]}, ([0,$1]);
      $new_harmonic = 1;
    }
    else 
      {  $new_harmonic = 0; }
  }
  
  close IN;

  # adds concentration stuff
  for( my $i=0; $i < scalar( @{$data{"harmonics"}} ); $i++ )
  {
    my $step = 1/( scalar( @{$data{"harmonics"}[$i]} ) - 1);
    for( my $j=0; $j < scalar( @{$data{"harmonics"}[$i]} ); $j++ )
      { $data{"harmonics"}[$i][$j][0] = $j * $step; }
  }
}

sub write_xml ()
{
  print "<?xml version=\"1.0\" standalone=no>\n";
  print "<LaDa>\n";

  printf "  <Jobs>\n";
  printf "    <Job type=\"read structures\"/>\n";
  printf "    <Job type=\"generate functionals\"/>\n";
  printf "    <--Job type=\"randomize\"/>\n";
  printf "    <Job type=\"evaluate\"/>\n";
  printf "    <--Job type=\"optimize\"/>\n";
  printf "    <--Job type=\"generate cells\"/>\n";
  printf "    <--Job type=\"print\"/>\n";
  printf "    <--Job type=\"convex hull\"/>\n";
  printf "  </Jobs>\n";

  
  printf "  <Axes>\n";
  printf "    <column x=\"%.8f\" y=\"%.8f\" z=\"%.8f\"/>\n",
         $data{"Axes"}[0][0], 
         $data{"Axes"}[0][1],
         $data{"Axes"}[0][2];
  printf "    <column x=\"%.8f\" y=\"%.8f\" z=\"%.8f\"/>\n",
         $data{"Axes"}[1][0],
         $data{"Axes"}[1][1],
         $data{"Axes"}[1][2];
  printf "    <column x=\"%.8f\" y=\"%.8f\" z=\"%.8f\"/>\n",
         $data{"Axes"}[2][0],
         $data{"Axes"}[2][1],
         $data{"Axes"}[2][2];
  printf "  </Axes>\n";

  printf "  <Lattice>\n";
  printf "    <column x=\"%.8f\" y=\"%.8f\" z=\"%.8f\"/>\n",
         $data{"Lattice"}{"cell"}[0][0],
         $data{"Lattice"}{"cell"}[0][1],
         $data{"Lattice"}{"cell"}[0][2];
  printf "    <column x=\"%.8f\" y=\"%.8f\" z=\"%.8f\"/>\n",
         $data{"Lattice"}{"cell"}[1][0],
         $data{"Lattice"}{"cell"}[1][1],
         $data{"Lattice"}{"cell"}[1][2];
  printf "    <column x=\"%.8f\" y=\"%.8f\" z=\"%.8f\"/>\n",
         $data{"Lattice"}{"cell"}[2][0],
         $data{"Lattice"}{"cell"}[2][1],
         $data{"Lattice"}{"cell"}[2][2];
  for ( my $i=0; $i < scalar( @{$data{"Lattice"}{"sites"}} ); $i++ )
  {
    printf "    <site x=\"%.8f\" y=\"%.8f\" z=\"%.8f\">\n",
           $data{"Lattice"}{"sites"}[$i][0],
           $data{"Lattice"}{"sites"}[$i][1],
           $data{"Lattice"}{"sites"}[$i][2];
    for( my $j=3; $j < scalar( @{$data{"Lattice"}{"sites"}[$i]} ); $j++ ) 
    {
      printf "      <atom type=\"%s\"/>\n",
             $data{"Lattice"}{"sites"}[$i][$j];
    }
    printf "    </site>\n",
  }
  printf "  </Lattice>\n";

  printf "  <Clusters>\n";
  for (my $cluster=0; $cluster < scalar (@{$data{"clusters"}});
       $cluster++ )
  {
    printf "    <Cluster eci=\"%.8f\">\n", $data{"clusters"}[$cluster]{"E"};
    if (exists $data{"clusters"}[$cluster]{"positions"} )
    {
      for ( my $i=0; $i < scalar( @{$data{"clusters"}[$cluster]{"positions"}} ); $i++ )
      {
        printf "      <spin x=\"%.8f\" y=\"%.8f\" z=\"%.8f\"/>\n",
               $data{"clusters"}[$cluster]{"positions"}[$i][0],
               $data{"clusters"}[$cluster]{"positions"}[$i][1],
               $data{"clusters"}[$cluster]{"positions"}[$i][2];
      }
    }
    printf "    </Cluster>\n";
  }
  printf "  </Clusters>\n";

  printf "  <CS attenuation=\"1.0\">\n";
  for( my $i=0; $i < scalar( @{$data{"harmonics"}}  ); $i++ )
  {
    printf "    <Harmonic rank=\"%i\">\n", $i;
    for( my $j=0; $j < scalar( @{$data{"harmonics"}[$i]} ); $j++ )
    { 
      printf "      <Point x=\"%.4f\" y=\"%.8f\"/>\n", 
             $data{"harmonics"}[$i][$j][0],
             $data{"harmonics"}[$i][$j][1];
    } 

    printf "    </Harmonic>\n", $i;
  }
  printf "  </CS>\n";

  printf "  <Volume min=\"2\" max=\"8\"/>\n";
  
  print "</LaDa>\n";
}

