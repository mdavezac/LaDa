#! /usr/bin/perl

exit if( not -e "../../darwin/ce_opt" );

system "rm calc/krossover_VA_gen:20.agr" if ( -e "calc/krossover_VA_gen:20.agr" );
system "cd calc; ln -s ../../../darwin/ce_opt " if ( not -e "calc/ce_opt" );
for( my $i = 0; $i < 1; $i++)
{
  system "cd calc; ./ce_opt -i ../input.xml";
  system "cat calc/out.xmg >> calc/krossover_VA_gen:20.agr";
}
system "rm calc/ce_opt" if ( -e "calc/ce_opt" );

