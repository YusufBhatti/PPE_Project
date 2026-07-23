#! /usr/bin/env perl
#-----------------------------------------------------------------------------

print " EXECUTION TIME [S]  MFLOP/S      FMA [%]  CPI     PEAK [%]  REGION\n";
print "-------------------------------------------------------------------\n"; 
#      1                    2            3        4       5        6

open HPM, "<echam5.hpm";

$label = "";

while (<HPM>) {
    if ( /Instrumented section/ ) {
	(@line) = split;
	$label = $line[5];
    }
    if ( /Execution time/ ) {
	(@line) = split;
	$time = $line[6];
    }
    if ( /Instructions per run cycle/ ) {
	(@line) = split;
	$cpi = 1.0 / $line[5];
    }
    if ( /Algebraic flop rate/ ) {
	(@line) = split;
	$flops = $line[7];
    }
    if ( /FMA percentage/ ) {
	(@line) = split;
	$fma = $line[3];
    }
    if ( /peak performance/ && $label ne "" ) {
	(@line) = split;
	$peak = $line[5];
	printf "%19.8f  %10.3f  %7.3f  %6.2f   %6.3f    %s\n", 
	  $time, $flops, $fma, $cpi, $peak, $label;	
	if ( $label eq "all" ) {
	    print "-------------------------------------------------------------------\n"; 
	}
	$label = "";
    }
}

print "-------------------------------------------------------------------\n"; 
close HPM;
