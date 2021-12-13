#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $info = "/home/users/xiaodong/Documents/Project/Seal/sample_list_v2/20210119/r2/greyseal.world.rmLow3andWSpecies.rmPlate1.v2.bamlist.info.csv";
my $bamlist = shift;
my $ibs = shift;


open (my $fh0, "<$info") or die $!;
my $header0 = <$fh0>;
my %hash;
while (<$fh0>) {
    chomp;
    my @row = split (",", $_);
    (my $file, my $pop) = @row[(1,4)];
    $hash{$file} = $pop;
}
close $fh0;

open (my $fh1, "<$bamlist") or die $!;
my %hash2;
my $idx = 0;
while(<$fh1>) {
    chomp;
    $hash2{$idx} = $hash{$_};
    ++$idx;
}
close $fh1;

my %hash3;
my @target = ("GREY","KOD","END","HOK","CAL","NEW","BIC","MET","NFL","GRE","ICE","SVA","LIS","ORK","ROL","SOT","WAS","WNL","LIE","ANH","MAA","ROD");
foreach my $pop (@target) {
    my @target_keys = grep { $hash2{$_} eq $pop } keys %hash2;
    $hash3{$pop} = [@target_keys];
}


open (my $fh2, "gunzip -c $ibs |" ) or die "gunzip $ibs: $!";
my $header = <$fh2>;
my %loci;
print join("\t",@target), "\n";

while (<$fh2>) {
    chomp;
    my @line = split /\t/, $_;
    my $no_miss = scalar(grep /-1$/, @line);
    next if $no_miss >=30;
    (my $chr,my $pos,my $major, my $mino) = @line[0..3];
    my @output ;
    
    foreach my $i (@target) {
	my @alleles = @line[map {$_+4} @{$hash3{$i}}];
	my $major_alleles = scalar(grep /^1$/, @alleles);
	my $minor_alleles = scalar(grep /0$/, @alleles);
	my $missing = scalar(grep /^-1$/, @alleles);
	my $out = join ",", ($major_alleles,$minor_alleles);
	push @output, $out;
    }
    print join("\t", @output), "\n";
}
