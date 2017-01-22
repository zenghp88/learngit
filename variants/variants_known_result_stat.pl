#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fkey,$outdir);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"k:s"=>\$fkey,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fIn and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my %freq_stat;
my @rank=(0, 0.5, 1, 5, 10, 50, 100);

my @sample_name;
my @variant;
open (I, $fIn) or die $!;
while (<I>) {
	chomp;
	next if(/^$/);
	if($_=~/^Chr/){
		my @unit = split /\t/, $_;
		@sample_name=@unit[11 .. $#unit];
		next;
	}
	my ($c, $p,$e,$r,$a, $g, $s, $mcds, $maa, $id, $des, @sample)=split /\t/, $_;
	my $var_info=join(",", ($c, $p,$e,$r,$a, $g, $s, $mcds, $maa, $id, $des));
	push @variant, $var_info;
	for (my $i=0; $i<@sample ; $i++) {
		my ($freq)=$sample[$i]=~/_(\S+)%/;
		next if($sample[$i] eq "No");
		if ( $sample[$i] eq "-") {
			$freq=0;
		}
		my $frank=&freq_rank($freq);
		$freq_stat{"sample"}{$sample_name[$i]}{$frank}++;
		$freq_stat{"variant"}{$var_info}{$frank}++;
	}
}
close (I);


open (V,">$outdir/$fkey.variant.stat") or die $!;
print V join("\t", "Variants",0,"(0,0.5]","(0.5,1]","(1,5]","(5,10]","(10,50]","(50,100]"),"\n";
foreach my $var (@variant) {
	my @info = ($var);
	for (my $i=0; $i<@rank ; $i++) {
		my $n = exists $freq_stat{"variant"}{$var}{$rank[$i]}? $freq_stat{"variant"}{$var}{$rank[$i]}: 0;
		push @info, $n;
	}
	if ($info[-1]+$info[-2] > (scalar @sample_name)*0.1 ) {
		print "Warn High Freq Variant: $var ",$info[-1]+$info[-2],"\n";
	}
	print V join("\t", @info),"\n";
}
close (V);

open (S,">$outdir/$fkey.sample.stat") or die $!;
print S join("\t", "SampleID",0,"(0,0.5]","(0.5,1]","(1,5]","(5,10]","(10,50]","(50,100]"),"\n";
foreach my $sample (@sample_name) {
	my @info = ($sample);
	for (my $i=0; $i<@rank ; $i++) {
		my $n = exists $freq_stat{"sample"}{$sample}{$rank[$i]}? $freq_stat{"sample"}{$sample}{$rank[$i]}: 0;
		push @info, $n;
	}
	print S join("\t", @info),"\n";
}
close (S);




#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub freq_rank{
	my ($f)=@_;
	for (my $i=0; $i<@rank ; $i++) {
		if ($f<=$rank[$i]) {
			return $rank[$i];
		}
	}
}


sub AbsolutePath
{		#获取指定目录或文件的决定路径
        my ($type,$input) = @_;

        my $return;
	$/="\n";

        if ($type eq 'dir')
        {
                my $pwd = `pwd`;
                chomp $pwd;
                chdir($input);
                $return = `pwd`;
                chomp $return;
                chdir($pwd);
        }
        elsif($type eq 'file')
        {
                my $pwd = `pwd`;
                chomp $pwd;

                my $dir=dirname($input);
                my $file=basename($input);
                chdir($dir);
                $return = `pwd`;
                chomp $return;
                $return .="\/".$file;
                chdir($pwd);
        }
        return $return;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:zeng huaping<huaping.zeng\@genetalks.com> 

Usage:
  Options:
  -i  <file>   Input file, forced
  -k  <str>    key of Output file, forced

  -od <dir>    outdir of output file, default ./
  -h         Help

USAGE
	print $usage;
	exit;
}

