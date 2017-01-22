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
my ($fvariant, $fvariant_cmp, $fcrowd, $fblack, $fkey,$outdir);
GetOptions(
				"help|?" =>\&USAGE,
				"iv:s"=>\$fvariant,
				"icmp:s"=>\$fvariant_cmp,
				"ic:s"=>\$fcrowd,
				"ib:s"=>\$fblack,
				"k:s"=>\$fkey,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fvariant and $fcrowd and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
my $min_1000eas = 0.001;

my %detect;
my %detect_old;
my %info;
my %freq_1000eas;

open (I,"$fvariant") or die $!;
my $head = <I>;
chomp $head;
my (@uhead)=split /\t/, $head;
my $info_head=join("\t",@uhead[0..($#uhead-1)]);
while (<I>) {
	chomp;
	next if(/^$/);
	my (@unit)=split /\t/, $_;
	my $info = join("\t",@unit[0..($#unit-1)]);
	my $det = $unit[-1];
	my $freq_1000eas = $unit[12];
	my ($id)=$_=~/ID=(\S+);GENE=/;

	$detect{$id}=$det.";".$freq_1000eas;
	$info{$id}=$info;
}
close (I) ;

my %black;
if (defined $fblack) {
	open (I,"$fblack") or die $!;
	while (<I>) {
		chomp;
		next if(/^$/);
		my ($id, $info)=split /\s+/, $_;
		$black{$id}=$info;
	}
	
	close (I) ;
}

if (defined $fvariant_cmp) {
	open (I,"$fvariant_cmp") or die $!;
	while (<I>) {
		chomp;
		next if(/^$/ || /Chr/);
		my (@unit)=split /\t/, $_;
		my $info = join("\t",@unit[0..($#unit-1)]);
		my $det = $unit[-1];
		my $freq_1000eas = $unit[12];
		my ($id)=$_=~/ID=(\S+);GENE=/;
		$detect_old{$id}=$det.";".$freq_1000eas;
		$info{$id}=$info;
	}
	close (I) ;
}

my %crowd_inter;
open (I,"$fcrowd") or die $!;
while (<I>) {
	chomp;
	next if(/^$/);
	my ($id)=(split /,/, $_)[9];
	my ($is_high, $inter)=(split /\t/, $_)[8,13];
	if ($is_high eq "Yes") {
		$crowd_inter{$id}=$inter;
	}
}
close (I) ;

my %filter;
#==================================================================
# filter: 两次是否一致
#==================================================================
my %compare;
if (defined $fvariant_cmp) {
	### compare
	my %stat;
	$stat{"NumNow"} = scalar keys %detect;
	$stat{"NumOld"} = scalar keys %detect_old;
	foreach my $var (keys %detect) {
		if (!exists $detect_old{$var}) {
			$filter{$var}="NoAgree";
			push @{$compare{"NoAgree"}{$var}}, ($detect{$var}, "-");
		}else{
			$stat{"Agree"}++;
			push @{$compare{"Agree"}{$var}}, ($detect{$var}, $detect_old{$var});
			delete $detect_old{$var};
		}
	}
	foreach my $var (keys %detect_old) {
		$filter{$var}="NoAgree";
		push @{$compare{"NoAgree"}{$var}}, ("-",$detect_old{$var});
	}

	### output stat
	open (S,">$outdir/$fkey.compare.stat") or die $!;
	print S "NumNow\tNumOld\tAgree\n";
	print S join("\t",$stat{"NumNow"}, $stat{"NumOld"}, exists $stat{"Agree"}? $stat{"Agree"}: 0),"\n";
	close (S) ;
}else{
	foreach my $var (keys %detect) {  ## 后期基于$compare{"Agree"} hash进行分析
		push @{$compare{"Agree"}{$var}},$detect{$var};
	}
}

#==================================================================
# filter
#==================================================================
foreach my $var (keys %{$compare{"Agree"}}) {
	my @filter;
	for (my $i=0; $i<@{$compare{"Agree"}{$var}} ; $i++) {
		my $det = $compare{"Agree"}{$var}->[$i];
		my ($dr, $dv, $freq, $f1000eas)=$det=~/(\d+),(\d+)_([1234567890.]+)%;([1234567890.]+)/;
		print $var,"\t",$det,"\t", join("\t",($dr, $dv, $freq, $f1000eas)),"\n";
		
		## filter black list
		if (exists $black{$var}) {
			$filter[$i] .= "BlackFilter:$black{$var}";
			next;
		}

		## filter high freq crowd
		if (exists $crowd_inter{$var}) {
			my ($min, $max)=split /-/, $crowd_inter{$var};
			if ($freq >= $min && $freq <= $max) {
				$filter[$i] .= "CrowdFilter";
				next;
			}else{
				$filter[$i] .= "HighCrowd_";
			}
		}

		## filter poly variants-- f>40% and freq in 1000g2014oct_eas >=0.001
		if ( $freq>40 && $f1000eas ne "." && $f1000eas>=$min_1000eas) {
			$filter[$i] .= "PolyFilter";
			next;
		}

		## filter_loose
		if ($dv<2 || $freq<0.5) {
			$filter[$i] .= "FilterLoose";
			next;
		}

		## filter strict
		if ($dv<2 || $freq < 0.8 ) {
			$filter[$i] .= "FilterStrict";
			next;
		}
		print "YES",$var,"\t",$det,"\t", join("\t",($dr, $dv, $freq, $f1000eas)),"\n";
		
		$filter[$i] .= "PASS";
	}

	## pass when exist one pass
	$filter{$var}=join(",",@filter);
}

## output
open (FR,">$outdir/$fkey.remain_repeat.txt") or die $!;
open (FP,">$outdir/$fkey.remain_report.txt") or die $!;
open (F,">$outdir/$fkey.filter.info") or die $!;
my %stat_filter;

print F $info_head,"\tDetect";
print FR $info_head,"\tDetect";
print FP $info_head,"\tDetect";
if (defined $fvariant_cmp) {
	print F "\tDetectOld";
	print FR "\tDetectOld";
	print FP "\tDetectOld";
}
print F "\tFilter\n";
print FR "\n";
print FP "\n";

foreach my $var (sort {$info{$a} cmp $info{$b}} keys %info) {
	$stat_filter{"Total"}++;
	my @det_info = exists $compare{"Agree"}{$var}? @{$compare{"Agree"}{$var}}: @{$compare{"NoAgree"}{$var}};
	for (my $i=0; $i<@det_info ; $i++) {
		my ($det,$f1000eas)=split /;/, $det_info[$i];
		$det_info[$i] = $det;
	}
	if (exists $compare{"NoAgree"}{$var}) {
		$stat_filter{$filter{$var}}++;
		print F $info{$var},"\t",join("\t",@det_info),"\t",$filter{$var},"\n";
		next;
	}
	print F $info{$var},"\t",join("\t",@det_info),"\t",$filter{$var},"\n";
	$filter{$var}=~s/HighCrowd_//g;
	$stat_filter{$filter{$var}}++;
	
	if ($filter{$var}!~/Filter/ ) { ## both must be pass
		$stat_filter{"Report"}++;
		$stat_filter{"Repeat"}++;
		print FR $info{$var},"\t",join("\t",@det_info),"\n";
		print FP $info{$var},"\t",join("\t",@det_info),"\n";
	}elsif($filter{$var}!~/FilterLoose/ && $filter{$var}!~/CrowdFilter/ && $filter{$var}!~/PolyFilter/){
		$stat_filter{"Repeat"}++;
		print FR $info{$var},"\t",join("\t",@det_info),"\n";
	}
}

close (FR) ;
close (FP) ;
close (F) ;

open (S,">$outdir/$fkey.filter.stat") or die $!;
print S "remain_report\t",exists $stat_filter{"Report"}? $stat_filter{"Report"}: 0,"\n";
print S "remain_repeat\t",exists $stat_filter{"Repeat"}? $stat_filter{"Repeat"}: 0,"\n";
print S "Total\t", exists $stat_filter{"Total"}? $stat_filter{"Total"}: 0,"\n";
foreach my $type (sort {$a cmp $b} keys %stat_filter) {
	next if($type eq "Report" || $type eq "Repeat" || $type eq "Total");
	print S $type,"\t",$stat_filter{$type},"\n";
}
close (S);
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
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

	v2: v>2 f>0.8%
	v3: Variants both detected by two times must be both v>2 f>0.8%
	v4: io change to icmp
	v5: filter poly variants-- f~50% and freq in 1000g2014oct_eas >=0.001
	v6: add filter variants in black list

Usage:
  Options:
  -iv    <file>   Input variant file, forced
  -ic    <file>   Input crowd file, forced
  -icmp  <file>   Input compared variant file, optional
  -ib    <file>   Input ID black list, optional

  -k  <str>    Key of output file, forced
  -od <dir>    Dir of output file, default ./
  -h         Help

USAGE
	print $usage;
	exit;
}
