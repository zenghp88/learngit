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
my $fvariant = "/database/zenghp/bin/annovar/humandb/hg19_49Gene_cosmic.txt";
my ($fanno,$sample,$outdir);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fanno,
				"v:s"=>\$fvariant,
				"s:s"=>\$sample,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fanno and $sample);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my %detect;

#Chr	Start	End	Ref	Alt	Func.refGene	Gene.refGene	GeneDetail.refGene	ExonicFunc.refGene	AAChange.refGene	avsnp144	49Gene_cosmic	Otherinfo
#chr1	27087403	27087403	T	C	exonic	ARID1A	.	synonymous SNV	ARID1A:NM_006015:exon5:c.T1977C:p.P659P,ARID1A:NM_139135:exon5:c.T1977C:p.P659P	.	.	247,2_0.80%
#chr1	27087415	27087415	A	G	exonic	ARID1A	.	synonymous SNV	ARID1A:NM_006015:exon5:c.A1989G:p.T663T,ARID1A:NM_139135:exon5:c.A1989G:p.T663T	.	.	277,1_0.36%
#chr1	27087439	27087439	A	G	exonic	ARID1A	.	synonymous SNV	ARID1A:NM_006015:exon5:c.A2013G:p.G671G,ARID1A:NM_139135:exon5:c.A2013G:p.G671G	.	.	327,1_0.30%
#

open (O,">$outdir/$sample.final.txt") or die $!;
open (I,"$fanno") or die $!;
my $head=<I>;
print O $head;
my @hunit = split /\t/, $head;
my $target_index = (scalar @hunit) -1-1;
while (<I>) {
	chomp;
	my @unit = split /\t/, $_;
	next if($unit[$target_index] eq "."); #
	$detect{$unit[$target_index]}=$unit[$target_index+1];
	print O $_,"\n";
}
close (I) ;
close (O) ;

open (O,">$outdir/$sample.variants_known.result") or die $!;
#1       27087413        27087413        A       T       ID=COSM4143737;GENE=ARID1A;STRAND=+;CDS=c.1987A>T;AA=p.T663S;CNT=1
#1       27087417        27087417        C       G       ID=COSM79135;GENE=ARID1A;STRAND=+;CDS=c.1991C>G;AA=p.S664*;CNT=2
#1       27087418        27087420        AGG     A       ID=COSM3358384;GENE=ARID1A;STRAND=+;CDS=c.1993_1994delGG;AA=p.G665fs*10;CNT=1
print O "Chr\tStart\tEnd\tRef\tAlt\tGene\tStrand\tMut-CDS\tMut-AA\tMutID\tMut-description\t$sample\n";
open (V,"$fvariant") or die $!;
while (<V>) {
	chomp;
	my ($chr, $start, $end, $ref, $alt, $info)=split /\t/, $_;
	my @unit=split /;/, $info;
	my %value;
	for (my $i=0;$i<@unit ;$i++) {
		my ($key,$value)=split /=/, $unit[$i];
		$value{$key}=$value;
		if ($key eq "CNT" || $key eq "SNP" || $key eq "TARGET") {
			push @{$value{"DES"}}, $unit[$i];
		}

	}
	
	my @des=exists $value{"DES"}? @{$value{"DES"}}: (".");
	my $detect = exists $detect{$info}? $detect{$info}: "-";
	my $gene = exists $value{"GENE"}? $value{"GENE"}: ".";
	my $strand = exists $value{"STRAND"}? $value{"STRAND"}: ".";
	my $cds = exists  $value{"CDS"}?  $value{"CDS"}: ".";
	my $aa = exists $value{"AA"}? $value{"AA"}: ".";
	my $id = exists $value{"ID"}? $value{"ID"}: ".";
	print O join("\t", "chr".$chr, $start,$end,$ref,$alt, $gene, $strand, $cds, $aa, $id, join(";",@des), $detect),"\n";
}

close (V) ;
close (O);

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

	note: target_annovar must be the last two column.

Usage:
  Options:
  -i  <file>   Input annovar result file(hg19_multianno.txt), forced
  -v  <file>   Input annovar variant known file(hg19_49Gene_cosmic.txt), [/database/zenghp/bin/annovar/humandb/hg19_49Gene_cosmic.txt]
  -s  <str>     sample name , forced
  -od <dir>    Dir of output file, default ./
  -h         Help

USAGE
	print $usage;
	exit;
}
