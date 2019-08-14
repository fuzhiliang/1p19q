#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename qw(basename dirname);
use Statistics::Descriptive;
use Math::Round;
use Data::Dumper;
my ($inhouse_snp,$path_ffpe,$path_blood,$sample,$out,$up,$down,$bam);
GetOptions(
                                "help|?" =>\&USAGE,
                                "snp:s"=>\$inhouse_snp,
                                "path_ffpe:s"=>\$path_ffpe,
				"path_bc:s"=>\$path_blood,
				"sampleid:s"=>\$sample,
				"up:s"=>\$up,
				"down:s"=>\$down,
				"path_bam:s"=>\$bam,
                                "outdir:s"=>\$out,
) or &USAGE;
&USAGE unless ($inhouse_snp and $path_ffpe and $path_blood and $out and $sample);
`mkdir $out` unless(-d $out);
#===========================================================================default value
$down||=0.4;
$up||=0.6;
#===========================================================================loh analysis
`mkdir -p $out/work_sh` unless(-d "$out/work_sh");
open (LOH,">$out/work_sh/$sample.loh.sh") or die $!;
print LOH "/lustre/rdi/user/licq/software/Anaconda3/bin/perl /lustre/rdi/user/licq/research/LOH/Script_V2/loh_analysis.oncotarget_v2.bak.pl  -snp $inhouse_snp -sampleid $sample -path_ffpe $path_ffpe -path_bc $path_blood -down $down -up $up -outdir $out\n";

`sh $out/work_sh/$sample.loh.sh`;

#===========================================================================call ctcnv
if (-d $bam ){
	my $bamfile = "$bam/$sample.sorted.rmdup_new.realign.bam";
	if (-f $bamfile){
		open (CNV,">$out/work_sh/$sample.ctcnv.sh") or die $!;
		print CNV "/lustre/rdi/user/yujn/software/Anaconda2/bin/python /lustre/rdi/user/yujn/tools/ctCNV.py autocall --inBam $bamfile --refCov /lustre/rdi/user/licq/project/LOH/second/120ratio/BC_negative/plot/second/ctCNV/Analysis/GBM_ref_COV.txt --refGCS /lustre/rdi/user/licq/project/LOH/second/120ratio/BC_negative/plot/second/ctCNV/Analysis/GBM_ref_GCS.txt --output $sample --output-dir $out/$sample --bed /lustre/rdi/user/licq/project/LOH/NO_BC_LOH/ctCNV/GBM.merge.1p19q.bed --fasta /lustre/rdi/user/hantc/tools/hg19/ucsc.hg19.fasta\n";
		close CNV;
	}
	`sh $out/work_sh/$sample.ctcnv.sh`;
}else{

	print "Please check if the sample bam file exists $bam directory \n";die;
}
#==========================================================================summary cnv and loh result

if (-f "$out/$sample/$sample\_SCNA_results.txt"){
	my %cnv;
	open (IN,"less $out/$sample/$sample\_SCNA_results.txt|awk '{print \$1\"\t\"\$NF}'|grep -E -w '1p|1q|19p|19q'|") or die $!;
	while (<IN>){
		chomp;
		next if (/^$/ || /^#/);
		my @lines = split(/\s+/,$_);
		$cnv{$lines[0]} = $lines[1];
	}	
	close IN;
	open (IN1,"$out/$sample/$sample.summary.result") or die $!;
	open (OUT,">$out/$sample/$sample.loh.cnv.result") or die $!;
	while (<IN1>){
		chomp;
		next if (/^$/ || /^#/);
		my @lines = split(/\t+/,$_);
		my $Chr1p = join("\t",@lines[0..2]);
		my $Chr1q = join("\t",@lines[3..4]);
		my $Chr19p = join("\t",@lines[6..7]);
		my $Chr19q = join("\t",@lines[8..9]);
		my $IDH12 =  join("\t",@lines[10..12]);
		if ($lines[0] eq "sample"){
			print OUT "$Chr1p\t1p_copy number\t$Chr1q\t1q_copy number\t$lines[5]\t$Chr19p\t19p_copy number\t$Chr19q\t19q_copy number\t$IDH12\n";
		}else{
			print OUT "$Chr1p\t$cnv{'1p'}\t$Chr1q\t$cnv{'1q'}\t$lines[5]\t$Chr19p\t$cnv{'19p'}\t$Chr19q\t$cnv{'19q'}\t$IDH12\n";
		}
	}
	close IN1;
	close OUT;
}

#===========================================================================Sub usage

sub USAGE {
        my $usage=<<"USAGE";
Description:
Usage: perl $Bin/loh_analysis.oncotarget.pl -snp snpfile -path_ffpe ffpe_varscan.snp.gz -path_bc bc_varscan.snp.gz -sampleid sampleid -outdir outdir
Options:
-snp            <file>  The SNP sites involved in the bed file [rs 1p/19q chr pos]
-path_ffpe	<file>  The SNP file corresponding to the tumor sample[ffpe.varscan.snv.vcf.gz]
-path_bc	<file>  The SNP file corresponding to the normal sample[bc.varscan.snv.vcf.gz]
-sampleid       <str>   sampleid
-up		<float> A measure of the upper limit of heterozygosity[optional,default:0.6]
-down           <float> A measure of the lower limit of heterozygosity[optioanl,defalut:0.4]
-path_bam	<float> The sample corresponds to the path of the bam file[optional]
-outdir         <file>  output directory
-h         Help
USAGE
        print $usage;
        exit;
}
