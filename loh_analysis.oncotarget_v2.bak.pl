#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename qw(basename dirname);
use Statistics::Descriptive;
use Math::Round;
use Data::Dumper;
my ($inhouse_snp,$path_ffpe,$path_blood,$sample,$out,$up,$down);
GetOptions(
                                "help|?" =>\&USAGE,
                                "snp:s"=>\$inhouse_snp,
                                "path_ffpe:s"=>\$path_ffpe,
				"path_bc:s"=>\$path_blood,
				"sampleid:s"=>\$sample,
				"up:s"=>\$up,
				"down:s"=>\$down,
                                "outdir:s"=>\$out,
) or &USAGE;
&USAGE unless ($inhouse_snp and $path_ffpe and $path_blood and $out and $sample);
`mkdir $out` unless(-d $out);
#===================================================snp on the 1p and 19q 
my %dist=('chr1'=>500000,'chr19'=>10000);
$down||=0.4;
$up||=0.6;

open(SNP,$inhouse_snp) or die "$!";
my %snp;
while(<SNP>){
	chomp;
	my $rs=(split(/\t/,$_))[0];
	my $type=(split(/\t/,$_))[1];
	$snp{$rs}=$type;
}
close(SNP);
my (@panelsnp,@sample);
`mkdir $out/$sample` unless (-d "$out/$sample");
my $outdir = "$out/$sample";
my $file_ffpe = $path_ffpe;
my $file_blood = $path_blood;
chomp($file_ffpe);chomp($file_blood);
#====================================================IDH1 and IDH2 mutation && annotation
my %depth_ffpe;
open(IN1,"zcat $file_ffpe|") or die "$!";
open(ANNOVAR,">$outdir/$sample.avinput");
while(<IN1>){
	next if($_=~/#/);
	my @data=split(/\t/,$_);
	my @info=split(/\:/,$data[9]);
	$depth_ffpe{$data[0]}{$data[1]}{$data[3]}{$data[4]}="$info[4]:$info[5]:$info[3]";
	if($data[0] eq 'chr2'){
		if($data[1]>=209113111 and $data[1]<=209113113){
			print ANNOVAR "$data[0]\t$data[1]\t$data[1]\t$data[3]\t$data[4]\t$info[5]/$info[3]\n";
		}
	}
	elsif($data[0] eq 'chr15'){
		if($data[1]>=90630405 and $data[1]<=90630407){
			print ANNOVAR  "$data[0]\t$data[1]\t$data[1]\t$data[3]\t$data[4]\t$info[5]/$info[3]\n";
		}
	}
}
close(IN1);
close(ANNOVAR);
unless (-z "$outdir/$sample.avinput"){
`perl /work/app/annovar_for_zhanggh_scripts/table_annovar.pl $outdir/$sample.avinput /work/app/annovar_for_zhanggh_scripts/humandb/ -buildver hg19 -out $outdir/$sample.avinput -remove -protocol refGene -operation g -nastring NA`;
}
#====================================================Mutation sites Shared by BC and FFPE samples were extracted
my %depth_blood;
open(IN2,"zcat $file_blood|") or die "$!";
while(<IN2>){
	next if($_=~/#/);
	my @data=split(/\t/,$_);
	my @info=split(/\:/,$data[9]);
	$depth_blood{$data[0]}{$data[1]}{$data[3]}{$data[4]}="$info[4]:$info[5]:$info[3]";
}
close(IN2);
open(OUT2,">$outdir/$sample.snv.txt");
foreach my $chr(keys %depth_ffpe){
	foreach my $pos(keys %{$depth_ffpe{$chr}}){
		foreach my $ref(keys %{$depth_ffpe{$chr}{$pos}}){
			foreach my $alt(keys %{$depth_ffpe{$chr}{$pos}{$ref}}){
				if(exists $depth_blood{$chr}{$pos}{$ref}{$alt}){
					print OUT2 "$chr\t$pos\t$pos\t$ref\t$alt\t$depth_ffpe{$chr}{$pos}{$ref}{$alt}\t$depth_blood{$chr}{$pos}{$ref}{$alt}\n";
				}
			}
		}
	}
}
close(OUT2);
`perl /work/app/annovar_for_zhanggh_scripts/table_annovar.pl $outdir/$sample.snv.txt /work/app/annovar_for_zhanggh_scripts/humandb/  -buildver hg19 -out $outdir/$sample.snv.txt -remove -protocol refGene,snp138 -operation g,f -nastring NA --otherinfo`;
#====================================================The loci on 1q and 19p were extracted
my $threshold = &threshold(%snp);

#====================================================The loci on 1p and 19q were extracted
open(IN3,"sort -k1,1 -k2,2n $outdir/$sample.snv.txt.hg19_multianno.txt|") or die "$!";
open(OUT,">$outdir/$sample.snv.1p19q.ratio.txt") or die "$!";
open(OUT4,">$outdir/$sample.all.snv.ratio.txt") or die "$!";
my ($num,$before_chr,$next_chr,$before_pos,$next_pos,%hash1,%hash2);
my $part = 1;
my $pos=0;
while(<IN3>){
	chomp;
	my @data=split(/\t/,$_);
	next if($data[0] ne 'chr1' and $data[0] ne 'chr19');
	next if($data[0] eq 'chr1' and $data[1] >125000000);
	next if($data[0] eq 'chr19' and $data[1] <26500000);
	my $rs;
	next if($data[10] eq 'NA');
	if($data[10]=~/,/){
		$rs=(split(/\,/,$data[10]))[0];
	}
	else{
		$rs=$data[10];
	}
	$rs=~s/Name=//;
	my $ratio;
	if (exists $snp{$rs}){
		my @ffpe=split(/:/,$data[11]);
		my @blood=split(/:/,$data[12]);
		my $allele2_ffpe=$ffpe[1]/$ffpe[2];
		my $allele1_ffpe=$ffpe[0]/$ffpe[2];
		my $allele2_blood=$blood[1]/$blood[2];
		my $allele1_blood=$blood[0]/$blood[2];
		$allele1_blood=0.00001 if($allele1_blood==0);
		$allele2_blood=0.00001 if($allele2_blood==0);
		$allele1_ffpe=0.00001 if($allele1_ffpe==0);
		$allele2_ffpe=0.00001 if($allele2_ffpe==0);
		next if($allele2_blood<0.2 or $allele2_blood>0.8);
		$ratio=($allele2_blood/$allele1_blood)/(($allele2_blood/$allele1_blood)+($allele2_ffpe/$allele1_ffpe));
		my $status='normal';
		$status='aberrant' if($ratio<$down or $ratio >$up);
		my $tmp="$sample\t$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$rs\t$ratio\t$status";
		print OUT4 "$tmp\n";
		$hash1{$part}{'normal'}=0 unless (exists $hash1{$part}{'normal'});
		$hash1{$part}{'aberrant'}=0 unless (exists $hash1{$part}{'aberrant'});
		$hash2{$part}{'normal'}="-" unless (exists $hash2{$part}{'normal'});
		$hash2{$part}{'aberrant'}="-" unless (exists $hash2{$part}{'aberrant'});
		$hash2{$part}{$status}=$tmp;
		$hash1{$part}{$status}=1 if ($hash1{$part}{$status}==0);
		if(($data[1] - $pos)<$dist{$data[0]}){
			$hash1{$part}{$status}++;
		}else{
			$part++;
		}
		$pos=$data[1];

	}
}

close(IN3);
close(OUT4);

foreach my $block (sort {$a cmp $b} keys %hash1){
	my $status;
	if ($hash1{$block}{'normal'} > $hash1{$block}{'aberrant'}){
		$status = 'normal';
		print OUT "$hash2{$block}{$status}\n";
	}else{
		$status = 'aberrant';
		print OUT "$hash2{$block}{$status}\n";
	}
}
close(OUT);
#=====================================================Plot the ratio value for loci on the chr1 and chr19 
`cat $outdir/$sample.snv.1p19q.ratio.txt $outdir/$sample.snv.1q19p.ratio.txt |sort -k2,2 -k3,3n>$outdir/$sample.snv.chr1_chr19.ratio.txt`;
open (RS,">$outdir/$sample.snv.R") or die $!;
print RS "library(ggplot2)\n";
print RS "data <- read.table(\"$outdir\/$sample.snv.chr1_chr19.ratio.txt\",header = F)\n";
print RS "data1 <- data[,c(2,3,8)]\n";
print RS "fish <- data.frame(chr = c(\"chr1\",\"chr1\", \"chr19\",\"chr19\"),  pos = c(0, 7.2, 43.4,51.4))\n";
print RS "pq <- data.frame(chr = c(\"chr1\", \"chr19\"),  pos = c(125,26.5))\n";
print RS "names(data1) <- c(\"chr\",\"pos\",\"ratio\")\n";
print RS "df <- data.frame(x = data1\$pos,y=data1\$ratio,chr=data1\$chr)\n";
print RS "p <- ggplot(df,aes(x=x/1000000,y=y)) + geom_point(size = 0.9) + facet_wrap(~chr,scales = \"free\",nrow = 2) + xlab(\"pos\") + ylab(\"ratio\") + labs(title=\"$sample\") + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = \"red\") + scale_y_continuous(limits = c(0,1))  + geom_vline(aes(xintercept = pos),color = \"red\", fish) + geom_vline(aes(xintercept = pos) , color = 'red', size = 0.2, linetype=2,pq)+ geom_hline(yintercept = 0.4 , color = 'red', size = 0.2, linetype=\"solid\") + geom_hline(yintercept = 0.6 , color = 'red', size = 0.2, linetype=\"solid\") \n";
print RS "ggsave(p, file = \"$outdir\/$sample.snv.ratio.png\", width = 12, height = 8, dpi = 100)\n";
close RS;
`/lustre/rde/user/rde_admin/bin/Rscript $outdir/$sample.snv.R`;
#=====================================================Summary result
my @pq_ratio = ("$outdir/$sample.snv.1p19q.ratio.txt","$outdir/$sample.snv.1q19p.ratio.txt");
for (my $i=0;$i<@pq_ratio;$i++){
open(RESULT,"sort -k2,2 -k3,3n $pq_ratio[$i]|") or die "$!";
my $pq_type = (split(/.ratio/,basename($pq_ratio[$i])))[0];
my ($num1,$aberrant_num);
my ($before_sample,$before_chr1,$next_sample,$next_chr1);
my (%loh,$per);
$aberrant_num=0;
while(<RESULT>){
	chomp;
	my($sample,$chr,$start,$end,$ref,$alt,$rs,$ratio,$status)=split(/\t/,$_);
	$num1++;
	if($num1==1){
		$before_sample=$sample;
		$before_chr1=$chr;
		$aberrant_num++ if($status eq 'aberrant');
	}
	else{
		$next_sample=$sample;
		$next_chr1=$chr;
		if($next_sample ne $before_sample or $next_chr1 ne $before_chr1){
			$num1=$num1 - 1;
			$per = $aberrant_num/$num1;
			if($aberrant_num>1){
				if($aberrant_num == $num1){
					$loh{$before_sample}{$before_chr1}="LOH ($aberrant_num/$num1)\t$per";
				}
				else{
					$loh{$before_sample}{$before_chr1}="partial LOH ($aberrant_num/$num1)\t$per";
				}
			}
			else{
				$loh{$before_sample}{$before_chr1}="NO LOH ($aberrant_num/$num1)\t$per";
			}
			$num1=1;
			$aberrant_num=0;
		}
		$aberrant_num++ if($status eq 'aberrant');
		$before_sample=$sample;
		$before_chr1=$chr;
	}
}
close(RESULT);
$per = $aberrant_num/$num1;
if($aberrant_num>1){
	if($aberrant_num == $num1){
		$loh{$before_sample}{$before_chr1}="LOH ($aberrant_num/$num1)\t$per";
	}
	else{
		$loh{$before_sample}{$before_chr1}="partial LOH ($aberrant_num/$num1)\t$per";
	}
}
else{
	$loh{$before_sample}{$before_chr1}="NO LOH ($aberrant_num/$num1)\t$per";
}
my @chr=('chr1','chr19');
foreach my $sample(keys %loh){
	open (SING,">$out/$sample/$pq_type.snv.result.txt") or die $!;
	if ($pq_type =~/1p19q/){
		print SING "sample\t1p\tratio_1p\t19q\tratio_19q\tIDH1\tIDH2\n$sample";
	}else{
		print SING "sample\t1q\tratio_1q\t19p\tratio_19p\tIDH1\tIDH2\n$sample";
	}
	for(my $i=0;$i<@chr;$i++){
		if(exists $loh{$sample}{$chr[$i]}){
			print SING "\t$loh{$sample}{$chr[$i]}";
		}
		else{
			print SING "\tNA\tNA";
		}
	}
	my ($idh1,$idh2);
	if (-f "$out/$sample/$sample.avinput.hg19_multianno.txt"){
		$idh1=`grep "IDH1" $out/$sample/$sample.avinput.hg19_multianno.txt`;
		$idh2=`grep "IDH2" $out/$sample/$sample.avinput.hg19_multianno.txt`;
		chomp($idh1);chomp($idh2);
	}
	if($idh1){
		my $mutation=(split(/\t/,$idh1))[9];
		my $mutation2=(split(/,/,$mutation))[0];
		print SING "\t$mutation2";
	}
	else{
		print SING "\tNA";
	}
	if($idh2){
		my $mutation=(split(/\t/,$idh2))[9];
		my $mutation2=(split(/,/,$mutation))[0];
		print SING "\t$mutation2";
	}
	else{
		print SING "\tNA";
	}
	print SING "\n";
	close SING;
}
}
open (SUM,"cat $outdir/$sample.snv.1p19q.snv.result.txt $outdir/$sample.snv.1q19p.snv.result.txt|grep -v 'ratio'|") or die $!;
open (SUMM,">$outdir/$sample.summary.result") or die $!;
print SUMM "sample\t1p\t1p_ratio\t1q\t1q_ratio\t1q/1p\t19p\t19p_ratio\t19q\t19q_ratio\t19p/19q\tIDH1\tIDH2\n";
while (<SUM>){
	chomp;
	next if (/^$/ || /^#/);
	my @lines = split(/\t/,$_);
	print SUMM "$lines[0]\t$lines[1]\t$lines[2]\t";
	my $q19p = <SUM>;
	chomp($q19p);
	my @lines_q19p = split(/\t/,$q19p);
	my ($p_1_q_1,$p_19_q_19,$p_1_q_1_2,$p_19_q_19_2);
	if ($lines[2] eq 0 || $lines[2] eq 'NA' || $lines_q19p[2] eq 'NA'){
		$p_1_q_1_2  = "-";
	}else{
		$p_1_q_1_2  = $lines_q19p[2]/$lines[2];
	}
	if ($lines[4] eq 0 || $lines[4] eq 'NA' || $lines_q19p[4] eq 'NA'){
		$p_19_q_19_2 = "-";
	}else{
		$p_19_q_19_2 = $lines_q19p[4]/$lines[4];
	}
	print SUMM "$lines_q19p[1]\t$lines_q19p[2]\t$p_1_q_1_2\t$lines_q19p[3]\t$lines_q19p[4]\t$lines[3]\t$lines[4]\t$p_19_q_19_2\t$lines[5]\t$lines[6]\n";
}
close SUM;
close SUMM;
#=====================================================Sub function
sub threshold{
	my (%snp) = @_;
	my %ratio;
	open(INFILE,"sort -k1,1 -k2,2n $outdir/$sample.snv.txt.hg19_multianno.txt|") or die "$!";
	open(OUTT,">$outdir/$sample.snv.1q19p.ratio.txt");
	while(<INFILE>){
        	chomp;
		my @data=split(/\t/,$_);
		next if($data[0] ne 'chr1' and $data[0] ne 'chr19');
		next if($data[0] eq 'chr1' and $data[1] <125000000);
		next if($data[0] eq 'chr19' and $data[1] >26500000);
		my $rs;
		next if($data[10] eq 'NA');
		if($data[10]=~/,/){
                	$rs=(split(/\,/,$data[10]))[0];
       		} else{
               		$rs=$data[10];
       		}
       		 $rs=~s/Name=//;
		my $ratio;
	        #if (exists $snp{$rs}){
           		my @ffpe=split(/:/,$data[11]);
               		my @blood=split(/:/,$data[12]);
                        my $allele2_ffpe=$ffpe[1]/$ffpe[2];
                        my $allele1_ffpe=$ffpe[0]/$ffpe[2];
                        my $allele2_blood=$blood[1]/$blood[2];
                        my $allele1_blood=$blood[0]/$blood[2];
                        $allele1_blood=0.00001 if($allele1_blood==0);
                        $allele2_blood=0.00001 if($allele2_blood==0);
                        $allele1_ffpe=0.00001 if($allele1_ffpe==0);
                        $allele2_ffpe=0.00001 if($allele2_ffpe==0);
                        if($allele2_blood<0.2 or $allele2_blood>0.8){
                              next;
                        }else{
                 #             if ($snp{$rs} eq "1q" or $snp{$rs} eq "19p"){
                                    my $status='normal';
                                    $ratio=($allele2_blood/$allele1_blood)/(($allele2_blood/$allele1_blood)+($allele2_ffpe/$allele1_ffpe));
				    $status='aberrant' if($ratio<$down or $ratio >$up);
				    print OUTT "$sample\t$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$rs\t$ratio\t$status\n";
				    if ($data[0] eq 'chr1'){
                                        push @{$ratio{'1q'}},$ratio;
				    }else{
					push @{$ratio{'19p'}},$ratio;
				    }
                  #            }
                        }
               # }
         }
         close INFILE;
	 close OUTT;
=cut
         my %flag;
         foreach my $type (keys %ratio){
		my @value = sort @{$ratio{$type}};
		#my $value = &mid(@{$ratio{$type}});	
		$flag{$type} = "$value[-1]_$value[0]";
		#$flag{$type} = "$value"; 
		print "$type\t$flag{$type}\n";
	 }
	 return \%flag;
          
=cut
}

sub mid {
	my @list = sort @_;
	my $count = @list;
	if ($count == 0){
		return 0;die;
	}
	my $up = $list[round(($count-1)*0.9)];
	my $down = $list[round(($count-1)*0.1)];
	my $th = "$up\_$down";
	return $th;
}
#====================================================Sub usage

sub USAGE {#
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
-outdir         <file>  output directory
-h         Help
USAGE
        print $usage;
        exit;
}
