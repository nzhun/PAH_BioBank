
use strict;
use warnings;
use lib "/home/local/ARCS/nz2274/Pipeline/NA_script/";
use INFO2map;
sub cell2array {
	my @arr=($_[0]);
	return \@arr;
}

print "This script is used to correct the genotype and info filed for each individuals when several ALT occured in the same site in the joint calling.\n ";
print "Input format: vcf_file\n";
my $file=$ARGV[0];
my %samples=();
open IN,"zcat $file|" or die "Error: $file cannot find!\n";
open OUT, ">$file.bed";
print "Please find out in $file.bed\n";
my @header=("QUAL","FILTER","INFO","FORMAT","END","BLOCKAVG_min30p3a","SNVSB","SNVHPOL","CIGAR","RU","REFREP","IDREP","Unphased","CSQT","CSQR","AF1000G","AA","GMAF","EVS","cosmic","clinvar","phastCons","RefMinor");
#my @header=("FROM","TO","FILTER","ABHet","ABHom","AF","BaseQRankSum","ClippingRankSum","DB","DP","DS","Dels","ExcessHet","FS","GQ_MEAN","GQ_STDDEV","HRun","HapMapV3","HaplotypeScore","InbreedingCoeff","MLEAC","MLEAF","MQ","MQ0","MQRankSum","NCC","NEGATIVE_TRAIN_SITE","OND","POSITIVE_TRAIN_SITE","QD","RAW_MQ","ReadPosRankSum","SOR","VQSLOD","culprit","ANNOVAR_DATE","gnomAD_exome_ALL","gnomAD_exome_AFR","gnomAD_exome_AMR","gnomAD_exome_ASJ","gnomAD_exome_EAS","gnomAD_exome_FIN","gnomAD_exome_NFE","gnomAD_exome_OTH","gnomAD_exome_SAS","gnomAD_genome_ALL","gnomAD_genome_AFR","gnomAD_genome_AMR","gnomAD_genome_ASJ","gnomAD_genome_EAS","gnomAD_genome_FIN","gnomAD_genome_NFE","gnomAD_genome_OTH","1000g2015aug_all","1000g2015aug_eur","1000g2015aug_amr","1000g2015aug_eas","1000g2015aug_afr","1000g2015aug_sas","ExAC_ALL","ExAC_AFR","ExAC_AMR","ExAC_EAS","ExAC_FIN","ExAC_NFE","ExAC_OTH","ExAC_SAS","SIFT_score","SIFT_converted_rankscore","SIFT_pred","Polyphen2_HDIV_score","Polyphen2_HDIV_rankscore","Polyphen2_HDIV_pred","Polyphen2_HVAR_score","Polyphen2_HVAR_rankscore","Polyphen2_HVAR_pred","LRT_score","LRT_converted_rankscore","LRT_pred","MutationTaster_score","MutationTaster_converted_rankscore","MutationTaster_pred","MutationAssessor_score","MutationAssessor_score_rankscore","MutationAssessor_pred","FATHMM_score","FATHMM_converted_rankscore","FATHMM_pred","PROVEAN_score","PROVEAN_converted_rankscore","PROVEAN_pred","VEST3_score","VEST3_rankscore","MetaSVM_score","MetaSVM_rankscore","MetaSVM_pred","MetaLR_score","MetaLR_rankscore","MetaLR_pred","M-CAP_score","M-CAP_rankscore","M-CAP_pred","CADD_raw","CADD_raw_rankscore","CADD_phred","DANN_score","DANN_rankscore","fathmm-MKL_coding_score","fathmm-MKL_coding_rankscore","fathmm-MKL_coding_pred","Eigen_coding_or_noncoding","Eigen-raw","Eigen-PC-raw","GenoCanyon_score","GenoCanyon_score_rankscore","integrated_fitCons_score","integrated_fitCons_score_rankscore","integrated_confidence_value","GERP++_RS","GERP++_RS_rankscore","phyloP100way_vertebrate","phyloP100way_vertebrate_rankscore","phyloP20way_mammalian","phyloP20way_mammalian_rankscore","phastCons100way_vertebrate","phastCons100way_vertebrate_rankscore","phastCons20way_mammalian","phastCons20way_mammalian_rankscore","SiPhy_29way_logOdds","SiPhy_29way_logOdds_rankscore","Interpro_domain","GTEx_V6_gene","GTEx_V6_tissue","CADD13_RawScore","CADD13_PHRED","cosmic70","genomicSuperDups","MCAP","REVEL","avsnp147","ALLELE_END","esp6500siv2_all","esp6500siv2_aa","esp6500siv2_ea","CLINSIG","CLNDBN","CLNACC","CLNDSDB","CLNDSDBID","dbscSNV_ADA_SCORE","dbscSNV_RF_SCORE","HRC_AF","HRC_AC","HRC_AN","HRC_non1000G_AF","HRC_non1000G_AC","HRC_non1000G_AN","Kaviar_AF","Kaviar_AC","Kaviar_AN","GnomAD_Genome_covMean","GnomAD_Genome_cov10","GnomAD_Genome_cov15","GnomAD_Genome_covMedian","GnomAD_Genome_cov20","GnomAD_Genome_cov25","BRAVO_MAF","SF","AC","AN");
my $wi=0;
my $lc=0;
#push(@header,"FILTER");
while(my $line=<IN>){
	chomp($line);
	#print $line."\n";
	if ($line =~ /^##INFO/){
		next;
		#my @sets=split(/ID=|,/,$line);
		#print $sets[1]."\n";
		#push(@header,$sets[1]);
#		@header=sort(@header);
	}elsif ($line =~ /^#C/){
		my @sets=split(/\s+/,$line);
		for(my $i=9;$i<@sets;$i++){
			$samples{$i}=$sets[$i];
		}
	    next;
	
	}elsif ($line =~ /^#/){next}
	else{
		
		if($wi==0){
		#	@header= (@header,"GeneName","Transcript","Exonloc","NucleotideChange","ProteinChange");

			print OUT "#CHROM\tFROM\tTO\tREF\tALT\t".join("\t",@header)."\tproband"."\tGenotype"."\n";
		}
		$wi+=1;
		my @sets=split(/\s+/,$line);
		if(@sets<9){next;}
		my $chr=$sets[0];
		my $start=$sets[1];
		my $end=$sets[1]+length($sets[3])-1;
		#print $sets[1]."\t".$sets[3]."\n";
		my $ref=$sets[3];
		my @alts=$sets[4];
		if($sets[4] =~/,/){
			@alts=split(",",$sets[4]);
			
		}
		my %names=%{INFO2map::loadINFO($sets[7])};
		my @fv=($sets[6]);
		$names{"FILTER"}=\@fv;
	    for(my $index=0; $index<@alts;$index=$index+1){
			my $alt=$alts[$index];	
			#print "inputt: $alt\t ".scalar(@alts)."\t".$index."\n";
			my %sub_names=%{INFO2map::extractINFO(\%names,scalar(@alts),$index)};
			my $nalt=$alt;
			my $nref=$ref;
			if(length($alt)>1 && length($ref)>1){
				if(length($alt)==length($ref)){
					$nref=substr($ref,0,1);
					$nalt=substr($alt,0,1);
				}elsif(length($alt)>length($ref)){
					$nref=substr($ref,0,1);
					$nalt=substr($alt,0,length($alt)-length($ref)+1)
				}else{
					
					$nalt=substr($alt,0,1);
					$nref=substr($ref,0,length($ref)-length($alt)+1)
					
				}
			}
			#if($nalt eq "*"||($sub_names{"AN"}>1000 && $sub_names{"AC"}/$sub_names{"AN"}>0.01)){next;}
			$end=$start+length($nref)-1;
			my $outstr="$chr\t$start\t$end\t$nref\t$nalt";
			foreach my $item(@header){
				if(exists($sub_names{$item})){
					$outstr=$outstr."\t".$sub_names{$item};
				}else{
					$outstr=$outstr."\t.";
				}
			}
			for(my $id=9;$id<@sets;$id++){
				my $k=($index+1);
			    if($sets[$id] =~ /\/$k/){
#					print $outstr."\n";
					my $v=$sets[$id];
					$v =~ s/\/$k/\/1/;
			    	print OUT $outstr."\t".$samples{$id}."\t".$v."\n";
			   }
			}
		}
	}
	
}


close IN;
close OUT;
