package INFO2map;

use strict;
use warnings;
my $loadINFO;
sub loadINFO{
	my %names=();
	my $part=0;
	my @NULL=('NA');
	my $fullkey="";
	my @sets=split(";",$_[0]);
	
	## load info first
	for(my $j=0;$j<@sets;$j++){
		my @temp=split(/=/,$sets[$j]);
		my $key=$temp[0];
		if(@temp<2){
			$names{$key}=\@NULL;
			next;
		}
		my $value=$temp[1];
		#print $sets[$j]."\t".$key."\t".$value."\n";
		my @items=split(/,/,$value);
		if(exists($names{$key})){
			my @temp_arr=@{$names{$key}};
			#push(@temp_arr,@items);
			@items=(@temp_arr,@items)
		}
	#	print "$key\t".$value."\t".@items."\n";
		$names{$key}=\@items;
	}
	

	return \%names;
}

sub extractINFO{
	my %names=%{$_[0]};
	my $len=$_[1];
	my $index=$_[2];
	my %part_info=();
	#print "searching\t".$len."\t".$index."\n";
	foreach my $tkey(keys %names){
		my @values=@{$names{$tkey}};
		if(@values>=$len){
			$part_info{$tkey}=$values[$index]	;
	#	print "a search\t".$tkey."\t".$values[$index]."\t".$index."\t".$len."\n";
		}else{
			$part_info{$tkey}=$values[0];
		}
	}
	return \%part_info;
}




1;
