#usr/bin/perl
##############################################################################
#
# File   :  perSampleProcess.pl
# Description: A simple script to generate the required matrix for DESeq2 input
# History:  27-May-2014 (sam) Code written
#
##############################################################################

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;

my %sampleInfo;
my %countTable;
my $sampleList="";
my $output="";
my $remove="";
my $help;
my $man;

GetOptions('help|?|h' => \$help,
			'man' => \$man, 
			'file=s' => \$sampleList,
			'remove' => \$remove,
			'out=s' => \$output) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if$man;
pod2usage("$0: No samples given.") if ($sampleList eq "" && (@ARGV == 0) && (-t STDIN));

if($output ne ""){
	#Output provided
	while(-e $output){
		print STDERR "Output file $output exists, replace it? (y or n)\t";
		my $replace = <STDIN>;
		chomp($replace);
		if($replace ne "y"){
			print STDERR "Please provide a new output file name (or Ctrl+C to terminate):\t";
			my $output = <STDIN>;
			chomp($output);
			unless(-e $output){
				last;
			}
		}
		else{
			print STDERR "Will replace the file\n";
			last;
		}
	}
}

if($sampleList eq ""){
	#obtain information from argument
	print STDERR "Getting information from argument\n";
	foreach(@ARGV){
		#print STDERR $_."\n";
		chomp($_);
		my @info = split(/:/, $_);
		if(@info < 2){
			print STDERR "Incorrect input format, argument input format should be <Sample Name>:<Count File>\n";
			print STDERR "$_\n";
			print STDERR "Terminate\n";
			exit -1;
		}
		#The input is correct
		if(exists($sampleInfo{$info[0]})){
			print STDERR "Duplicated sample: $info[0] !";
			print STDERR "Terminate\n";
			exit -1;
		}
		#Format is correct, and there is no duplicated samples, now check if the file input exists
		if(-e $info[1]){
			$sampleInfo{$info[0]} = $info[1];
			print STDERR "$info[0]\t$info[1]\tadded\n";
		}
		else{
			print STDERR "$info[1] not found, please check your input\n";
			print STDERR "Terminate\n";
			exit -1;
		}
	}
}
else{
	#obtain information from the file
	if(-e $sampleList){
		open(FILE, $sampleList);
		my $lineNum = 0;
		for my $line (<FILE>){
			$lineNum = $lineNum+1;
			#The first should be the sample name, whereas the second should be the file lcation
			chomp($line);
			my @info = split(/\t/, $line);
			if(@info < 2){
				print STDERR "Incorrect input format, input should contain at least two column and should be tab delimited\n";
				print STDERR "Error at line:$lineNum\n";
				print STDERR "Terminate\n";
				exit -1;
			}
			if(exists($sampleInfo{$info[0]})){
				print STDERR "Duplicated sample: $info[0] !";
				print STDERR "Terminate\n";
				exit -1;
			}
			#Format is correct, and there is no duplicated samples, now check if the file input exists
			if(-e $info[1]){
				$sampleInfo{$info[0]} = $info[1];
				print STDERR "$info[0]\t$info[1]\tadded\n";
			}
			else{
				print STDERR "$info[1] not found, please check your input\n";
				print STDERR "Terminate\n";
				exit -1;
			}
		}
		close(FILE);
		#Now check if we have more than 2 samples
	}
	else{
		print STDERR "$sampleList not found, please check if you have the correct input\n";
	}
}
if(keys( %sampleInfo ) < 2){
	print STDERR "You only have one sample, doesn't require to build matrix\n";
	print STDERR "Terminate\n";
	exit -1;
}

#Now start building the matrix using hash of array
my $sampleID = -1;
foreach my $key (keys %sampleInfo){
	$sampleID = $sampleID+1;
	print STDERR "Reading $sampleInfo{$key} for $key\n";
	open(FILE, $sampleInfo{$key});
	for my $line (<FILE>){
		#Now read the information, should be tab delim, with two column
		chomp($line);
		my @info = split(/\t/, $line);
		if(@info != 2){
			print STDERR "Number of column isn't 2, please check if you have the correct htseq-count output\n";
			print STDERR "$line\n";
			print STDERR "Terminate\n";
			exit -1;
		}
		#elsif($remove && ($info[0] eq "no_feature" || $info[0] eq "ambiguous" || $info[0] eq "too_low_aQual" || $info[0] eq "not_aligned" || $info[0] eq "alignment_not_unique")){
			#Ignore field
		#}
		else{
			#Format is correct
			if(exists($countTable{$info[0]})){
				#This is not the first sample
				@{$countTable{$info[0]}}[$sampleID] = $info[1];
			}
			elsif($sampleID == 0){
				#This is the first sample, so it is normal for it to be novel
				my @count = (-1) x keys(%sampleInfo);
				$count[$sampleID] = $info[1];
				$countTable{$info[0]} = \@count;
			}
			else{
				#The is not the first, something is wrong.
				print STDERR "Cannot find $info[0] in other files, please use the same GTF file for htseq-count\n";
				print STDERR "Terminate\n";
				exit -1;
			}
		}
	}
	close(FILE);
}

#output information
if($output ne ""){
	open FILE, ">", $output or die "$0: open $output: $!";
	
	foreach my $key (keys %sampleInfo){
		print FILE ",$key";
	}
	print FILE "\n";
	foreach my $key (keys %countTable){
		print FILE "$key";
		foreach my $next (@{$countTable{$key}}){
			print FILE ",$next";
		}
		print FILE "\n";
	}
	close FILE, $output or die "$0: close $output: $!";
}
else{
	foreach my $key (keys %sampleInfo){
		print STDOUT ",$key";
	}
	print STDOUT "\n";
	foreach my $key (keys %countTable){
		print STDOUT "$key";
		foreach my $next (@{$countTable{$key}}){
			print STDOUT ",$next";
		}
		print STDOUT "\n";
	}
}
__END__


=head1 NAME

generate Matrix - Form matrix from the given files

=head1 SYNOPSIS

perl perlSampleProcess.pl [options] [Sample:Sample_count_file ...]

Options:

  -file     Provide file input of samples instead of through command line argument
  -remove   Remove the last 5 fields from the count matrix
  -out      output file name (default: STDOUT)
  -help		brief help message 
  -man		full documentation
  
=head1 DESCRIPTION

B<This programme> will read in the count files for each sample and generate a count matrix in csv format

=head1 OPTIONS

=over 8 

=item B<-file>

Provide the file containing the sample information instead of through the command line argument.
File should be of the following format (tab delimited)

Sample_Name		Count_File

Please also note that all the count file should be generated using the same GTF file, or else the gene name might not be compatable.

=item B<-out>

Output file name. (default: STDOUT)

=item B<-remove>

If specified, will remove the last 5 lines from the count matrix:

no_feature, ambiguous, too_low_aQual, not_aligned, alignment_not_unique

Will affect the size factor normalization used in DESeq2

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exists.

=back

=head1 ACKNOWLEDGEMENT

Thank you

=head1 AUTHOR

Choi Shing Wan - L<https://github.com/choishingwan>

=head1 LICENSE

GPLv3

=cut