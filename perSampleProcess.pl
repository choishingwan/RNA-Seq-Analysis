#/usr/bin/perl
##############################################################################
#
# File   :  perSampleProcess.pl
# History:  27-May-2014 (sam) Code written
#
##############################################################################

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);


sub escape { $_[0] =~ s/([^a-zA-Z0-9_])/\\$1/g; return $_[0]; }

#This function is used to see if a certain programme is located in the path
sub location{
	my $tool= $_[0];
	for my $path (split /:/, $ENV{PATH}){
		if( -f "$path/$tool" && -x _){
			$path = escape($path);
			#print "$tool found in $path\n";
			#print "Will set $path/$tool as the default $tool location\n";
		return "$path/$tool";
		}
	}
	return "";
};

#This function is used to check the memory input format
sub checkFormat{
	my $memory = shift;
	my $typeLong = shift;
	my $typeShort = shift;
	my $longLoc = index($memory, $typeLong);
	my $shortLoc = index($memory, $typeShort);
	if($longLoc != -1){
		if(length($memory)-2 ne $longLoc){
			return -1;
		}
		my $num = substr $memory, 0, -2;
		looks_like_number($num) ? return $num : return -1;
	}
	elsif($shortLoc != -1){
		if(length($memory)-1 ne $shortLoc){
			return -1;
		}
		my $num = substr $memory, 0, -1;
		looks_like_number($num) ? return $num : return -1;
	}
	else{
		return -1;
	}
}

#This function will check if both reads are of the same format and will return whether if gz is required
sub checkSample{
	my $read1 = shift;
	my $read2 = shift;
	#Check it contain file type info
	if(index($read1, ".") == -1 || index($read2, ".") == -1){
		print STDERR "The read files doesn't contain file type info, cannot proceed\n";
		return -1;
	}
	my @oneDecom = split(/\./, $read1);
	my @twoDecom = split(/\./, $read2);
	if($oneDecom[@oneDecom-1] eq $twoDecom[@twoDecom-1]){
		if(lc($oneDecom[@oneDecom-1]) ne "fastq" &&  lc($oneDecom[@oneDecom-1]) ne "fq" && lc($oneDecom[@oneDecom-1]) ne "gz"){
			print STDERR "The read files are not of the correct format (allowed = \x27*.fastq\x27, \x27*.fq\x27, \x27*.gz\x27)\n";
			return -1;
		}
		lc($oneDecom[@oneDecom-1]) eq "gz" ? return 1 : return 0;
	}
	else{
		print STDERR "The read files are of different format, cannot proceed: $read1\t$read2\n";
		return -1;
	}
	
}

sub checkSuffix{
	my $input = shift;
	my $longForm = shift;
	my $shortForm = shift;
	my @inList = split(/\./,$input);
	if(lc($inList[@inList-1]) ne $longForm && lc($inList[@inList-1]) ne $shortForm){
		print STDERR "Incorrect format: $inList[@inList-1] for $input\n";
		return 0;
	}
	return 1;
}

#The only thing that we are going to do here is to write the file
sub produceFile{
	my $sampleName = shift;
	my $read1 = shift;
	my $read2 = shift;
	my $sam = shift;
	my $star = shift;
	my $htseq = shift;
	my $lib = shift;
	my $readLength = shift;
	my $thread = shift;
	my $samThread = shift;
	my $memory = shift;
	my $reference=  shift;
	my $gtf = shift;
	my $index = shift;
	my $gz =shift;
	
	my $outName = $sampleName.".sh";

	
	while(-e $outName){
		print STDERR "Output file $outName exists, replace it? (y or n)\t";
		my $replace = <>;
		chomp($replace);
		if($replace ne "y"){
			print STDERR "Please provide a new output file name:\t";
			my $outName = <>;
			chomp($outName);
			unless(-e $outName){
				last;
			}
		}
		else{
			print STDERR "Will replace the file\n";
			last;
		}
	}
	
	open(OUTFILE,">", $outName);
	my $information =
	"
	Sample Information             $sampleName
	First Read                     $read1
	Second Read                    $read2
	Stranded                       $lib
	Read Length                    $readLength
	Index Directory                $index
	GTF File                       $gtf";
	if(defined($reference)){
		$information = $information."\n	Reference Fasta                $reference";
	}
	if(defined($gz)){
		$information = $information."
	gunzip input                   yes";
	}
	else{
		$information = $information."
	gunzip input                   no";
	}
	$information = $information."
	STAR                           $star
	htseq-count                    $htseq
	samtools                       $sam";
	if(defined($samThread)){
		$information = $information."
	samtools Thread                $samThread
	samtools max mem               $memory";
	}
	$information = $information."
	thread                         $thread\n";
	print OUTFILE ": <<\x27END\x27\n";
	print OUTFILE $information;
	print OUTFILE "END\n";
	unless(-e "$index/Genome"){
		#Need to build index;
		print OUTFILE "echo \"\$(date) Build index\"\n";
		my $junctionHang = $readLength-1;
		print OUTFILE "$star --runMode genomeGenerate --genomeDir $index --genomeFastaFiles $reference --runThreadN $thread --sjdbGTFfile $gtf --sjdbOverhang $junctionHang\n";
		print OUTFILE "#===================================================================================\n\n";
	}
	my $alignment = "";
	if($gz){
		$alignment = "awk -F \"_\" \x27{print \"$star --genomeDir $index --readFilesIn $read1 $read2 --runThreadN $thread genomeLoad=LoadAndRemove --readFilesCommand zcat --outFileNamePrefix $sampleName. -- outSAMunmapped Within\"}\x27";
	}
	else{
		$alignment = "awk -F \"_\" \x27{print \"$star --genomeDir $index --readFilesIn $read1 $read2 --runThreadN $thread genomeLoad=LoadAndRemove --outFileNamePrefix $sampleName. -- outSAMunmapped Within\"}\x27";
	}
	print OUTFILE "echo \"\$(date) Performing alignment on $sampleName, will produce $sampleName.sam\"\n";
	print OUTFILE $alignment."\n";
	print OUTFILE "#===================================================================================\n\n";
	
	#Perform the sorting, depending on the samThread option, will consider whether if multithread is used.
	
	my $samSort="";
	if(defined($samThread)){
		#Use multi-thread
		$samSort = "$sam view -bSh $sampleName.sam | $sam sort -n -m $memory -@ $samThread - $sampleName.sorted";
	}
	else{
		$samSort = "$sam view -bSh $sampleName.sam | $sam sort -n - $sampleName.sorted";
	}
	print OUTFILE "echo \"\$(date) Performing Sorting, will produce $sampleName.sorted.bam\"\n";
	print OUTFILE $samSort."\n";
	print OUTFILE "#===================================================================================\n\n";
	
	#Get the counts
	my $counting = "$sam view $sampleName.sorted.bam | $htseq -s $lib - $gtf > $sampleName.counts";
	print OUTFILE "echo \"\$(date) Performing counting, will produce $sampleName.counts\"\n";
	print OUTFILE $counting."\n";
	
	
	close(OUTFILE);
	
	print $information;
	print "Completed\n";
}

#Performing the documentation
my $man = 0;
my $help = 0;
my $sampleName = "";
my $readLength = 101;
my $lib = "yes";
my $thread = 1;
my $star ="";
my $htseq = "";
my $read1 = "";
my $read2 = "";
my $sam="";
my $memory="20gb";
my $reference = "";
my $gtf = "";
my $index ="";  #Genome
#BEGIN PARSE COMMANDS
GetOptions('help|?|h' => \$help,
			'man' => \$man, 
			'sampleName|name=s' => \$sampleName, 
			'length=i' => \$readLength, 
			'strand=s'=> \$lib, 
			'thread|t=i' => \$thread, 
			'star=s' => \$star,
			'htseq=s' => \$htseq,
			'samtools=s' => \$sam,
			'first=s' => \$read1,
			'second=s' => \$read2,
			'memory=s' => \$memory,
			'ref=s' => \$reference,
			'gtf=s' => \$gtf,
			'index=s' => \$index) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if$man;
#pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my $error = '';

if($star eq ""){
	$star = location("STAR");
	
	if($star eq ""){
		print STDERR "Please provide the required STAR programme location\n";
		$error =1;
	}
}
if($htseq eq ""){
	$htseq = location("htseq-count");
	if($htseq eq ""){
		print STDERR "Please provide the required htseq-count programme location\n";
		$error =1;
	}
}

if($sam eq ""){
	$sam = location("samtools");
	if($sam eq ""){
		print STDERR "Please provide the required htseq-count programme location\n";
		$error =1;
	}
}

print STDERR "\n";
if($sampleName eq "" && ($read1 eq "" || $read2 eq "")){
	print STDERR "Please provide either the sample name or directly provide the location of the reads\n";
	print STDERR "sampleName or reads not provided\n";
	$error = 1;
}
if($readLength < 1){
	print STDERR "Read length cannot be smaller than 0\n";
	$error = 1;
}
if($thread < 1){
	print STDERR "Must use at least 1 thead\n";
	$error =1;
}
if(lc($lib) ne "yes" && lc($lib) ne "no" && lc($lib) ne "reverse"){
	print STDERR "Strand option specify whether if the sequencing is strand specific.\n";
	print STDERR "Can only be \x27yes\x27, \x27no\x27 or \x27reverse\x27\n";
	$error =1;
}

if($index eq ""){
	print STDERR "Please provide the index directory\n";
	$error = 1;
}
if($gtf eq ""){
	print STDERR "Please provide the GTF file\n";
	$error = 1;
}
else{
	#check file suffix
	if(!checkSuffix($gtf, "gtf", "gff")){
		$error = 1;
	}
}
unless( -e "$index/Genome"){
	if($reference eq ""){
		print STDERR "Didn't find the Genome file within $index\n";
		print STDERR "Please provide the reference file to build the reference\n";
		$error = 1;
	}
	else{
		#check file suffix
		if(!checkSuffix($reference, "fa", "fasta")){
			$error = 1;
		}
	}
}



#Check the version of samtools
my @testVersion = `$sam sort 2>&1`;
my $samThread = "";
foreach (@testVersion){
	if(index($_, "@") != -1){
		$samThread = $thread;
		last;
	}
}
if(!$samThread){
	print STDERR "The current samtools doesn't support threading, will disable threading for samtools\n";
}

if($samThread ne ""){
	#Can run threading in samtools, will require to calculate the required memory
	#check whether if the memory input is correct
	#Levels are gb, mb, kb
	$memory = lc($memory);
	my $format;
	if(($format = checkFormat($memory, "gb", "g")) != -1){
		#print $format."\n";
		#correct format, now divide the memory to each thread and round up to nearest integer
		$memory = int $format/$samThread;
		if($memory eq 0){
			#Too small, need to go down a level
			$memory = int $format*1024/$samThread;
			if($memory eq 0){
				#Go down to the last level
				$memory = int $format*1024*1024/$samThread;
				if($memory eq 0){
					print STDERR "Insufficient memory for samtools threading, will not perform threading when performing samtools sort\n";
					$samThread = "";
				}
				else{
					$memory = $memory."K";
				}
			}
			else{
				$memory = $memory."M";
			}
		}
		else{
			$memory = $memory."G";
		}
	}
	elsif(($format = checkFormat($memory, "mb", "m")) != -1){
		#correct format, now divide the memory to each thread and round up to nearest integer
		$memory = int $format/$samThread;
		if($memory eq 0){
			#Too small, need to go down a level
			$memory = int $format*1024/$samThread;
			if($memory eq 0){
				#Not enough memory
				print STDERR "Insufficient memory for samtools threading, will not perform threading when performing samtools sort\n";
				$samThread = "";
			}
			else{
				$memory = $memory."K";
			}
		}
		else{
			$memory = $memory."M";
		}
	}
	elsif(($format = checkFormat($memory, "kb", "k")) != -1){
		#correct format, now divide the memory to each thread and round up to nearest integer
		$memory = int $format/$samThread;
		if($memory eq 0){
			#Not enough memory
			print STDERR "Insufficient memory for samtools threading, will not perform threading when performing samtools sort\n";
			$samThread = "";
		}
		else{
			$memory = $memory."K";
		}
	}
	else{
		print STDERR "Incorrect memory input format: $memory\n";
		$error = 1;
	}
}



if(($read1 eq "" || $read2 eq "") && $sampleName ne ""){
	#Try to make the file
	my @list = qx{ls};
	foreach(@list){
		if(index($_, $sampleName) != -1 && (index(lc($_), ".fq") != -1 || index(lc($_), ".fastq") != -1)){
			#Contain the sample name, check if it contain the fq or fastq in the file name
			if(index($_, 1) != -1){
				#Likely be the first sample
				$read1 = $_;
				$read1 =~ s/^\s+|\s+$//g;
				print STDERR "Will use $read1 as the first read, correct? (y or n) \n";
				my $ok=<>;
				chomp($ok);
				if($ok ne "y"){
					print STDERR "Please provide the file directly\n";
					$error = 1;
					last;
				}
				#$read1 =~ s/[\$#@~!&*()\[\];.,:?^ `\\\/]+//g;
			}
			elsif(index($_, 2) != -2){
				$read2 = $_;
				$read2 =~ s/^\s+|\s+$//g;
				print STDERR "Will use $read2 as the second read, correct? (y or n) \n";
				my $ok=<>;
				chomp($ok);
				if($ok ne "y"){
					print STDERR "Please provide the file directly\n";
					$error = 1;
					last;
				}
				#$read2 =~ s/[\$#@~!&*()\[\];.,:?^ `\\\/]+//g;
			}
		}
	}
	if($read1 eq "" || $read2 eq ""){
		print STDERR "Cannot identify possible fastq files in current directory using $sampleName\n";
		print STDERR "You might want to directly provides the sample directory\n";
		$error = 1;
	}
}


my $gz;
if($read1 ne "" && $read2 ne ""){
	my $gz = checkSample($read1, $read2);
	if($gz eq -1){
		$error = 1;
	}
}


print STDERR "\n";


pod2usage(1) if (($error) && (-t STDIN));


#END PARSE COMMANDS

if($sampleName eq ""){
	my @token = split(/\./, $read1);
	if(index($token[0], "_") != -1){
		my @temp = split(/_/, $token[0]);
		$sampleName = $temp[0];
	}
	else{
		$sampleName = $token[0];
	}
}

#print "$sampleName\t$readLength\t$memory\t$samThread\n";

produceFile($sampleName, $read1, $read2, $sam, $star, $htseq, $lib, $readLength, $thread, $samThread, $memory, $reference, $gtf, $index,  $gz);


__END__


=head1 NAME

perSampleProcess - Perform basic sample processing

=head1 SYNOPSIS

perl perlSampleProcess.pl [options] 

Options:

  -sampleName
  -first
  -second
  -ref
  -gtf
  -index
  -star
  -htseq
  -lib
  -length
  -thread
  -memory
  -help		brief help message
  -man		full documentation

=head1 DESCRIPTION

B<This programme> will read the given input and produce the script for performing the basic RNA Sequencing analysis. 

=head1 OPTIONS

=over 8 

=item B<-sampleName>

The sample name, if reads were not given, will try to search the local directory for <sampleName>_1.fq, <sampleName>_2.fq or <sampleName>_1.fq.gz, <sampleName>_2.fq.gz
If reads were given, then the output will use sampleName as the prefix

=item B<-first>

The first fastq file. 
The programme will check whether if the file contain fastq or fq and if it is gz.
If sample name was not provided, the sample name will be obtained from the fastq. 
For example, for Sample_1.fq, the sample name will be defined as Sample

=item B<-second>

The second fastq file.
The programme will check whether if the file contain fastq or fq and if it is gz.

=item B<-ref>

The reference fasta file.
Only require if index isn't built.

=item B<-gtf>

The GTF file containing the splicing information. Recommend to use the Ensemble gtf file.

=item B<-index>

The reference index directory.
Will not build the reference if it contains the Genome file (assume built).

=item B<-star>

The STAR aligner programme. 
Will use the STAR programme in <PATH> as default (if available)

=item B<-htseq>

The htseq-count programme. 
Will use the htseq-count in <PATH> as default (if available)

=item B<-strand>

whether the data is from a strand-specific assay. 
Specify 'yes', 'no', or 'reverse' (default: yes).
'reverse' means 'yes' with reversed strand


=item B<-length>

The read length, used for STAR alignment (default: 101).

=item B<-thread>

The number of thread used for alignment (default: 1).

=item B<-memory>

The allowed memory for the analysis(default: 20gb).
Used for calculating the available threads for samtools

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exists.

=back

=head1 CAVEATS

Problems

=head1 ACKNOWLEDGEMENT

Thank you

=head1 AUTHOR

Choi Shing Wan - L<https://github.com/choishingwan>

=head1 LICENSE

GPLv3

=cut