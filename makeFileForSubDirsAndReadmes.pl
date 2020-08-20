#! /usr/local/bin/perl

# script for taking mysql dump and generating the following:
#
# 1. CSV file of subdirectory name (author_year_PMID_#) and corresponding PCL file names
# 2. Readmes for each author_year_PMID_#
# 
# 1. create a subdirectory of the author_year_PMID_#
# 2. cp the appropriate pcl files into the subdirectory
# 3. Add a readme file with the following info:
# 	 	at the top: citation, PMID, GEO ID, Full Description
#	 	one line for each pcl file - filename, short description, #conditions, tags
# 4. compress readme and PCL files together
# 5. make a zip file with all the PCL and readme files

#########################################
use strict;

#use Getopt::Mixed;

use constant TEST => 0;

#Getopt::Mixed::init('i=s input>i', 'o=s output>o');
#Getopt::Mixed::getOptions();

my %fileDataHash;
my %pclDataHash;

#our $opt_i;
#our $opt_o;

my $dataDir = "/var/www/html/data/expression/";

die ("Need input filename and output directory\n") if (! $ARGV[0] || ! $ARGV[1]); #$opt_i || ! $opt_o);

my $opt_i = $ARGV[0];
my $opt_o = $ARGV[1];

my $fileForCMS = $opt_o."microarray_expression_file_for_CMS_PROD.txt";


die ("Can't open $opt_i, file doesn't exist\n") if (! -e "$opt_i");

open(INPUT, "$opt_i");

# make hashes of data #

while (<INPUT>) {

    chomp;
    my ($pmid, $pclFilename, $shortDesc, $fullDesc, $numConds, $auth, $allAuths, $title, $journal, $year, $tags) = split("\t", $_);
    
    my @firstAuthArray = split(" ", $auth);
    my $firstInit = pop @firstAuthArray;  # remove first and middle initial
    
    my $lastName = join("_", (@firstAuthArray)); 

    my $subdir = $lastName."_".$year."_PMID_".$pmid;

    $fileDataHash{$subdir}{'pmid'} = $pmid unless ($fileDataHash{$subdir}{'pmid'});
    $fileDataHash{$subdir}{'desc'} = $fullDesc unless ($fileDataHash{$subdir}{'desc'});
    $fileDataHash{$subdir}{'cit'} = $allAuths.". ".$year.". ".$title." ".$journal unless ($fileDataHash{$subdir}{'cit'});
    push (@{$fileDataHash{$subdir}{'pcl_files'}}, $pclFilename);

    $pclDataHash{$pclFilename}{'short_desc'} = $shortDesc;
    $pclDataHash{$pclFilename}{'num_conds'} = $numConds;
    $pclDataHash{$pclFilename}{'tags'} = $tags;

}

close INPUT;

# for each key in %fileDataHash, make a subdiretory and then cp pcl files over #
# make a readme for that directory #

# check to see if the directory exists or make it

open (CMSFILE, ">$fileForCMS");

my $readmeURL = "http://downloads.yeastgenome.org/expression/microarray/";

for my $dir (sort {lc($a) cmp lc($b)} keys %fileDataHash) { ## need to sort alphabetically
    
    my $fullPath = $opt_o."/".$dir."/";
    my $readMePath = $readmeURL.$dir."/README";

    print CMSFILE "../expression/microarray/".$dir."\tDatasets (PCL file format) loaded into SPELL from ".$fileDataHash{$dir}{'cit'}."\t".$readMePath."\n";
    
    unless (-d $fullPath) {

	mkdir $fullPath or die ("Can't make subdirectory: $fullPath\n");
    }

    
    ## make README FILE ##
    my $readme = $fullPath."README";
    
    ## get GEO num ##
	
    my ($geo, $set) = split("_", $fileDataHash{$dir}{'pcl_files'}[0], 2);

    open (README, ">$readme");

    print README "Citation: ".$fileDataHash{$dir}{'cit'}."\n\n";
    print README "Full Description: ".$fileDataHash{$dir}{'desc'}."\n\n";
    print README "PMID: ".$fileDataHash{$dir}{'pmid'}."\n";

    if ($geo =~ /^GSE/) {
	print README "GEO ID: ".$geo."\n\n";
    } else {
	print README "GEO ID: N/A\n\n";
    }

    print README join("\t", ('PCL filename', 'short description', '# conditions', 'tags'))."\n";

    for my $filename (@{$fileDataHash{$dir}{'pcl_files'}}) {
	my $origFile = $dataDir.$filename;
	my $newFile = $fullPath.$filename;

	system("cp", "$origFile", "$newFile") == 0 or die("Can't copy $origFile to $newFile: $?\n");
	
	print README $filename."\t";
	print README $pclDataHash{$filename}{'short_desc'}."\t";
	print README $pclDataHash{$filename}{'num_conds'}."\t";
	print README $pclDataHash{$filename}{'tags'}."\n";
    }

    die if (TEST > 0);
}

close CMSFILE;

### 

1;  ## to keep perl happy
