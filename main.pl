#!/usr/bin/perl

use strict;
use warnings;
use Cwd qw(getcwd);
use Time::Local;
use Getopt::Long;

my $current_directory;

BEGIN{
    die "Your version of perl is old, please update it ASAP! \n" unless
    perl -> VERSION => 5.30;
    $current_directory = getcwd(); #We put this here to get the workign directory at initialization
}

#We get the directory where the file is located
print "The directory where binaries are located is: " . $current_directory . "\n";
use lib "$current_directory";



#We load the modules required for the analysis
use Modules::QualityControl;
use Modules::ChaosMatrixGenerator;
use Modules::Seurat_analysis;

#First, we gather all the variables involved in the program

my %variables;
my %genome_links;
my $read_settings;
my @split_variables;
my $med_variable;
my $line = 0;
my @genome_names;
my $configuration_file;
my $help;
#First we get from the options the file that will be used as ini

GetOptions("config=s" => \$configuration_file,
            "help" => \$help); 

if($configuration_file){
open(SETTINGS, $configuration_file);
while($read_settings=<SETTINGS>){
$line = $line + 1; #Update the line in which the setting is stablished, adding track for the settings value which is modified
    if($read_settings =~ /^\w./){ #If the line starts by a letter, it contains a variable to save on the hash
        #print $read_settings . "\n";
        @split_variables = split(/\s=\s?/, $read_settings); #We split the varaibles trough the equal symbol. we put the second \s in ? because the user might delete it accidentally
        $med_variable = $split_variables[1];
        chomp($med_variable);
        $variables{"$split_variables[0]"} = $med_variable;#Save the variables with the key that is the first part of the setting
    }
    else{
        next;
    }

}

close(SETTINGS); #We close the file after reading all the variables

#And we add the current working directory to the variables hash

$variables{"Binaries_directory"} = $current_directory;

#Check if all vars have been stored correctly and we pass them through the screen to see
#The execution parameters

#print "\nThe execution parameters are...\n";
#foreach my $keys (keys %variables){
#    print $keys . "=" . $variables{$keys} . "\n";
#}


#Now we load the links for all the refernece files for the genomes

open(LINKS, "links.cfg");

while($read_settings = <LINKS>){ #Although this is not the settings file, we are gonna recycle the $read_settings variable
    if($read_settings =~ /^\w./){
        #print $read_settings . "\n";
        @split_variables = split(/\s=\s/, $read_settings);
        $genome_links{$split_variables[0]} = $split_variables[1];
        push @genome_names, $split_variables[0]; #We save the name for the genomes and pass it to the cellranger command, this is ordered
    }
    else{
        next;
    }
}

close(LINKS);


#We determine the day and month in which the analysis is being done
my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
my ($sec,$min,$hour,$month_day,$month,$year,$weekday,$yearday,$isdst) = localtime();

#This way we can attribute each of the parameters of localtime to a variable
#and now we can generate the date for the folder. We have to calculate the
#year by adding 1900 to the value in $year, as this value keeps the years that has
#passed since 1900

my $actual_year = 1900 + $year;

$variables{'Actual_date'} = "$month_day\_$months[$month]\_$actual_year";
#We create the log file
my $log_name = $variables{"Working_directory"}."/SCIPA_log\-".$variables{'Actual_date'}.".txt";

#And we open the log file

unless(open LOGTXT, '>'.$log_name){
    #Die if we cant open it
    die "\nUnable to create file $log_name\n";
}

##We check that both hashes contain all the information needed:

print LOGTXT "\n----- USER SETTINGS -----\n";

foreach my $key (keys %variables){
    print LOGTXT "$key = $variables{$key}\n";
}

print LOGTXT "\n----- AVAIABLE GENOMES -----\n";

#foreach my $key (keys %genome_links){
#    print "$key web link: $genome_links{$key}\n";
#}
foreach (@genome_names){
    print LOGTXT $_ . "\n";
}
print LOGTXT "\n----- Initiating analysis -----\n";

##And now we can call the different modules on the sofware using the 
##Stored variables determined for each one of the parts

##We get the date and time in which the analysis is being performed

my $date_and_time_1 = localtime();

print "\n\nInitiating quality control analysis [$date_and_time_1] \n";
print LOGTXT "\n\nInitiating quality control analysis [$date_and_time_1] \n";

QualityControl(%variables);

my $date_and_time_2 = localtime();

print "\n\nQuality control analysis ended at [$date_and_time_2].\n\n";
print LOGTXT "\n\nQuality control analysis ended at [$date_and_time_2].\n\n";
print "\n-------------------------------------------------------------------------------\n";
print "\n\nInitiating Cellranger analysis at [$date_and_time_2]\n\n";
print LOGTXT "\n\nInitiating Cellranger analysis at [$date_and_time_2]\n\n";

$variables{"Cellranger_output"} = ChaosMatrixGeneration(\%variables, \%genome_links, \@genome_names); #We pass them as reference

my $date_and_time_3 = localtime();
print "\n\nCellranger analysis ended at [$date_and_time_3]\n\n";
print LOGTXT "\n\nCellranger analysis ended at [$date_and_time_3]\n\n";
print "Output for Cellranger located at: $variables{'Cellranger_output'}";
print LOGTXT "Output for Cellranger located at: $variables{'Cellranger_output'}";
print "\n--------------------------------------------------------------------------------------------\n";
print "\n\nInitiating Seurat based analysis of data [$date_and_time_3]\n\n";
print LOGTXT "\n\nInitiating Seurat based analysis of data [$date_and_time_3]\n\n";

SeuratAnalysis(%variables);

my $date_and_time_4 = localtime();

print "\n-------------------------------------------------------------------------------------------------------\n";
print "Analysis ended at [$date_and_time_4]. Exiting...\n";
print LOGTXT "Analysis ended at [$date_and_time_4]. Exiting...\n";

close LOGTXT;

##We assign this function to the output because we will be retourning the path to the output
##of cellranger to the key Cellranger_output on the variables hash. This way we can access
##this location on the subroutine for R analysis and perform the analysis at that position
}
elsif($help){
    print "\nSCIPA HELP\n";
    print "SCIPA pipeline is designed to automate the analysis of a set of FASTQ files from a single cell RNA seq experiment. Optimized for the use of cellranger as the main FASTQ analysis software, and use of Seurat for gene expression analysis. To execute SCIPA, please use the following command:\n\nmain.pl --config=config.ini \n";
}
else{
    print "\nNO PARAMETERS DETECTED.\n";
    print "\nSCIPA HELP\n";
    print "SCIPA pipeline is designed to automate the analysis of a set of FASTQ files from a single cell RNA seq experiment. Optimized for the use of cellranger as the main FASTQ analysis software, and use of Seurat for gene expression analysis. To execute SCIPA, please use the following command:\n\nmain.pl --config=config.ini \n";
}