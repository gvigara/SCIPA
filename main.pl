#!/usr/bin/perl

use strict;
use warnings;
use Cwd qw(getcwd);
use Time::Local;

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

open(SETTINGS, "config.ini");
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

##We check that both hashes contain all the information needed:

print "\n----- USER SETTINGS -----\n";

foreach my $key (keys %variables){
    print "$key = $variables{$key}\n";
}

print "\n----- AVAIABLE GENOMES -----\n";

#foreach my $key (keys %genome_links){
#    print "$key web link: $genome_links{$key}\n";
#}
foreach (@genome_names){
    print $_ . "\n";
}
print "\n----- Initiating analysis -----\n";

##We get the date and time in which the analysis is being performed

my $date_and_time_1 = localtime();

##And now we can call the different modules on the sofware using the 
##Stored variables determined for each one of the parts

print "\n\nInitiating quality control analysis [$date_and_time_1] \n";

QualityControl(%variables);

my $date_and_time_2 = localtime();

print "\n\nQuality control analysis ended at [$date_and_time_2].\n\n";
print "\n-------------------------------------------------------------------------------\n";
print "\n\nInitiating Cellranger analysis at [$date_and_time_2]\n\n";

$variables{"Cellranger_output"} = ChaosMatrixGeneration(\%variables, \%genome_links, \@genome_names); #We pass them as reference

my $date_and_time_3 = localtime();
print "\n\nCellranger analysis ended at [$date_and_time_3]\n\n";
print "Output for Cellranger located at: $variables{'Cellranger_output'}";
print "\n--------------------------------------------------------------------------------------------\n";
print "\n\nInitiating Seurat based analysis of data [$date_and_time_3]\n\n";

SeuratAnalysis(%variables);

##We assign this function to the output because we will be retourning the path to the output
##of cellranger to the key Cellranger_output on the variables hash. This way we can access
##this location on the subroutine for R analysis and perform the analysis at that position