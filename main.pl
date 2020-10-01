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
my $read_settings;
my @split_variables;
my $med_variable;
my $line = 0;

open(SETTINGS, "config.ini");
while($read_settings=<SETTINGS>){
$line = $line + 1; #Update the line in which the setting is stablished, adding track for the settings value which is modified
# GENERAL SETTINGS

    if($read_settings =~ /^(Working_directory)./){ #Read the working directory on the config file
        @split_variables = split(/\s=\s/, $read_settings); #Split the line by the =
        if($split_variables[1] =~/\/\w/){ #If it starts by a / (indicating linux or mac directory)
            $med_variable = $split_variables[1]; 
            chomp($med_variable);
            $variables{"Working_directory"} = $med_variable; #Save it to the directory hash key
        }
        else{ #In the case that there is not specified directory or wrong, the program dies
            die "\nNo $split_variables[0] specified under section [General], line $line. Please check config. ini \n";

        }
    }
    elsif($read_settings =~ /^(Cellranger_cores)./){ #Read the core threads designated for analysis
        @split_variables = split(/\s=\s/, $read_settings);
        if($split_variables[1] =~ /[0-9]+/ && $split_variables[1] > 0){ #If it matches a number
            $med_variable = $split_variables[1];
            chomp($med_variable);
            $variables{"Core_threads"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [General], line $line. Please check config. ini \n";

        }
    }
    elsif($read_settings =~ /^(Cellranger_mem)./){ #Read the RAM memory amount designated for analysis
        @split_variables = split(/\s=\s/, $read_settings);
        if($split_variables[1] =~/[0-9]+/ && $split_variables[1] > 0){
            $med_variable = $split_variables[1];
            chomp($med_variable);
            $variables{"RAM_memory"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [General], line $line. Please check config. ini \n";
        }
    }

#QUALITY CHECK SETTINGS

    elsif($read_settings =~ /^(Output_QC)./){ #Read the output directory for FASTQC analysis
        @split_variables = split(/\s=\s/, $read_settings);
        if($split_variables[1] =~ /\.\/\w/){ #If the setting starts with ./ indicating an extension of the WD
            $split_variables[1] =~ s/\.//g; #Remove the 
            $med_variable = $split_variables[1]; #Intermediate varaible to modify the string
            chomp($med_variable); #Remove the /n from the string
            $variables{"Output_QC"} = $variables{"Working_directory"} . $med_variable;
        }
        elsif($split_variables[1] =~ /^\/.+/){ #If it starts with the / directory meaning another output
            $med_variable = $split_variables[1];
            chomp($med_variable);
            $variables{"Output_QC"} = $med_variable;
        }
        else{ #There is no wrong thing on not putting a directory on this section, it will use the working directory as output
            $variables{"Output_QC"} = $variables{"Working_directory"};
        }
    }
    elsif($read_settings =~ /^(Multi_QC)./){ #Read if the user wants to perform a MultiQC analysis
        @split_variables = split(/\s=\s/, $read_settings);
        if($split_variables[1] =~ /y[es]?/s or $split_variables[1] =~ /no?/s){
            $med_variable = $split_variables[1];
            chomp($med_variable);
            $variables{"MultiQC_analysis"} = $med_variable;
        }
        else{
            die "\nIncorrect $split_variables[0] specified under section [QC analysis], line $line. Please check config. ini \n";
        }
    }

#CELLRANGER SETTINGS

    elsif($read_settings =~ /^(Cellranger_path)./){ #We get the path where cellranger is installed
        @split_variables = split(/\s=\s/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable eq ""){ #With this loop we check if the parameter is specified. If not, the program stops. 
            die "\nNo $split_variables[0] specified under section [Cellranger], line $line. Please check config. ini \n";
        }
        else{
            $variables{"Cellranger_path"} = $med_variable;
        }
    }
    elsif($read_settings =~ /^(Cellranger_id)./){ #Information about the id for the analysis
        @split_variables = split(/\s=\s/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable eq ""){
            die "\nNo $split_variables[0] specified under section [Cellranger], line $line. Please check config. ini \n";
        }
        else{
            $variables{"Cellranger_id"} = $med_variable;
        }
    }
    elsif($read_settings =~ /^(Cellranger_organism)./){ #We gather information about the organism. 
        @split_variables = split(/\s=\s/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable eq ""){
            die "\nNo $split_variables[0] specified under section [Cellranger], line $line. Please check config. ini \n";
        }
        elsif($med_variable =~ /[H-h]uman/ || $med_variable =~ /[M-m]ouse/){ #Only for human or mouse, if not, program stops. 
            $variables{"Cellranger_org"} = $med_variable;
        }      
        else{
            die "\nOrganism in line $line not recognized. Please check config.ini \n";
        } 
    }elsif($read_settings =~ /^(Cellranger_reference)./){ #Reference genome. If empty, will be downloaded. 
        @split_variables = split(/\s=\s/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable eq ""){
            print "\nNo $split_variables[0] specified under section [Cellranger], line $line.\n";
            $variables{"Cellranger_ref"} = $med_variable;
        }
        else{
            $variables{"Cellranger_ref"} = $med_variable;
        }       
    }
    elsif($read_settings =~ /^(Expected_cells)./){ #The amount of expected cells for the analysis
        @split_variables = split(/\s=\s/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable eq ""){ #In this case, this parameter is not mandatory. Thus, if it has no input, the program will continue. 
            print "\nNo $split_variables[0] specified under section [Cellranger], line $line.\n";
            $variables{"Cellranger_cells"} = $med_variable;
        }
        else{
            $variables{"Cellranger_cells"} = $med_variable;
        } 
    }

# SEURAT SETTINGS

    elsif($read_settings =~ /^(Seurat_project)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\w/){
            $variables{"Seurat_project"} =$med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(Output_directory)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\/\w./){
            $variables{"Seurat_output_directory"} = $med_variable;
        }
        else{
            die "\n$split_variables[0] not specified, please check configuration under [Seurat], line $line...\n";
        }
    }
    elsif($read_settings =~ /^(Min_cells)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\d+/ && $med_variable > 0){
            $variables{"Seurat_min_cells"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(Min_features)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\d+/ && $med_variable > 0){
            $variables{"Seurat_min_features"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(Minimum_subset_features)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\d+/ && $med_variable > 0){
            $variables{"Seurat_minimum_subset_features"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(Maximum_subset_features)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\d+/ && $med_variable > 0){
            $variables{"Seurat_maximum_subset_features"} = $med_variable;
        }
        else{
            die "\nNo maximum subset features specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(Maximum_subset_mitodna)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\d+/ && $med_variable > 0){
            $variables{"Seurat_max_subset_mito_dna"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(Normalization_method)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\w/){
            $variables{"Seurat_normalization_method"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(Scale_factor)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\d+/ && $med_variable > 0){
            $variables{"Seurat_scale_factor"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(VarFeature_selection_method)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\w/){
            $variables{"Seurat_varfeature_selection_method"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(VarFeature_number)/){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\d+/ && $med_variable > 0){
            $variables{"Seurat_varfeature_number"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(Top_variable_features)/){
                @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\d+/ && $med_variable > 0){
            $variables{"Seurat_topvariable_features"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(Print_dimensions)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\d+/ && $med_variable > 0){
            $variables{"Print_dimensions"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(Print_features)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\d+/ && $med_variable > 0){
            $variables{"Print_features"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(VDL_dimensions)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\d+/ && $med_variable > 0){
            $variables{"VLD_dimensions"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(DHM_dimensions)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\d+/ && $med_variable > 0){
            $variables{"DHM_dimensions"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(DHM_cells)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\d+/ && $med_variable > 0){
            $variables{"DHM_cells"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(JackStraw_replicates)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\d+/ && $med_variable > 0){
            $variables{"Jackstraw_replicates"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(JackStraw_score_dimensions)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\d+/ && $med_variable > 0){
            $variables{"JackStraw_score_dimensions"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(JackStrawPlot_dimensions)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\d+/ && $med_variable > 0){
            $variables{"JackStrawPlot_dimensions"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(PC_chosen_Neighbors)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\d+/ && $med_variable > 0){
            $variables{"PC_chosen_Neighbors"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(cluster_resolution)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\d+/ && $med_variable > 0){
            $variables{"cluster_resolution"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    elsif($read_settings =~ /^(Reduction_method)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\w+/){
            $variables{"cluster_reduction_method"} = $med_variable;
        }
        else{
            die "\nNo $split_variables[0] specified under section [Seurat], line $line. Please check config. ini \n";
        }
    }
    else{}
}

close(SETTINGS); #We close the file after reading all the variables

#Check if all vars have been stored correctly and we pass them through the screen to see
#The execution parameters

print "\nThe execution parameters are...\n";
foreach my $keys (keys %variables){
    print $keys . "=" . $variables{$keys} . "\n";
}


#Now we load the links for all the refernece files for the genomes

open(LINKS, "links.cfg");

while($read_settings = <LINKS>){ #Although this is not the settings file, we are gonna recycle the $read_settings variable
    if($read_settings =~ /^(Cellranger_GRCh38)./){
        @split_variables = split(/\s=\s/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        #print $med_variable . "\n";
        $variables{"GRCH38_link"} = $med_variable;
        #and now we will save the name of the .tar.gz file
        @split_variables = split (/\//, $med_variable);
        foreach(@split_variables){
            if($_ =~ /.+\.tar\.gz"$/){
                chop $_;
                $_ =~ s/\.tar\.gz//g; #remove the .tar.gz from the file 
                $variables{"GRCH38_file"} = $_;
            }
            else{
            }
        }
    }
    elsif($read_settings =~ /^(Cellranger_hg19)./){
        @split_variables = split(/\s=\s/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        #print $med_variable . "\n";
        $variables{"hg19_link"} = $med_variable;
        @split_variables = split (/\//, $med_variable);
        foreach(@split_variables){
            if($_ =~ /.+\.tar\.gz"$/){
                chop $_;
                $_ =~ s/\.tar\.gz//g; #remove the .tar.gz from the file 
                $variables{"hg19_file"} = $_;
            }
            else{
            }
        }
    }
    elsif($read_settings =~ /^(Cellranger_mouse)./){
        @split_variables = split(/\s=\s/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        #print $med_variable . "\n";
        $variables{"mouse_link"} = $med_variable;
        @split_variables = split (/\//, $med_variable);
        foreach(@split_variables){
            if($_ =~ /.+\.tar\.gz"$/){
                chop $_; #remove the ""
                $_ =~ s/\.tar\.gz//g; #remove the .tar.gz from the file 
                $variables{"mouse_file"} = $_;
            }
            else{
            }
        }
    }
    else{}
}

close(LINKS);

#We get the date and time in which the analysis is being performed

my $date_and_time_1 = localtime();

#And now we can call the different modules on the sofware using the 
#Stored variables determined for each one of the parts

print "\n\nInitiating quality control analysis [$date_and_time_1] \n";

#QualityControl(%variables);

my $date_and_time_2 = localtime();

print "\n\nQuality control analysis ended at [$date_and_time_2].\n\n";
print "\n\nInitiating Cellranger analysis at [$date_and_time_2]\n\n";

$variables{"Cellranger_output"} = ChaosMatrixGeneration(%variables);

my $date_and_time_3 = localtime();
print "\n\nCellranger analysis ended at [$date_and_time_3]\n\n";
print "Output for Cellranger located at: $variables{'Cellranger_output'}";
print "\n\nInitiating Seurat based analysis of data [$date_and_time_3]\n\n";

SeuratAnalysis($current_directory, %variables);

#We assign this function to the output because we will be retourning the path to the output
#of cellranger to the key Cellranger_output on the variables hash. This way we can access
#this location on the subroutine for R analysis and perform the analysis at that position

