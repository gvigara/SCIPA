#!/usr/bin/perl

use strict;
use warnings;
use lib '/home/gvigast/TFM/TFM/Scripts';

BEGIN{
    die "Your version of perl is old, please update it ASAP! \n" unless
    perl -> VERSION => 5.30;
}

#We load the modules required for the analysis

use Modules::QualityControl;
use Modules::ChaosMatrixGenerator;

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
    elsif($read_settings =~ /^(Core_threads)./){ #Read the core threads designated for analysis
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
    elsif($read_settings =~ /^(RAM_memory)./){ #Read the RAM memory amount designated for analysis
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
    elsif($read_settings =~ /^(Analysis_files_path)./){
        @split_variables = split(/\s?=\s?/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable =~ /\/\w./){
            $variables{"Previous_analysis_seurat"} = $med_variable;
        }
        else{
            print "\n$split_variables[0] not specified, initializing de novo analysis...\n";
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
            $variables{"Seurat_max_subset_mito_dna"} = $med_variable;
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
}

close(SETTINGS); #We close the file after reading all the variables

#Check if all vars have been stored correctly and we pass them through the screen to see
#The execution parameters

print "\nThe execution parameters are...\n";
foreach my $keys (keys %variables){
    print $keys . "=" . $variables{$keys} . "\n";
}

#And now we can call the different modules on the sofware using the 
#Stored variables determined for each one of the parts

#QualityControl(%variables);

ChaosMatrixGeneration(%variables);