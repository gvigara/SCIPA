#!/usr/bin/perl

use strict;
use warnings;
use Cwd qw(getcwd); 

#######################################
# CELLRANGER ANALYSIS FOR SCRNA FASTQ #
#######################################

#First we declare the variables

my %variables;
my $read_settings;
my @split_variables; #intermediate list for splitted strings
my $med_variable; #intermediate variable for splitting strings
my $human_genome_selection; #variable storing the option the user chooses for human genome
my $reference_folder; #folder where the referecene genome is gonna be downloaded
my $cellranger; #stores the cellranger command
my $localmem;
my $localcores;

#Now we open the config file and read the settings for cellranger 
#The settings will be read from the main pl when this programms are
#configured as Perl Modules

open(SETTINGS, "config.ini");
while($read_settings=<SETTINGS>){
    if($read_settings =~ /^(Working_directory)./){ #This is the path to FASTQ files
        @split_variables = split(/\s=\s/, $read_settings);
        if($split_variables[1] =~ /\/\w/){ #If there is a directory specified, save it to variables hash
            $med_variable = $split_variables[1];
            chomp($med_variable); 
            $variables{"Working_directory"} = $med_variable;
        }
        else{
            #In case there is nothing written under that section, the program dies with this error message
            die "Working directory under [General] section not specified or incorrect. \n";
        }
        #print $split_variables[1];   
    }
    elsif($read_settings =~ /^(Cellranger_path)./){ #We get the path where cellranger is installed
        @split_variables = split(/\s=\s/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable eq ""){ #With this loop we check if the parameter is specified. If not, the program stops. 
            die "$split_variables[0] parameter not specified. Please check config.ini";
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
            die "$split_variables[0] parameter not specified. Please check config.ini";
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
            die "$split_variables[0] parameter not specified. Please check config.ini";
        }
        elsif($med_variable =~ /[H-h]uman/ || $med_variable =~ /[M-m]ouse/){ #Only for human or mouse, if not, program stops. 
            $variables{"Cellranger_org"} = $med_variable;
        }      
        else{
            die "\n Organism not recognized. Please check config.ini \n";
        } 
    }elsif($read_settings =~ /^(Cellranger_reference)./){ #Reference genome. If empty, will be downloaded. 
        @split_variables = split(/\s=\s/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable eq ""){
            print "\n $split_variables[0] parameter not specified.\n";
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
            print "\n $split_variables[0] parameter not specified. \n";
            $variables{"Cellranger_cells"} = $med_variable;
        }
        else{
            $variables{"Cellranger_cells"} = $med_variable;
        } 
    }
    elsif($read_settings =~ /^(Cellranger_cores)./){ #Number of CPU cores. If not specified, program will choose
        @split_variables = split(/\s=\s/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable eq ""){ #In this case, this parameter is not mandatory. Thus, if it has no input, the program will continue. 
            print "\n $split_variables[0] parameter not specified.\n";
            $variables{"Cellranger_cores"} = "AUTOMATIC SELECTION BY SOFTWARE";
        }
        else{
            $variables{"Cellranger_cores"} = $med_variable;
        } 
    }
    elsif($read_settings =~ /^(Cellranger_mem)./){
        @split_variables = split(/\s=\s/, $read_settings);
        $med_variable = $split_variables[1];
        chomp($med_variable);
        if($med_variable eq ""){ #In this case, this parameter is not mandatory. Thus, if it has no input, the program will continue. 
            print "\n $split_variables[0] parameter not specified.\n";
             $variables{"Cellranger_mem"} = "AUTOMATIC SELECTION BY SOFTWARE";
        }
        else{
            $variables{"Cellranger_mem"} = $med_variable;
        } 
    }
    else{}
}

close(SETTINGS);

#We check that the input of all parameters is correct

print "\n----- CELLRANGER PARAMETERS -----\n\n";

print "Working directory: " . $variables{"Working_directory"} . "\n";
print "Cellranger installation path: " . $variables{"Cellranger_path"} . "\n";
print "ID for the cellranger project: " . $variables{"Cellranger_id"} . "\n";
print "Organism: " . $variables{"Cellranger_org"} . "\n";
print "Path to reference genome: " . $variables{"Cellranger_ref"} . "\n";
print "Expected cells for the analysis: " . $variables{"Cellranger_cells"} . "\n";
print "CPU cores for the analysis: " . $variables{"Cellranger_cores"} . "\n";
print "RAM memory assigned for the analysis: " . $variables{"Cellranger_mem"} . "\n";

#If Path to reference genome is empty, we will print a message asking for the user to download it.
#First we want to check wether wget or curl is present on the system

my $wget_path = qx/which wget/;
#print "\nWget Path = " . $wget_path;
my $curl_path = qx/which curl/;
#print "\nCurl Path = " . $curl_path;

#And we locate the system on the working directory
chdir $variables{"Cellranger_path"};

GENOME_SELECTION: if($variables{"Cellranger_ref"} eq ""){ #If there is no input on the genome, we use this taggged loop
    $reference_folder = $variables{"Cellranger_path"} . "/Reference_genomes";
    mkdir($reference_folder); #We generate a new reference folder
    chdir $reference_folder; #We move to the reference folder
    print "\nREFERENCE GENOME MISSING, Downloading...\n";
    if($variables{"Cellranger_org"} =~ /[H-h]uman/){
        print "\nPlease, indicate which human reference genome you want to be downloaded:\n";
        print "\n1) Human GRCh38\n2) Human hg19\n\n";
        print "Specify 1 or 2: ";
        $human_genome_selection = <STDIN>; #we allow the user to choose
        if(length($wget_path) != 0){ #If wget is absent, the system will download the genome with curl
            if ($human_genome_selection == 1){ #If the first option is selected, we donwload the selected genome
                if(-d "./refdata-cellranger-GRCh38-3\.0\.0"){#This checks if the process has been already done although
                #the user didnt specify it on the options file
                    $variables{"Cellranger_ref"} = getcwd() . "/refdata-cellranger-GRCh38-3\.0\.0";
                    print "\nThe selected genome has already been downloaded and decompressed on: " . $variables{"Cellranger_ref"} . "\n";
                    print "\nSkipping...\n";
                }
                else{
                    print "\nDownloading Human GRCh38 genome...\n";
                    #We download the genome. -C is for continue the donwload the file
                    #In case an incomplete download took place, this time the genome will donwload again. 
                    qx/wget -c http:\/\/cf\.10xgenomics\.com\/supp\/cell-exp\/refdata-cellranger-GRCh38-3\.0\.0\.tar\.gz/;
                    print "\nDecrompressing the genome\n";
                    #Uncompress the genome in the predefined folder
                    qx/tar -xzvf refdata-cellranger-GRCh38-3\.0\.0\.tar\.gz/;
                    #Stablish the reference genome variable to the current folder
                    $variables{"Cellranger_ref"} = getcwd() . "/refdata-cellranger-GRCh38-3\.0\.0";
                    print "\n---------------------------------------------------------------------------------------------\n";
                    print "\nThe reference genome for the analysis has been stored on:\n";
                    print $variables{"Cellranger_ref"};
                }
            }
            elsif ($human_genome_selection == 2){
                if(-d "./refdata-cellranger-hg19-3\.0\.0"){
                    $variables{"Cellranger_ref"} = getcwd() . "/refdata-cellranger-hg19-3\.0\.0";
                    print "\nThe selected genome has already been downloaded and decompressed on: " . $variables{"Cellranger_ref"} . "\n";
                    print "\nSkipping...\n";
                }
                else{
                    print "\nDownloading Human hg19 genome...\n";
                    #We download the genome
                    qx/wget -c http:\/\/cf\.10xgenomics\.com\/supp\/cell-exp\/refdata-cellranger-hg19-3\.0\.0\.tar\.gz/;
                    print "\nDecrompressing the genome\n";
                    #Uncompress the genome in the predefined folder
                    qx/tar -xzvf refdata-cellranger-hg19-3\.0\.0\.tar\.gz/;
                    #Stablish the reference genome variable to the current folder
                    $variables{"Cellranger_ref"} = getcwd() . "/refdata-cellranger-hg19-3\.0\.0";
                    print "\n---------------------------------------------------------------------------------------------\n";
                    print "\nThe reference genome for the analysis has been stored on:\n";
                    print $variables{"Cellranger_ref"};
                }
            }
            else{
                print "\nInvalid option, please select again...\n";
                goto GENOME_SELECTION; #In the case we select other genome, we go back to the human selection loop
            }
        }
        elsif(length($curl_path) != 0){
            if ($human_genome_selection == 1){ #If the first option is selected, we donwload the selected genome
                if(-d "./refdata-cellranger-GRCh38-3\.0\.0"){
                    $variables{"Cellranger_ref"} = getcwd() . "/refdata-cellranger-GRCh38-3\.0\.0";
                    print "\nThe selected genome has already been downloaded and decompressed on: " . $variables{"Cellranger_ref"} . "\n";
                    print "\nSkipping...\n";
                }
                else{
                    print "\nDownloading Human GRCh38 genome...\n";
                    #We download the genome
                    qx/curl -O "http:\/\/cf\.10xgenomics\.com\/supp\/cell-exp\/refdata-cellranger-GRCh38-3\.0\.0\.tar\.gz"/;
                    print "\nDecrompressing the genome\n";
                    #Uncompress the genome in the predefined folder
                    qx/tar -xzvf refdata-cellranger-GRCh38-3\.0\.0\.tar\.gz/;
                    #Stablish the reference genome variable to the current folder
                    $variables{"Cellranger_ref"} = getcwd() . "/refdata-cellranger-GRCh38-3\.0\.0";
                    print "\n---------------------------------------------------------------------------------------------\n";
                    print "\nThe reference genome for the analysis has been stored on:\n";
                    print $variables{"Cellranger_ref"};
                }
            }
            elsif ($human_genome_selection == 2){
                if(-d "./refdata-cellranger-hg19-3\.0\.0"){
                    $variables{"Cellranger_ref"} = getcwd() . "/refdata-cellranger-hg19-3\.0\.0";
                    print "\nThe selected genome has already been downloaded and decompressed on: " . $variables{"Cellranger_ref"} . "\n";
                    print "\nSkipping...\n";
                }
                else{
                    print "\nDownloading Human hg19 genome...\n";
                    #We download the genome
                    qx/curl -O "http:\/\/cf\.10xgenomics\.com\/supp\/cell-exp\/refdata-cellranger-hg19-3\.0\.0\.tar\.gz"/;
                    print "\nDecrompressing the genome\n";
                    #Uncompress the genome in the predefined folder
                    qx/tar -xzvf refdata-cellranger-hg19-3\.0\.0\.tar\.gz/;
                    #Stablish the reference genome variable to the current folder
                    $variables{"Cellranger_ref"} = getcwd() . "/refdata-cellranger-hg19-3\.0\.0";
                    print "\n---------------------------------------------------------------------------------------------\n";
                    print "\nThe reference genome for the analysis has been stored on:\n";
                    print $variables{"Cellranger_ref"};
                }
            }
            else{
                print "\nInvalid option, please select again...\n";
                goto GENOME_SELECTION; #In the case we select other genome, we go back to the human selection loop
            }
        }
        else{
            die "\nCannot find Curl nor Wget, please download them and try again!\n"
        }
    }
    elsif($variables{"Cellranger_org"} =~ /[M-m]ouse/){
        if(-d "./refdata-cellranger-mm10-3\.0\.0"){
            $variables{"Cellranger_ref"} = getcwd() . "/refdata-cellranger-mm10-3\.0\.0";
            print "\nThe selected genome has already been downloaded and decompressed on: " . $variables{"Cellranger_ref"} . "\n";
            print "Skipping...\n";
        }
        else{
            print "\nDownloading mouse reference genome...\n";
            if(length($wget_path) != 0){ #If wget is absent, the system will download the genome with curl
            qx/curl -O "http:\/\/cf\.10xgenomics\.com\/supp\/cell-exp\/refdata-cellranger-mm10-3\.0\.0\.tar\.gz/;
            #Now we uncompress the genome in the predefined folder
            qx/tar -xzvf refdata-cellranger-mm10-3\.0\.0\.tar\.gz/;
            #STablish the reference genome variable to the current folder
            $variables{"Cellranger_ref"} = getcwd() . "/refdata-cellranger-mm10-3\.0\.0";
            }
        }
    }
    else{
        die "The genome for the organism " . $variables{"Cellranger_org"} .  " is not avaiable. Exiting...";
    }
}
elsif($variables{"Cellranger_ref"} =! ""){ #In case the reference genome path is already specified on the config.ini file
    print "\n---------------------------------------------------------------------------------------------\n";
    print "\nThe reference genome for the analysis is stored on:\n";
    print $variables{"Cellranger_ref"};
}
else{} #Close the loop 



#After we have our reference genome we just have to launch the Cellranger Analysis

print "\n-------------------------------------------------------------------------------------------\n";
print "\t\tINITIATING CELLRANGER ANALYSIS";
print "\n-------------------------------------------------------------------------------------------\n";

#Now we change the directory back to our working directory

chdir $variables{"Working_directory"};

#And then we use all the variables to execute the programm

$cellranger = "cellranger count " . "--id=" . $variables{"Cellranger_id"} . " --transcriptome=" . $variables{"Cellranger_ref"}
            . " --fastqs=" . $variables{"Working_directory"} . " --expect-cells=" . $variables{"Cellranger_cells"};

#now we add the flags for the number of CPU cores and memory if they were selected. If not
#The flag will not appear on the command

$localmem = $variables{"Cellranger_mem"};
$localcores = $variables{"Cellranger_cores"}; #We assign the hash to variables to concatenate them, if not
#the concatenation will return a "position" rather than the value

#Remember when cores or memory values are not assigned by user the value is filled with the string
#"AUTOMATIC SELECTION BY SOFTWARE"

if($variables{"Cellranger_cores"} ne "AUTOMATIC SELECTION BY SOFTWARE" && $variables{"Cellranger_mem"} ne "AUTOMATIC SELECTION BY SOFTWARE"){
    $cellranger = $cellranger . " --localcores=" . $localcores . " --localmem=" . $localmem;
}
elsif($variables{"Cellranger_cores"} ne "AUTOMATIC SELECTION BY SOFTWARE" && $variables{"Cellranger_mem"} eq "AUTOMATIC SELECTION BY SOFTWARE"){
    $cellranger = $cellranger . " --localcores=" . $localcores;
}
elsif($variables{"Cellranger_mem"} ne "AUTOMATIC SELECTION BY SOFTWARE" && $variables{"Cellranger_cores"} eq "AUTOMATIC SELECTION BY SOFTWARE"){
    $cellranger = $cellranger . " --localmem=" . $localmem;
}
else{}

#And now we execute the command!!

print "\nExecuting:\n" . $cellranger . "\n";

system($variables{"Cellranger_path"} . "/" . $cellranger);

