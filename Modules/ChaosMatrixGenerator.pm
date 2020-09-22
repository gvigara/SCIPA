########################################################
# MODULE FOR GENERATION OF THE CHAOS MATRIX FROM FASTQ #
#       GENERATED BY 10X scRNA SEQUENCING              #
#       author: Gonzalo Vigara Astillero               #
########################################################

package Modules::ChaosMatrixGenerator;
use strict;
use warnings;
use Cwd qw(getcwd);

use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(ChaosMatrixGeneration); #This is the module name

sub ChaosMatrixGeneration{
    my %variables = @_;

    #We define the individual variables we need for the use of the software

    my $human_genome_selection; #variable storing the option the user chooses for human genome
    my $reference_folder; #folder where the referecene genome is gonna be downloaded
    my $cellranger; #stores the cellranger command
    my $localmem;
    my $localcores;

    #And begin with the proper analysis script

    print "\n----- CELLRANGER PARAMETERS -----\n\n";

    print "Working directory: " . $variables{"Working_directory"} . "\n";
    print "Cellranger installation path: " . $variables{"Cellranger_path"} . "\n";
    print "ID for the cellranger project: " . $variables{"Cellranger_id"} . "\n";
    print "Organism: " . $variables{"Cellranger_org"} . "\n";
    print "Path to reference genome: " . $variables{"Cellranger_ref"} . "\n";
    print "Expected cells for the analysis: " . $variables{"Cellranger_cells"} . "\n";
    print "CPU cores for the analysis: " . $variables{"Core_threads"} . "\n";
    print "RAM memory assigned for the analysis: " . $variables{"RAM_memory"} . "\n";

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
    else{
        die "\nUnable to locate the reference genome: $variables{'Cellranger_ref'}. Please, check it and execute the software again. \n"
    } #Close the loop 



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

    $localmem = $variables{"RAM_memory"};
    $localcores = $variables{"Core_threads"}; #We assign the hash to variables to concatenate them, if not
    #the concatenation will return a "position" rather than the value

    #Remember when cores or memory values are not assigned by user the value is filled with the string
    #"AUTOMATIC SELECTION BY SOFTWARE"

    if($variables{"Core_threads"} ne "AUTOMATIC SELECTION BY SOFTWARE" && $variables{"RAM_memory"} ne "AUTOMATIC SELECTION BY SOFTWARE"){
        $cellranger = $cellranger . " --localcores=" . $localcores . " --localmem=" . $localmem;
    }
    elsif($variables{"Core_threads"} ne "AUTOMATIC SELECTION BY SOFTWARE" && $variables{"RAM_memory"} eq "AUTOMATIC SELECTION BY SOFTWARE"){
        $cellranger = $cellranger . " --localcores=" . $localcores;
    }
    elsif($variables{"RAM_memory"} ne "AUTOMATIC SELECTION BY SOFTWARE" && $variables{"Core_threads"} eq "AUTOMATIC SELECTION BY SOFTWARE"){
        $cellranger = $cellranger . " --localmem=" . $localmem;
    }
    else{
        print STDERR "\nUnable to determine the amount of RAM and number of cores specified\n";
    }

    #And now we execute the command!!

    print "\nExecuting:\n" . $cellranger . "\n";

    system($variables{"Cellranger_path"} . "/" . $cellranger);
}