########################################################
#     MODULE FOR QUALITY CONTROL OF ALL FASTQ FILES    #
#       GENERATED BY 10X scRNA SEQUENCING              #
#       author: Gonzalo Vigara Astillero               #
########################################################

package Modules::QualityControl; #The package we have to load
use strict;
use warnings;
use Cwd qw(getcwd);

use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(QualityControl); #This is the module name

sub QualityControl{
    my %variables=@_; #We define the parameter input for the module
    #Then declare the variables we are going to use
    my $fastq_files;
    my @files_on_dir;
    my @fastq_files;
    my $analysis_list;

    #FIRST: we wan to check if FASTQC is installed on the system

    my $fastqc_path = qx/which fastqc/; #With qx we execute bash commands

    if($fastqc_path eq ""){
        die "\n FASTQC HAS NOT BEEN FOUND ON THE SYSTEM. PLEASE INSTALL IT AND TRY AGAIN. \n"
    }
    else{
        print "\n----- Initiating FASTQC analysis -------\n\n"; 
    }

    #Now we are going to set the working path to the value in $variables{"Working_directory"}
    chdir $variables{"Working_directory"}; #Change the directory to the one specified in the var WD
    print "Working directory has been set to: " . getcwd() . "\n"; #We check the working directory

    #After setting the working directory we are going to open it and check for FASTQ files
    #FASTQ files can be either. fastq or fastq.gz

    opendir(WORKINGDIR, $variables{"Working_directory"}) || die "Unable to open" . $variables{"Working_directory"}. "\n";
    #Now we read the files contained in the directory and list them on the output
    @files_on_dir = readdir(WORKINGDIR);
    foreach my $i (@files_on_dir){ #We get the FASTQ files on the directory
        print $i . "\n";
        if($i =~ /.*\.fastq$/ or $i =~ /.*\.fastq\.gz$/ or $i =~ /.*\.fq$/ or $i =~ /.*\.fq\.gz$/ or $i=~ /.*\.fastq\.lane\.clean$/ or $i=~ /.*\.fq\.bz2$/ or $i=~ /.*\.fastq\.bz2$/){
            #print $i . "\n";
            push(@fastq_files, $i);
        }
    }
    closedir(WORKINGDIR);

    #Now we check if the Output directory exists. If not, we create it. 
    if(-d $variables{"Output_QC"}){
        print "\n" . $variables{"Output_QC"} . " already exists. \n";
    }
    

    #We check if there are fastq files on the list. If not, the program stops. 
    if(scalar(@fastq_files) == 0){
        die "\nTHERE ARE NO FASTQ FILES ON THE DIRECTORY. PLEASE CHECK CONFIG.INI \n";
    }
    else{
        print "\nInitiating FASTQC analysis\n\n";
        print "\nThe FASTQ files contained in " . $variables{"Working_directory"} . " are: \n\n";
        foreach my $i (@fastq_files){
            print "- " . $i . "\n";
            #qx/fastqc -o $variables{"Output_QC"} -t $variables{"Threads"} -f fastq --noextract  $i/;
        }
        print "\n";
    }

    #Now we perform the analysis on all the files without limit 
    $analysis_list = join(' ', @fastq_files); #Concatenate all the fastq files for analysis
    #print $analysis_list;
    #Now we execute the FASTQ quality analysis, previous check if there are already analysis files

    if(-d $variables{"Output_QC"}){ #If the output directory exists on
        print "\n" . $variables{"Output_QC"} . " already exists. \n";
        chdir($variables{"Output_QC"}); #Change to the analysis directory
        #Now we use glob to check all the files on the directory with the extension .html (output for fastqc)
        my @QC_files_on_dir = glob "*.html";
        if(length(@QC_files_on_dir != 0)){  #If there are files on the directory that correspond to html
            print "\nFASTQC analysis already performed for files on " . $variables{"Working_directory"} . ". Skipping...\n";
        }
        else{ #if the file is empty of html files we have to perform the analysis for QC control
           chdir($variables{"Working_directory"}); #Change again to the working directory
           qx/fastqc $analysis_list -o $variables{"Output_QC"} -t $variables{"QC_threads"} --noextract/;
        }
    }
    else{
        print "\nCreating output directory for FASTQC analysis... \n\n"; 
        qx/mkdir $variables{"Output_QC"}/; #We create the output directory
        qx/fastqc $analysis_list -o $variables{"Output_QC"} -t $variables{"QC_threads"} --noextract/; #And execute the analysis
    }
  
    print "----- QC ANALYSIS FINISHED -----";

    ####################
    # MULTIQC ANALYSIS #
    ####################

    if($variables{"Multi_QC"} eq "yes" or $variables{"Multi_QC"} eq "Yes"){
        #MultiQC generates a single report file from different FASTQC files, allowing a fast view of all of them

        #We need first the username of whom it's executing the script

        my $logname = $ENV{LOGNAME} || $ENV{USERNAME} || $ENV{USER};
        #print "The user is $logname";

        #If MULTIQC has been already used, it should be on the path stored on $multiqc_path
        my $multiqc = "/home/".$logname."/.local/bin/multiqc";
        my $multiqc_git = "git clone https://github.com/ewels/MultiQC.git";

        if(-e $multiqc){ #If it has been already installed and exists on the path
                print "\n\n----- MULTIQC FILE GENERATION -----\n";
                #We change the working directory to the one containing the QC analysis HTML file
                chdir $variables{"Output_QC"};
                #And we exectue the multiqc analysis using the installed variable
                qx/$multiqc ./;
        }
        else{
            my $python3_pip = qx/which pip3/; #Check if python3-pip is installed on the system
            if($python3_pip eq ""){
            die "\nPython 3 pip is not installed. Please install it before \n
            generating multiqc file; (Ubuntu: sudo apt-get install python3-pip)\n\n";}
            else{
                #print $multiqc_path;
                print "\nMultiQC not found on the system, installing...\n";
                qx/$multiqc_git/;
                qx/cd ".\/MultiQC" && pip3 install \./;
                qx/rm -rf .\/MultiQC/;
                print "\n\n----- MULTIQC FILE GENERATION -----\n";
                #We change the working directory to the one containing the QC analysis HTML file
                chdir $variables{"Output_QC"};
                #And we exectue the multiqc analysis using the installed variable
                qx/$multiqc ./;
            }
        }

        print "\nMultiQC containing all analysis generated.\n";
    }
    elsif($variables{"Multi_QC"} eq "no" or $variables{"Multi_QC"} eq "No"){
        print "\nThe user has decided no to perform a MultiQC analysis. Skipping... \n";
    }
    else{
        die "\nMultiQC analysis value on config.ini incorrect. Please check config.ini and try again...\n"
    }
    return "OK";
}