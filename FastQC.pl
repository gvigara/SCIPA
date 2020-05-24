#!/usr/bin/perl

# FASTQC QUALITY ANALYSIS SCRIPT

use strict;
use warnings;
use Cwd qw(getcwd); #This function will tell us the current working directory

my %variables;
my $read_settings;
my @split_variables;
my $med_variable;
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

#Open the configuration file to obtain the variables
open(SETTINGS,"config.ini");
while($read_settings=<SETTINGS>){
    if($read_settings =~ /^(Working_directory)./){
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
    elsif($read_settings =~ /^(Output_QC)./){
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
        else{
            $variables{"Output_QC"} = $variables{"Working_directory"};
        }
    }
    elsif($read_settings =~ /^(QC_threads)./){ #Get the number of threads for QC analysis
        @split_variables = split(/\s=\s/, $read_settings);
        print "Core threads for analysis: " . $split_variables[1] . "\n\n";
        $variables{"Threads"} = $split_variables[1];
    }
}
close(SETTINGS);

#print "Working path containing fastqc files: " . $variables{"Working_directory"}. "\n";
#print "Output path for FASTQC analysis: " . $variables{"Output_QC"} . "\n";

#Now we are going to set the working path to the value in $variables{"Working_directory"}
#print getcwd() . "\n"; #The directory we are in
chdir $variables{"Working_directory"}; #Change the directory to the one specified in the var WD
print "Working directory has been set to: " . getcwd() . "\n"; #We check the working directory

#After setting the working directory we are going to open it and check for FASTQ files
#FASTQ files can be either. fastq or fastq.gz

opendir(WORKINGDIR, $variables{"Working_directory"}) || die "Unable to open" . $variables{"Working_directory"}. "\n";
#Now we read the files contained in the directory and list them on the output
@files_on_dir = readdir(WORKINGDIR);
foreach my $i (@files_on_dir){ #We get the FASTQ files on the directory
    if(-f $i && $i =~ /\w\.fastq\.?/){ # if it is a file and it's a fastq file
            #print $i . "\n";
            push(@fastq_files, $i);
    }
    else{}
}
#Now we check if the Output directory exists. If not, we create it. 
if(-d $variables{"Output_QC"}){
    print "\n" . $variables{"Output_QC"} . " already exists. \n";
}
else{
    print "\nCreating output directory for FASTQC analysis... \n\n";
    qx/mkdir $variables{"Output_QC"}/;
}
closedir(WORKINGDIR);

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

qx/fastqc $analysis_list -o $variables{"Output_QC"} -t $variables{"Threads"} -f fastq/;

print "----- QC ANALYSIS FINISHED -----";

####################
# MULTIQC ANALYSIS #
####################

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

print "\nMultiQC containing all analysis generated.\n";


