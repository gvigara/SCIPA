########################################################
#   MODULE FOR SEURAT ANALYSIS OF Al PROCESSED FILES   #
#       GENERATED BY 10X scRNA SEQUENCING              #
#       author: Gonzalo Vigara Astillero               #
########################################################

package Modules::Seurat_analysis; #The package we have to load
use strict;
use warnings;
use Cwd qw(getcwd);
use Statistics::R;

use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(SeuratAnalysis); #This is the module name

sub SeuratAnalysis{
    #We retrieve the arguments passed to the subroutine
    my %variables = @_;
    my $current_directory = getcwd;
    #For the analysis we create a directory in the current directory with the
    #directory name specified in the variable

      
    my $output_directory = $variables{"Working_directory"}. "/Seurat_analysis\-". $variables{'Actual_date'};

    #We open the log file

    my $log_name = $variables{'Working_directory'}."/Seurat_Log_".$variables{'Actual_date'}.".txt";
    open(QCLOG, ">".$log_name);

    #Create the output directory

    print "\nGenerating output directory = $output_directory \n";
    print QCLOG "\nGenerating output directory = $output_directory \n";
    system(qq'mkdir $output_directory');

    #Check the current directory from where the script is being executed
    my $R_modules_directory = $variables{"Binaries_directory"} . "/Rscripts";

    print "The current directory for R modules is: ".$R_modules_directory . "\n";
    print QCLOG "The current directory for R modules is: ".$R_modules_directory . "\n";

    print "\nStarting R instance for analysis...\n";
    print QCLOG "\nStarting R instance for analysis...\n";

    #Start an R instance and execute everything
    my $R = Statistics::R->new();

    print ("Checking dependencies:\n- RCurl\n- httr\n- plotly\n- Seurat\n- patchwork\n- dply\n- scCATCH \n- limma");
    print QCLOG ("Checking dependencies:\n- RCurl\n- httr\n- plotly\n- Seurat\n- patchwork\n- dply\n- scCATCH \n- limma");

    my $output_dependencies = $R -> run_from_file("$R_modules_directory/Dependencies.R");

    print QCLOG "\n\n" . $output_dependencies . "\n\n";

    print("\nDependencies...OK; loading functions! \n");
    print QCLOG ("\nDependencies...OK; loading functions! \n");

    my $output_functions = $R -> run_from_file("$R_modules_directory/Scripts_functions.R");

    print QCLOG "\n\n $output_functions \n\n";
    print "\n\n $output_functions \n\n";

    #We add inverted commas to the needed variables, input directory and output directory

    my $Seurat_input_directory = "\"$variables{'Cellranger_output'}\"";
    my $Seurat_output_directory = "\"$output_directory\"";
    my $Seurat_project_name = "\"$variables{'Seurat_project'}\"";
    my $Seurat_normalization_method = "\"$variables{'Normalization_method'}\"";
    my $Seurat_varfeature_selection = "\"$variables{'VarFeature_selection_method'}\"";
    my $Seurat_reduction_method = "\"$variables{'Reduction_method'}\"";
    my $scCATCH_organism = "\"$variables{'Cellranger_organism'}\"";
    #my $scCATCH_tissue = "\"$variables{'scCATCH_tissue'}\"";


    print "\nEXECUTING SEURAT ANALYSIS...\n";
    print QCLOG "\nEXECUTING SEURAT ANALYSIS...\n";

    chdir($output_directory); #We change to the output directory for all the files

    #wE execute each of the analysis subroutines individually and print on screen each of the outputs
    my $setup_and_qc = $R -> run("seurat_object <- setup_qc($Seurat_input_directory, $Seurat_output_directory, $Seurat_project_name,
                                $variables{'Min_cells'}, $variables{'Min_features'}, $variables{'Minimum_subset_features'}, 
                                $variables{'Maximum_subset_features'},$variables{'Maximum_subset_mitodna'}, $scCATCH_organism)");

    print "Setup and quality control of the Seurat object... Done! - ". localtime() ."\n"; 
    print "\n$setup_and_qc\n";
    print QCLOG "\n$setup_and_qc\n";

    my $normalization = $R -> run("normalized_object <- Normalization(seurat_object, $Seurat_normalization_method, $variables{'Scale_factor'}, 
                                                        $Seurat_varfeature_selection, $variables{'VarFeature_number'},
                                                        $variables{'Top_variable_features'},  $variables{'QC_threads'})");
    
    print "Normalization of the Seurat object... Done! - ". localtime() ."\n"; 
    print QCLOG "\n$normalization\n";
    print "\n$normalization\n";

    my $PCA = $R -> run("reduced_object <- PCA(normalized_object, $variables{'Print_dimensions'}, 
                                                                $variables{'Print_features'}, $variables{'VDL_dimensions'}, 
                                                                $variables{'DHM_dimensions'}, $variables{'DHM_cells'}, $variables{'QC_threads'})");

    print "Principal component analysis of the Seurat object... Done! - ". localtime() ."\n"; 
    print QCLOG "\n$PCA\n";
    print "\n$PCA\n";

    my $dimensionality = $R -> run("dimensioned_object <- dimensionality(reduced_object, $variables{'JackStraw_replicates'}, 
                                                            $variables{'JackStraw_score_dimensions'}, $variables{'JackStrawPlot_dimensions'}, $variables{'QC_threads'})");

    print "Dimensionality calculation of the Seurat object... Done! - ". localtime() ."\n"; 
    print QCLOG "\n$dimensionality\n";
    print "\n$dimensionality\n";

    my $clustering = $R -> run("clustered_object <- clustering_and_NLR(dimensioned_object, $variables{'PC_chosen_Neighbors'}, 
                                                    $variables{'cluster_resolution'}, $Seurat_reduction_method,  $variables{'QC_threads'})");

    print "Clustering of the Seurat object... Done! - ". localtime() ."\n"; 
    print  QCLOG"\n$clustering\n";
    print  "\n$clustering\n";
    
    my $scCATCH_annotation = $R -> run ("scCATCH_annot_object <- scCATCH_annotation(clustered_object, $scCATCH_organism, $variables{'scCATCH_tissue'}, $variables{'Violin_plot_cluster_genes'}, $variables{'Heatmap_cluster_genes'}, $variables{'QC_threads'})");
    
    print "Annotation of clusters with scCATCH... Done! - ". localtime() ."\n"; 
    print QCLOG "\n$scCATCH_annotation\n";
    
    my $cluster_typing = $R -> run ("annotation_clusters <- cluster_marking(clustered_object, scCATCH_annot_object, $Seurat_reduction_method)");
    
    print "Typing of clusters... Done! - ". localtime() ."\n"; 
    print QCLOG "\n$cluster_typing\n";

    close QCLOG;
    #print "\n$seurat_analysis\n"; 
}
