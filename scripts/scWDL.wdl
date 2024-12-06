version 1.0

workflow scRNA_seq_analysis {
  input {
    File healthy_file
    File mild_file
    File severe_file
  }

  call RunAnalysis {
    input:
      rscript = "scRNA_seq_analysis.R",
      healthy_file = healthy_file,
      mild_file = mild_file,
      severe_file = severe_file
  }

  output {
    Array[File] plots = RunAnalysis.plots
    File best_markers_csv = RunAnalysis.best_markers_csv
  }
}

task RunAnalysis {
  input {
    File rscript
    File healthy_file
    File mild_file
    File severe_file
  }
  command {
    Rscript ${rscript} ${healthy_file} ${mild_file} ${severe_file} # Run the R script with the input file
  }


  output {
    Array[File] plots = glob("*.png")
    File best_markers_csv = "best_markers.csv"
  }
}