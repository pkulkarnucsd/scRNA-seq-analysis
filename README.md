# scRNA-seq-analysis
Single-cell RNA-seq analysis using Cromwell and WDL

# Single-cell RNA-seq Analysis with Cromwell

This project demonstrates single-cell RNA-seq analysis using Cromwell with WDL, analyzing BALF samples from a healthy control, and patients with mild and severe Covid 19. The workflow is defined in the `scWDL.wdl` file and can be executed with the provided `inputs.json` file.

## Files:
- `scRNA_seq_analysis.R`: R script for analyzing single-cell RNA-seq data.
- `scWDL.wdl`: Workflow description file for Cromwell.
- `inputs.json`: JSON file containing input parameters for the Cromwell workflow.
- `data/`: Contains the input `.h5` files for the analysis.
- `outputs/`: Contains the output generated by Cromwell after running the workflow.

## Instructions:

1. **Set up Cromwell**:
   - Make sure you have Cromwell Downloaded
    
2. **Run the Analysis**:
   - To run the WDL workflow, use the following command:
     ```bash
     java -jar cromwell-85.jar run scWDL.wdl --inputs inputs.json
     ```
