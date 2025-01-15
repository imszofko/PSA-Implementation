# NeedleWunsch and SmithWaterman

Two programmes created during the duration of the course 'Project in Bioinformatics.'

## Installation Instructions

The code can be downloaded from **https://github.com/imszofko/PSA-Implementation.git** or the **assignment page** in canvas. 
Once the zip folder is downloaded, you can unzip the files and save contents into a folder.
- Recommended you open the folder in the editor of your choice.
- Recommended you have **Biopython** package installed

The following files should be present in the folder:
- needleWunsch.py
- smithWaterman.py
- Variant700202.fasta
- Variant713680.fasta

## Running the Programmes

You open either the needleWunch.py file or the smithWaterman.py file and run the programme. 

Immediately you will be asked to enter an ENSEMBL sequence ID. In this case you will enter: **ENSG00000139618.**

After you press enter, it will ask you for the name/path to the FASTA file that is saved in the folder; you enter: **Variant700202.fasta** or **Variant713680.fasta**. 
- To use with different FASTA file, you just enter the name (and the path of in a different working directory) of the FASTA file.

The code will run on its own from here on out. It takes approximately 25 minutes for each programmes to run. 

## Output Expectations

When the programme is finished running, the terminal will read the message: **Highlighted alignment written to '_Alignment.txt**. 
- Depending on the programme, it might be NW_Alignment.txt or SW_Alignment.txt.

The alignment TXT files are located in the same folder. The files contain: the sequence ID, alignment, match string, and the highlighted regions. 
- The highlighted regions are denoted by the lowercase base pairs. These are base pairs that did not match in the alignment. The purpose of the feature is to make it easier to identify nucleotides that differ in an alignment.

