# GutReferenceSet

Installation:
make sure you have all required external softwares installed, and fill the path to the relevant executables in:
GutReferenceSet/Utils/config.py

fill in the working directory (or the Demo directory) in:
GutReferenceSet/Build_Species_Set/config.py

if trying the Demo - correct the index of GutReferenceSet/Demo/full_metadata.csv to the full path of your system

Run and tested on:
python 3.7.4 on CentOS Linux 7.9

Creating a new reference set:

1. Create a list of files and their metadata
   A file named full_metadata.csv with index witch is the full path of an assembled genome fasta file
   and at least the following columns:
   Source, Method, AssemblyName, SampleName, RegistrationCode, DoNotTake
   Where:?
   Source - source of data, so that Source + RegistrationCode is a unique identifier of an individual from which
            assemblies where created
   Method - assembly creation method (MAG / isolate / nanopore...)
   AssemblyName - a unique identifier of the assembly
   SampleName - an identifier of the sample the assembly was created from (so as to identify
                assemblies originating from the same sample)
   RegistrationCode - an identifier of the individual the sample was taken from (so as to identify
                      assemblies originating from the same individual)
   DoNotTake - a columns which is either empty or includes reason not to consider the assembly

2. Calculate the qualities of the genomes by running checkM:
   GutReferenceSet/Build_Species_Set/quality.py

3. Filter by quality and by criteria of not taking same species from same person twice:
   GutReferenceSet/Build_Species_Set/filter.py

4. Mash all vs all into a memory mapped file using:
   GutReferenceSet/Build_Our_Set/distance.py

5. Cluster based of memory mapped distances. Using:
   GutReferenceSet/Build_Species_Set/hierarchical_clustering.py

6. Choose representatives:
   GutReferenceSet/Build_Species_Set/choose_representatives.py

7. Name representatives with GTDB:
   GutReferenceSet/Build_Species_Set/naming.py

8. Compare to Segata and build tree structure:
   GutReferenceSet/Utils/phyphlan.py


