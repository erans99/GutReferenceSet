import os
from GutReferenceSet.Build_Species_Set.config import build_dir

checkm = "XXXX"  # path of checkm executable
prodigal = "XXXX"  # directory of prodigal executable
hmmer = "XXXX"  # directory of hmmer executable
pplacer = "XXXX"  # directory of pplacer executable

mash_exe = "XXXX"  # path of mash executable

gtdb_exe = "XXXX"  # path of gtdbtk executable
gtdb_info = "XXXX"  # directory of gtdbtk databases

phy_path = os.path.join(build_dir, "PhyPhlan")
phy_names = os.path.join(phy_path, "df_out_PhyPhlan.csv")
python_exe_dir = "XXXX"  # directory of python3 executable (bin directory)
phylophlan_path = "XXXX"  # directory of PhyloPhlan databases (so as not to download again and again)
Segata_desc = "XXXX"  # directory to "SupplementaryTable8-SGBsDescription.csv" of Pasolli et al.'s paper

