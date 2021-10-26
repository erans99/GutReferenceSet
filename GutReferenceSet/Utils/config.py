import os
from GutReferenceSet.Build_Species_Set.config import build_dir

checkm = "XXXX"  # path of checkm executable
prodigal = "XXXX"  # path of prodigal executable
hmmer = "XXXX"  # path of hmmer executable
pplacer = "XXXX"  # path of pplacer executable

mash_exe = "XXXX"  # mash executable

gtdb_exe = "XXXX"  # gtdbtk executable
gtdb_info = "XXXX"  # path of gtdbtk databases

phy_path = os.path.join(build_dir, "PhyPhlan")
phy_names = os.path.join(phy_path, "df_out_PhyPhlan.csv")
python_exe_dir = "XXXX"  # path of python3 executable (bin directory)
phylophlan_path = "XXXX"  # path of PhyloPhlan databases (so as not to download again and again)
Segata_desc = "XXXX"  # path to "SupplementaryTable8-SGBsDescription.csv" of Pasolli et al.'s paper

