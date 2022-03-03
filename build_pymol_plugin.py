import os
import tarfile

# compress the contents of a directory given by 'path' into a file 'tar_name'. and store it in the current directory
def tardir(path, tar_name):
    with tarfile.open(tar_name, "w:gz") as tar_handle:
        for root, dirs, files in os.walk(path):
            for file in files:
                tar_handle.add(os.path.join(root, file))

software_dir = "/home/ryan/Software"
ffeatools_dir = os.environ['FFEA_SRC'] + "/ffeatools"
tar_dir = software_dir + "/FFEAplugin"
pymol_build_dir = os.environ['FFEA_BUILD'] + "/share/ffea/plugins/pymol"
pymol_files =  ["FFEA_analysis/pymol_plugin/Blob.py",
    "FFEA_analysis/pymol_plugin/__init__.py",
    "FFEA_analysis/pymol_plugin/mtTkinter.py",
    "FFEA_analysis/pymol_plugin/system_names_pokemon.txt",
    "FFEA_analysis/pymol_plugin/system_names_dbzcharacters.txt",
    "FFEA_analysis/pymol_plugin/system_names_greekletters.txt",
    "modules/FFEA_trajectory.py", 
    "modules/FFEA_turbotrajectory.py",
    "modules/FFEA_measurement.py", 
    "modules/FFEA_frame.py",
    "modules/FFEA_universe.py",
    "modules/FFEA_binding_sites.py", 
    "modules/FFEA_script.py", 
    "modules/FFEA_pin.py", 
    "modules/FFEA_springs.py", 
    "modules/FFEA_material.py", 
    "modules/FFEA_topology.py", 
    "modules/FFEA_surface.py",
    "modules/FFEA_node.py", 
    "modules/FFEA_stokes.py", 
    "modules/FFEA_vdw.py", 
    "modules/FFEA_lj.py", 
    "modules/FFEA_pdb.py", 
    "modules/FFEA_beads.py", 
    "modules/FFEA_ctforces.py", 
    "modules/FFEA_exceptions.py", 
    "modules/FFEA_kinetic_map.py",
    "modules/FFEA_kinetic_rates.py",
    "modules/FFEA_kinetic_states.py",
    "modules/FFEA_rod.py",
    "modules/FFEA_skeleton.py"]

try:
    os.mkdir(tar_dir)
except OSError:
    pass
for f in pymol_files:
    os.system("cp -v " + ffeatools_dir + "/" + f + " " + tar_dir + "/")
    
tardir("." + tar_dir , tar_dir + ".tar.gz")
os.system("mv " + tar_dir + ".tar.gz " + pymol_build_dir + "/")
os.system("rm -r " + tar_dir)