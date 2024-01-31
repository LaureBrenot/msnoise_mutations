# -*- coding: utf-8 -*-
"""
SET UP MSNOISE WATCH
    create new environment python 3.7.16
    activate the new environment
    conda install -c conda-forge flask-admin flask-wtf markdown folium pymysql logbook pywct
    conda install -c conda-forge obspy
    pip install msnoise
    It will install msnoise 1.6.3 version
    pip install tables

add in .../site-package/msnoise/msnoise_mutations.py with get_ref and get_extension
modify in .../site-package/msnoise/scripts/msnoise.py for msnoise_save.py
add in .../site-package/msnoise/scripts/msnoise.py (the new one)
modify in .../site-package/msnoise/scripts/default.py for default_save.py
add in .../site-package/msnoise/scripts/default.py (the new one)

Create new_folder testx with compile_msnoise_C2.py inside
"""

import shutil
import os

######## WHAT DO YOU WANT? #########
## create new environment python 3.7.16 and install packages
env_name = "msnoise_watch_env"
python_version = "3.7.16"

############ CODE ###############

print(f"*** Creating a new Python environment: {env_name} with Python version {python_version} ***\n")
os.system(f"conda create -n {env_name} python={python_version}")
#os.system(f"python -m venv {env_name} python={python_version}")

print(f"*** Install missing package in {env_name} ***\n")
os.system(f"conda activate {env_name} && conda install -c conda-forge flask-admin flask-wtf markdown folium pymysql logbook")
print(f"*** Install obspy in {env_name} ***\n")
os.system(f"conda activate {env_name} && conda install -c conda-forge obspy")
print(f"*** Install msnoise in {env_name} ***\n")
os.system(f"conda activate {env_name} && pip install msnoise")
print(f"*** Install tables in {env_name} ***\n")
os.system(f"conda activate {env_name} && pip install tables")
print(f"*** Install pywct in {env_name} ***\n")
os.system(f"conda activate {env_name} && pip install pywct")

## Move files to the appropriate location
# msnoise_muations.py
mut_file = "msnoise_mutations.py"
env_path = os.environ['CONDA_PREFIX']

msn_dir = os.path.join(env_path, "envs", env_name, "Lib\site-packages\msnoise")

if not os.path.exists(msn_dir):
    print("Error: Path does not exist -", msn_dir)
else:
    print(f"*** Copying {mut_file} to {msn_dir} ***\n")
    shutil.copy(os.path.join(os.getcwd(), mut_file), msn_dir)

# msnoise.py
old_file_path = os.path.join(env_path,"envs", env_name, "Lib\site-packages\msnoise\scripts\msnoise.py")
dir_path = os.path.dirname(old_file_path)
new_file_name = "msnoise_save.py"
new_file_path = os.path.join(dir_path, new_file_name)

print(f"*** Renaming {old_file_path} to {new_file_path} and coping the new msnoise.py to {dir_path} ***\n")
os.rename(old_file_path, new_file_path)
shutil.copy(os.path.join(os.getcwd(), "msnoise.py"), dir_path)

# default.py
os.rename(os.path.join(msn_dir, "default.py"),os.path.join(msn_dir, 'default_save.py'))
shutil.copy(os.path.join(os.getcwd(), "default.py"), msn_dir)

## Create new folder and copy compile_msnoise_B2.py inside
dirs = next(os.walk('.'))[1]

# find the highest numbered test folder
highest_num = 0
for d in dirs:
    if d.startswith('test'):
        num = int(d[4:])
        if num > highest_num:
            highest_num = num

# create the new folder with the next number
new_folder = f'test{highest_num + 1}'
os.mkdir(new_folder)

file_name = 'compile_msnoise_C2.py'
comp_path = os.path.join(os.getcwd(), file_name)
dst_path = os.path.join(os.getcwd(), new_folder, file_name)
print(f"*** Copying {file_name} to {dst_path} ***")
shutil.copy(comp_path, dst_path)

#import subprocess
#file_path = f"{new_folder}/compile_msnoise_B2.py"
#subprocess.run(["nano", file_path])
#print("*** File has been edited. ***")

print(f"\n \n *** Don't forget to modify data_location, data_output, start_date and end_date in {new_folder}/compile_msnoise_C2.py ***")
print("Run the new script with the following command:")
print(f"conda activate {env_name} && python {new_folder}/compile_msnoise_C2.py")


## Run the new script
#os.system(f"conda activate {env_name} && python {new_folder}/compile_msnoise_B2.py")

