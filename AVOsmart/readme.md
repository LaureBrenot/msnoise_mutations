#SETUP PYTHON AND MSNOISE
1) create new environment python 3.7.16
2) activate the new environment 
3) conda install -c conda-forge flask-admin flask-wtf markdown folium pymysql logbook
4) conda install -c conda-forge obspy
5) pip install msnoise \
It will install msnoise 1.6.3 version
6) pip install tables

add in .../site-package/msnoise/msnoise_mutations.py with get_ref and get_extension\
modify in .../site-package/msnoise/scripts/msnoise.py for msnoise_save.py\
add in .../site-package/msnoise/scripts/msnoise.py (the new one)\
Create new_folder testx with compile_msnoise_B2.py inside\

Modify data_location, data_output, start_date and end_date\
Run the code

** OR **
run setup_for_avosmart.py with the data folder, compile_msnoise_B2.py, msnoise.py, and msnoise_mutations.py all in the same folder. setup_for_avosmart.py will execute the above steps until "Create new_folder testx with compile_msnoise_B2.py inside" (included).
