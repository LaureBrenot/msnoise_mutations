# Setup Python and MSNoise
1) In a prompt run "pyhon setup_for_avosmart.py" with the data folder, compile_msnoise_C2.py, default.py, msnoise.py, and msnoise_mutations.py all in the same folder.\
2) Modify data_location, data_output, start_date and end_date in compile_msnoise_C2.py\
3) Run compile_msnoise_C2.py

setup_for_avosmart.py will execute the steps below until "Create new_folder testx with compile_msnoise_C2.py inside" (included).\
This setup uses SQLite. It means that for each run of setup_for_avosmart.py, a "test" folder will be created associated with one database. The code compile_msnoise_B2.py inside the "test" folder will be specific to the associated database.\
The output will be in the testx folder. The final result is in file "tt.csv". The intermediate results are in folders CC, STACKS, MWCS, DTT, and WCT ... forthcoming).

# Or do it manually
1) create a new environment python 3.7.16
2) activate the new environment 
3) conda install -c conda-forge flask-admin flask-wtf markdown folium pymysql logbook pywct
4) conda install -c conda-forge obspy
5) pip install msnoise \
It will install msnoise 1.6.3 version
6) pip install tables

add in .../site-package/msnoise/msnoise_mutations.py with get_ref and get_extension\
modify in .../site-package/msnoise/scripts/msnoise.py for msnoise_save.py\
add in .../site-package/msnoise/scripts/msnoise.py (the new one)\
modify in .../site-package/msnoise/scripts/default.py for default_save.py
add in .../site-package/msnoise/scripts/default.py (the new one)
Create new_folder testx with compile_msnoise_C2.py inside\

Modify data_location, data_output, start_date and end_date\
Run the code
