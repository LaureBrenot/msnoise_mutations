#SETUP PYTHON AND MSNOISE
1) create new environment python 3.7.16
2) activate the new environment 
3) conda install -c conda-forge flask-admin flask-wtf markdown folium pymysql logbook
4) conda install -c conda-forge obspy
5) pip install msnoise \
It will install msnoise 1.6.3 version
6) pip install tables

add in .../site-package/msnoise/msnoise_mutations.py\
modify in .../site-package/msnoise/scripts/msnoise.py for msnoise_save.py\
add in .../site-package/msnoise/scripts/msnoise.py (the new one)\

