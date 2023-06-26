# SET UP
cd mariadb\
mkdir data\
cd ..\
mkdir tmp\
vi .my.cnf *configuration settings for MariaDB*\
```
[mysqld]
port=5050
bind-address=0.0.0.0
max_connections=1000
basedir=/home/ulb/gtime/lbrenot/mariadb
tmpdir=/home/ulb/gtime/lbrenot/tmp
datadir=/home/ulb/gtime/lbrenot/mariadb/data
socket=/home/ulb/gtime/lbrenot/mariadb/socket.sock

[mysql]
port=5050
socket=/home/ulb/gtime/lbrenot/mariadb/socket.sock
```
cd mariadb/ \
touch socket.sock 

cd bin/ \
mysql_install_db --defaults-file=~/.my.cnf --auth-root-authentication-method=normal    *install and initialize the MariaDB database + create all the necessary tables in the "data" directory*\
/home/ulb/gtime/lbrenot/mariadb/bin/mysqld_safe --datadir='/home/ulb/gtime/lbrenot/mariadb/data'  *Start the MariaDB server in safe mode, specifying the data directory*\

cd\
vi run_mariadb.sh
```
ifconfig ib0 | sed -En 's/127.0.0.1//;s/.*inet (addr:)?(([0-9]*\.){3}[0-9]*).*/\2/p' > ~/MARIADB.ini
cd mariadb/
./bin/mysqld_safe --defaults-file=~/.my.cnf --skip-grant-tables &     *Start the MariaDB server in safe mode, skipping the grant tables (which control user access), and run it in the background*
```

vi stop_mariadb.sh
```
./mariadb/bin/mysqladmin -h $( cat MARIADB.ini ) --port 5050 -u root shutdown   *connect to the MariaDB server and initiate a shutdown*
```
chmod +x *.sh     *Change the permissions of all files with the ".sh" extension to make them executable*\
mysql
```
(show databases ;)
(DROP DATABASE database_name;)
create database databasename
use databasename
(show tables ;)
quit 
```

cd Msnoise_para\
source activate msnoise_env *Activate the "msnoise_env" virtual environment*\
msnoise db init *in the good folder*\
```
10.47.1.1 :5050
databasename 
laure
Msnoise / noise
[] *easier if empty*
```
msnoise test

# SUBMISSION
vi msnoise_data_folder.py 
```
# -*- coding: utf-8 -*-
from msnoise import api
from msnoise import s000installer

db = api.connect() #connect to database

volcano = "Gareloi"
data_location = "/CECI/trsf/ulb/gtime/lbrenot/Alaska_volcanoes/"+volcano+"_long_term_csn_ASY/data_"+volcano
api.update_config(db, 'data_folder', data_location) 
```
python msnoise_data_folder.py\
msnoise config get data_folder\
msnoise populate\
msnoise info\
vi msnoise_scan_archive.sh
```
#!/bin/sh
#SBATCH --job-name=msn_scan
#SBATCH --output=res_msn_scan.txt
#
#SBATCH --ntasks=1
#SBATCH --time=30:00
#SBATCH --mem-per-cpu=500
#
#SBATCH --cpus-per-task=8
source /home/ulb/gtime/lbrenot/.bashrc
cd ./Msnoise_para/
conda activate msnoise_env
srun msnoise -t 8 scan_archive --init
```
sbatch msnoise_scan_archive.sh			*=msnoise scan_archive –init*

tail -f res_msn_scan.txt\
msnoise db execute "select count(*) from data_availability" 	ou "… a_data_availability"\
msnoise db execute "select sta, count(*) from data_availability group by sta"


msnoise db update_loc_chan\
vi msnoise_config.py
*Too long, see at the end*\
python msnoise_config.py

msnoise db execute "update a_stations set used_location_codes='--'"\
msnoise db execute "select * from a_stations"

msnoise info

msnoise new_jobs –init

msnoise db execute "insert into a_filters (ref, low, mwcs_low, high, mwcs_high, rms_threshold, mwcs_wlen, mwcs_step, used) values (1, 0.1, 0.1, 1.0, 1.0, 0.0, 12.0, 4.0, 1)"

vi msnoise_cccompute.sh 
```
#!/bin/sh
#SBATCH --job-name=msn_computecc
#SBATCH --output=res_msn_computecc.txt
#
#SBATCH --ntasks=1
#SBATCH --time=30:00
#SBATCH --mem-per-cpu=2000
#
#SBATCH --array=1-250
#
#SBATCH --mail-user=laure.brenot@ulb.be
#SBATCH --mail-type=ALL
#
source /home/ulb/gtime/lbrenot/.bashrc
conda activate msnoise_env
srun msnoise compute_cc
```
sbatch msnoise_cccompute.sh 

squeue --me\
watch -n 10 "msnoise info -j"\
(need to restart ? « msnoise reset CC » before)\

msnoise plot ccftime AV.GAEA AV.GANE\
msnoise plot interferogram AV.GAEA AV.GANE\
…\
./stop_mariadb.sh





```
# -*- coding: utf-8 -*-

from msnoise import api
from msnoise import s000installer

scan 		= False
newjobs 	= False
compute_cc 	= True
stacks 		= True
compute_dtt 	= True

db = api.connect() #connect to database

#api.reset_jobs(db, 'CC', alljobs=True)
#api.reset_jobs(db, 'MWCS', alljobs=True)
#api.reset_jobs(db, 'DTT', alljobs=True)

################################## SCAN ARCHIVE
start 	= '2019-01-01'
end 	= '2019-12-31'
if scan:
    api.update_config(db, 'startdate', start) #set startdate
    api.update_config(db, 'enddate', end) #set enddate
##ready to scan archive

################################# NEWJOBS
if newjobs:
    api.update_station(db, used_locations_codes='--')
    #api.update_station(db,'NZ','WIZ',177.1894,-37.5265,40,coordinates='DEG') # enter location of each station
    #api.update_station(db,'NZ','WSRZ',177.178,-37.518,290,coordinates='DEG')

    api.update_config(db, 'components_to_compute', 'ZZ') #for station-pair correlation functions
    #api.update_config(db, 'components_to_compute_single_station', 'EN,EZ,NZ')
##ready new_jobs init

################################ COMPUTE CC
if compute_cc:
    api.update_config(db, 'cc_sampling_rate', '20') #set frequency to resample to, default=20
    api.update_config(db, 'preprocess_lowpass', '8') #preprocess lowpass filter frequency, default=8
    api.update_config(db, 'preprocess_highpass', '0.01') #preprocess highpass filter frequency, default=0.01
    api.update_config(db, 'maxlag', '120') #maximum lag of cross-correlation functions, default=120
    api.update_config(db, 'corr_duration', '1800') #length of time slices (seconds) to cross-correlate, default=1800
    api.update_config(db, 'windsorizing', '3') #temporal normalization method (3 = clip 3 times RMS)
    api.update_config(db, 'whitening_type', 'B') #type of whitening (B = whiten to amplitude of 1.0)
    api.update_config(db, 'clip_after_whiten', 'Y') #Boolean to perform temporal norm after spectral whitening
    api.update_config(db, 'stack_method', 'linear') # stack method, linear or pws (phase-weighted)
    #api.update_config(db, 'remove_response', 'Y') #remove instrument response

    api.update_filter(db, 1, 0.1, 0.15, 1.0, 0.95, 0, 12, 2, True)
    api.update_filter(db, 2, 1.0, 1.05, 2.0, 1.95, 0, 12, 2, True)

#ready to compute_cc

#################################### STACK
#to record dvv : compare current and reference stacks of ccfs. Parameters to decide of the size stacks
if stacks:
    api.update_config(db, 'mov_stack', '5,10') #day length windows
    api.update_config(db, 'ref_begin', start)#define the reference stack
    api.update_config(db, 'ref_end', end)

#ready stack and compute_mwcs

################################### COMPUTE DTT
if compute_dtt:
    api.update_config(db, 'dtt_lag', 'static') #set lag time window to be the same for all CCFs
    api.update_config(db, 'dtt_minlag', '10') #set minimum lag time
    api.update_config(db, 'dtt_width', '30') #set width of window to compute dv/v , in sec
    api.update_config(db, 'dtt_sides', 'both') #choose which side(s) of CCF to compute dv/v
    api.update_config(db, 'dtt_maxerr', '0.2') #set maximum error
    api.update_config(db, 'dtt_maxdt', '0.2') #set maximum delay time
    api.update_config(db, 'dtt_mincoh', '0.6') #set minimum coherence
```
