# All the steps in [Setup_hpc_for_MSNoise](https://github.com/LaureBrenot/msnoise_mutations/blob/main/Parallelized_on_CECI/Setup_hpc_for_MSNoise.md)

# SUMMARY ONCE IT WORKS
./run_mariadb.sh\
mysql
```
create database databasename ;
quit 
```
cd Msnoise_para/volcanoname\
source activate msnoise_env\
msnoise db init (10.47.1.1 :5050)\
python 1msnoise_data_folder.py\
python 3msnoise_config.py\
msnoise populate			\
sbatch 2msnoise_scan_archive.sh\
msnoise db execute "select sta, count(*) from data_availability group by sta"\
msnoise db execute "update stations set used_location_codes='--'"\
msnoise config set enddate=2017-12-31\
msnoise db execute "select * from stations"\
msnoise info

(python 3msnoise_config.py)\
msnoise new_jobs –-init\
sbatch 4msnoise_cccompute.sh \
msnoise reset CC \
sbatch 5… 3 steps: ref: I, reset: T, done : D\
sbatch	6 : long\
./stop_mariadb.sh
