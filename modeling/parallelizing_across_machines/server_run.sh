module load mpi/openmpi-x86_64
cd /data/desmosome/newrun_main || exit

if [ ! -f ./logs/server_output_list ]; then
    echo server "$HOSTNAME" >> ./logs/server_output_list
fi

if [ ! -d ./all_outputs/output$1 ]; then
    mkdir ./all_outputs/output$1
fi

~/imp-custom/build/setup_environment.sh mpirun -np 8 python -u representation_sampling.py all_outputs/output$1 1> logs/$1_stream1.log 2> logs/$1_stream2.log
echo $1 >> ./logs/server_output_list

cd ~ || exit

