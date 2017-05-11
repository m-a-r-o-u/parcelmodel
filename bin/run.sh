#!/bin/bash

FOLDER_NAME=change_updraft_all
W=(0.1 0.2 0.4 0.6 0.8 1.0 1.2 1.4)
E=(0.17 0.0833 0.00)
#E=(0.0833)

#FOLDER_NAME=change_radiation
#E=(1.0 0.5 0.333 0.17 0.0833 0.00)
#W=(1.0)

for E_I in ${E[*]};
do
for W_I in ${W[*]}; 
do
INPUT_FILE_NAME=config_W${W_I}_E${E_I}.yaml
cat > ${INPUT_FILE_NAME} <<EOF
output: 
    logger:
        - type: NetCDFLogger
          file_name: w_${W_I}E_${E_I}
globals:
    file_path: ./output/${FOLDER_NAME}/
    executer: numpy
initial_conditions:
    groups: 50000
    dt: 0.2
    output_step: 1
    t_max: 1400
    T: 281.7
    p: 89880.
    qv: 0.0077
    l: 50.
    w: ${W_I}
    z0: 1000.
    particle_distribution:
        type: double_log_normal
        ratio: 0.6
        mu1: 2.e-8
        sigma1: 1.4 
        mu2: 7.5e-8
        sigma2: 1.6
        total: 100.e+6
    radiation_schema:
        type: stefan_boltzmann
        factor: ${E_I}
        #type: no_radiation
    turbulence_schema:
        type: markov_schema
        epsilon: 50.e-4
        #type: no_turbulence
    atmosphere_schema:
        #type: afglus
        #atmosphere_file: ./input/afglus.dat
        type: linear_and_hydrostatic
EOF

python box.py ${INPUT_FILE_NAME} &
sleep 2
rm ${INPUT_FILE_NAME}

done
done
