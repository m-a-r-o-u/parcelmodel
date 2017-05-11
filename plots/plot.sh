#!/bin/bash

# note: if this is set to -gt -1 the /etc/hosts part is not recognized ( may be a bug )
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -t|--type)
    PLOT_TYPE=$PLOT_TYPE" $2"
    shift # past argument
    ;;
    -o|--output)
    FILE_NAME=$FILE_NAME" $2"
    shift # past argument
    ;;
    -n|--netcdfs)
    NETCDFS=${NETCDFS}" $2"
    shift # past argument
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 $1
    exit
fi

echo plot_type: "${PLOT_TYPE}"
echo outout_filename: "${FILE_NAME}"
echo netcdf_files: "${NETCDFS}"

SWD=$HOME/Daten/Microphysics/Boxmodel/box_tobi/plots/
FILE_NAME=($FILE_NAME)
COUNTER=0

for P in $PLOT_TYPE;
do
    PLT=${SWD}${P}.py
    python $PLT $NETCDFS ${FILE_NAME[$COUNTER]}
    #echo $PLT $NETCDFS ${FILE_NAME[$COUNTER]}
    echo finished: ${P} plot
    let COUNTER=$COUNTER+1
done
