#!/bin/bash
#-----------
# Purpose: Test the sort command

function print_separator {
    echo "------------------------------------------------------------"
}

Red='\033[0;31m'
Green='\033[0;32m'
NC='\033[0m'
function check_process_success {
    err=$?
    if [ $err -eq 0 ]; then
        echo -e "${Green}Process success${NC}"
    else
        echo -e "${Red}Process failed:${NC} $err"
    fi
    print_separator
}

ppp -h
check_process_success
ppp process -h
check_process_success
ppp process -s ./data/particles.nc --topo ./data/topo.nc -o ./output.nc -O
check_process_success
ppp process -s ./data/particles.nc --topo ./data/topo.nc -o ./output.nc --sort id -O
check_process_success
ppp process -s ./data/particles.nc --topo ./data/topo.nc -o ./output.nc --sort all id -O
check_process_success
ppp process -s ./data/particles.nc --topo ./data/topo.nc -o ./output.nc --sort all state -O
check_process_success
ppp process -s ./data/particles.nc --topo ./data/topo.nc -o ./output.nc --sort id state -O
check_process_success

