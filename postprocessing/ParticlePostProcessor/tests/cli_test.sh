#!/bin/bash
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
# ------------------------------------------------------------------
# Basic help test
ppp -h
check_process_success
# Should fail (wrong command)
ppp -Q
check_process_success
# Process command test
ppp process -h
check_process_success
ppp process -s ./data/particles.nc --topo ./data/topo.nc -o ./output.nc
check_process_success
# Should fail (file exists)
ppp process -s ./data/particles.nc --topo ./data/topo.nc -o ./output.nc
check_process_success
ppp process -O -s ./data/particles.nc --topo ./data/topo.nc -o ./output.nc
check_process_success