# do all QTLs from list of config files
# send off separate snakemake job to HPC 

# this will allow me to run all analyses in parallel

script=snakejob_HPC
# use bash getopt to use short flags
config='config.yaml'
snakefile='Snakefile'
dryrun=""
mode="eQTL"

print_usage() {
  printf "Usage: ./run_all_QTLs.sh [-h] [-c] config list [-s] Snakefile\nOptions:\n\t-c path to config file; OR path to file listing paths to config files\n\t-s Snakefile\n\t-n dry run mode\n\t-m mode [eQTL/sQTL]\n"
}

while getopts 'c:s:m:nh' flag; do
  case "${flag}" in
    c) configList="${OPTARG}" ;;
    s) snakefile="${OPTARG}" ;;
    m) mode="${OPTARG}" ;;
    n) dryrun="-n" ;;
    h) print_usage
       exit 1 ;;
  esac
done

echo config list used is $configList
echo snakefile used is $snakefile 
echo mode is $mode

if [ ! -z "$dryrun" ];then 
    echo "dry run mode" 
fi 
 
if [ ! -e $configList ]; then 
    echo "error: config file does not exist"  
    exit 0 
fi 
 
if [ ! -e $snakefile ]; then 
    echo "error: snakefile does not exist" 
    exit 0 
fi 
 


# CONFIG can be either a file containing a list of configs or just a single config.yaml

#conda activate snakemake


if grep -q "yaml" <<< "$configList"; then
    echo using $configList
    if [ $mode == "eQTL" ]; then
        bsub -P acc_als-omics -W 24:00 -n 1 -q premium -o cluster/snakejob_HPC.stdout -e cluster/snakejob_HPC.stderr -L /bin/bash  "sh $script -s Snakefile -c $configList -m eQTL"
    fi
    if [ $mode == "sQTL" ]; then
        bsub -P acc_als-omics -W 24:00 -n 1 -q premium -o cluster/snakejob_HPC.stdout -e cluster/snakejob_HPC.stderr -L /bin/bash  "sh $script -s Snakefile -c $configList -m sQTL"
    fi
    exit 0 
elif [ -f $configList ]; then

    for config in $(cat $configList);do
        echo $config
            if [ $mode == "eQTL" ]; then
                bsub -P acc_als-omics -W 24:00 -n 1 -q premium -o cluster/snakejob_HPC.stdout -e cluster/snakejob_HPC.stderr -L /bin/bash  "sh $script -s Snakefile -c $config -m eQTL"
            fi
            if [ $mode == "sQTL" ]; then
                bsub -P acc_als-omics -W 24:00 -n 1 -q premium -o cluster/snakejob_HPC.stdout -e cluster/snakejob_HPC.stderr -L /bin/bash  "sh $script -s Snakefile -c $config -m sQTL"
            fi
    done

fi
