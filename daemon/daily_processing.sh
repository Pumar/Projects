#!/bin/bash
source /opt/root/bin/thisroot.sh

export MODEXP_BASE_DIR=/home/process_modulation_joranversion

export MODEXP_RAW_DATA_DIR=/home/process_modulation_joranversion/rawdata

export MODEXP_PROCESSED_DATA_DIR=/media/sf_Modulation/Modulation_data/Processed

export MODEXP_ANALYSIS_DIR=$MODEXP_BASE_DIR/analysis-master

export MODEXP_DAQANA_DIR=$MODEXP_BASE_DIR/daqana-joran

export MODEXP_TEMP_SCRIPTS_DIR=$MODEXP_BASE_DIR/temp-scripts

export MODEXP_CALIBRATION_DATA_DIR=$MODEXP_PROCESSED_DATA_DIR/calibration

export MODEXP_ANALYSIS_DATA_DIR=$MODEXP_PROCESSED_DATA_DIR/analysis

export ROOTSYS=/opt/root

raw_data_dir=/media/sf_Modulation/Modulation_data/Commissioning/Data/Collaboration/Sources
proc_data_dir=/media/sf_Modulation/Modulation_data/Processed
plot_dir=/media/sf_Modulation/Modulation_data/plots/pumar.github.io/daily_plot
cd $raw_data_dir
yest_full=$(date -d "1 day ago" "+%A, %d %B %Y")
yest=$(date -d "1 day ago" "+%Y%m%d")
recent_file=$(find . -type d -name "*${yest}*")

cd /home/process_modulation_joranversion/daqana-joran
if [ -z "$recent_file" ]
then
    echo $yest_full >> /home/process_modulation_joranversion/daqana-joran/daily_processing_logs/no_raw_data_dates.txt
else
    ./daqprocessor_single_calibrate.py $raw_data_dir/$recent_file/ -l -s -f -p0 > /home/process_modulation_joranversion/daqana-joran/daily_processing_logs/file_output.txt
    cd /opt/processing_pipe/plotter
    python slowdata.py $proc_data_dir/$recent_file/
    python fastdata.py $proc_data_dir/$recent_file/ > /home/process_modulation_joranversion/daqana-joran/daily_processing_logs/plotter_output.txt
    cd $plot_dir
    git add --all
    git commit -m "daily plot update"
    git push
fi


