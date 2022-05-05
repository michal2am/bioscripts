#!/bin/bash
# conda install -c conda-forge pyngrok
# cron -e
# */1 * * * * /home/mm/ssh_autorun.sh

pkill ngrok
sleep 10s
nmcli radio wifi off
sleep 10s
nmcli radio wifi on
sleep 20s
/home/mm/anaconda3/envs/37_omen/bin/python /home/mm/repos/bioscripts/stockholm/mdanalysis/ngrok_auto.py
