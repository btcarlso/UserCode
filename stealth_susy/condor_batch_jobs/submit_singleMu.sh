#!/bin/bash

out="--output_dir /eos/uscms/store/user/btcarlso/trees/Jan30/"
#echo ./submitCondor_AnalyzeSusyNtuple.py --prefix root://xrootd.rcac.purdue.edu/  --filelist_input filelist_singleMu.txt --inputFolders /store/user/yiiyama/ $out --njobs 100 --o singleMu --json Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt

./submitCondor_AnalyzeSusyNtuple.py --prefix root://xrootd.rcac.purdue.edu/  --filelist_input filelist_singleMu.txt --inputFolders /store/user/yiiyama/ $out --njobs 100 --o singleMu --json Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt 
