#!/bin/bash

out="--output_dir /eos/uscms/store/user/btcarlso/trees/Jan30/"

#echo "./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/WWWJets $out --njobs 1 --o WWWJets"

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/WWWJets $out --njobs 1 --o WWWJets

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/WWZJets $out --njobs 1 --o WWZJets

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/ZZZJets $out --njobs 1 --o ZZZJets
