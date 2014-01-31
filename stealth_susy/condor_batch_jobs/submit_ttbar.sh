#!/bin/bash

out="--output_dir /eos/uscms/store/user/btcarlso/trees/Jan30/"

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/ttJetsSemiLept_v1:/pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/ttJetsSemiLept_ext1_v1:/pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/ttJetsSemiLept_ext2_v1 $out --njobs 50 --o ttJetsSemiLept_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/ttJetsFullLept_v1  $out --njobs 20 --o ttJetsFullLept_trigger