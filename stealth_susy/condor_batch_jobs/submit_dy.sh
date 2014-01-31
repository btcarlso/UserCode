#!/bin/bash
out="--output_dir /eos/uscms/store/user/btcarlso/trees/Jan30/"

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/dyJetsToLL_v1 $out --njobs 20 --o dyJetsToLL_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/dy1JetsToLL_v1  $out --njobs 20 --o dy1JetsToLL_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/dy2JetsToLL_v1  $out --njobs 20 --o dy2JetsToLL_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/dy3JetsToLL_v1  $out --njobs 20 --o dy3JetsToLL_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/dy4JetsToLL_v1  $out --njobs 20 --o dy4JetsToLL_trigger
