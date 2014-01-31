
#!/bin/bash

out="--output_dir /eos/uscms/store/user/btcarlso/trees/Jan30/"
#echo "./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/wJetsToLNu_v1 $out --njobs 20 --o wJetsToLNu_trigger"

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/wJetsToLNu_v1 $out --njobs 20 --o wJetsToLNu_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/w1JetsToLNu_v1 $out --njobs 20 --o w1JetsToLNu_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/w2JetsToLNu_v1 $out --njobs 20 --o w2JetsToLNu_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/w3JetsToLNu_v1 $out --njobs 20 --o w3JetsToLNu_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/w4JetsToLNu_v1  $out --njobs 20 --o w4JetsToLNu_trigger
