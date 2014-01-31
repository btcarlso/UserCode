
#!/bin/bash

out="--output_dir /eos/uscms/store/user/btcarlso/trees/Jan30/"

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/RPV_UDD212_300 $out --njobs 1 --o UDD_300_trigger
./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/RPV_LQD221_M200 $out --njobs 1 --o RPV_200_trigger

./submitCondor_AnalyzeSusyNtuple.py --inputFolders /eos/uscms/store/user/btcarlso/RPV/RPV300/ $out --njobs 1 --o RPV_300_trigger
./submitCondor_AnalyzeSusyNtuple.py --inputFolders /eos/uscms/store/user/btcarlso/RPV/RPV400/ $out --njobs 1 --o RPV_400_trigger
./submitCondor_AnalyzeSusyNtuple.py --inputFolders /eos/uscms/store/user/btcarlso/RPV/RPV500_rerun/ $out --njobs 1 --o RPV_500_trigger
./submitCondor_AnalyzeSusyNtuple.py --inputFolders /eos/uscms/store/user/btcarlso/RPV/RPV600/ $out --njobs 1 --o RPV_600_trigger
./submitCondor_AnalyzeSusyNtuple.py --inputFolders /eos/uscms/store/user/btcarlso/RPV/RPV700/ $out --njobs 1 --o RPV_700_trigger
./submitCondor_AnalyzeSusyNtuple.py --inputFolders /eos/uscms/store/user/btcarlso/RPV/RPV800/ $out --njobs 1 --o RPV_800_trigger
./submitCondor_AnalyzeSusyNtuple.py --inputFolders /eos/uscms/store/user/btcarlso/RPV/RPV900/ $out --njobs 1 --o RPV_900_trigger
./submitCondor_AnalyzeSusyNtuple.py --inputFolders /eos/uscms/store/user/btcarlso/RPV/RPV1000/ $out --njobs 1 --o RPV_1000_trigger
