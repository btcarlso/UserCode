
#!/bin/bash
out="--output_dir /eos/uscms/store/user/btcarlso/trees/Jan30/"

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/WWJetsTo2L2Nu $out --njobs 5 --o WWJetsTo2L2Nu_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/WZJetsTo2L2Q/ $out --njobs 5 --o WZJetsTo2L2Q_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/WZJetsTo2QLNu $out --njobs 5 --o WZJetsTo2QLNu_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/WZJetsTo3LNu $out --njobs 5 --o WZJetsTo3LNu_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/ZZJetsTo2L2Q $out --njobs 5 --o ZZJetsTo2L2Q_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/ZZJetsTo2L2Nu $out --njobs 2 --o ZZJetsTo2L2Nu_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/ZZJets4L $out --njobs 10 --o ZZJetsTo4L_trigger
