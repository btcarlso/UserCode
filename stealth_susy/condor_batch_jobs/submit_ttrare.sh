#!/bin/bash

out="--output_dir /eos/uscms/store/user/btcarlso/trees/Jan30/"

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/ttJetsSemiLept_v1:/pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/ttJetsSemiLept_ext1_v1:/pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/TTTT $out --njobs 1 --o TTTT

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/TTWJets $out --njobs 1 --o TTWJets

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/TTWJets $out --njobs 1 --o TTWJets

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/TTWWJets $out --njobs 1 --o TTWWJets

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/TTZJets $out --njobs 1 --o TTZJets