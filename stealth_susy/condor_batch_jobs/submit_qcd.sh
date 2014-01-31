#!/bin/bash

out="--output_dir /eos/uscms/store/user/btcarlso/trees/Jan30/"

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/QCD_MuEnriched_pt30-50  $out --njobs 20 --o QCD_30-50_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/QCD_MuEnriched_pt50-80  $out --njobs 20 --o QCD_50-80_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/QCD_MuEnriched_pt80-120  $out --njobs 20 --o QCD_80-120_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/QCD_MuEnriched_pt120-170  $out --njobs 20 --o QCD_120-170_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/QCD_MuEnriched_pt170-300  $out --njobs 20 --o QCD_170-300_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/QCD_MuEnriched_pt300-470  $out --njobs 20 --o QCD_300-470_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/QCD_MuEnriched_pt470-600  $out --njobs 20 --o QCD_470-600_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/QCD_MuEnriched_pt600-800  $out --njobs 20 --o QCD_600-800_trigger

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/QCD_MuEnriched_pt800-1000  $out --njobs 20 --o QCD_800-1000_trigger
