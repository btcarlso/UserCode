#!/bin/bash

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/QCD4Jets_Pt100-180  --njobs 2 --o QCD4Jets_Pt100-180

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/QCD4Jets_Pt180-250  --njobs 2 --o QCD4Jets_Pt180-250

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/QCD4Jets_Pt250-400  --njobs 2 --o QCD4Jets_Pt250-400

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/QCD4Jets_Pt400-5600  --njobs 2 --o QCD4Jets_Pt400-5600

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/QCD6Jets_Pt100-180  --njobs 2 --o QCD6Jets_Pt100-180

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/QCD6Jets_Pt180-250  --njobs 2 --o QCD6Jets_Pt180-250
er

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/QCD6Jets_Pt250-400  --njobs 2 --o QCD6Jets_Pt250-400

./submitCondor_AnalyzeSusyNtuple.py --prefix dcap:// --inputFolders /pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/QCD6Jets_Pt400-5600  --njobs 2 --o QCD6Jets_Pt400-5600
