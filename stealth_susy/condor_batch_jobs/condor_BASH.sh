#!/bin/bash

export X509_USER_PROXY x509up_u46545

echo "Grid Init: " 
ls -l -h x509up_u46545

JOB_NUMBER=$1
WORK_DIR=`pwd`

tar -xzf fileLists.tgz

export VO_CMS_SW_DIR=/uscmst1/prod/sw/cms
source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc462
scramv1 project CMSSW CMSSW_5_3_8_patch3

cd CMSSW_5_3_8_patch3/src/
mv $WORK_DIR/condor_src.tgz .
tar -xzf condor_src.tgz
cd SusyAnalysis/SusyNtuplizer/macro
eval `scramv1 runtime -sh`
make

mv $WORK_DIR/ANALYZER .
mv $WORK_DIR/filelist_$JOB_NUMBER .
mv $WORK_DIR/filelist_outputdir.txt .
mv $WORK_DIR/JSON .
mv $WORK_DIR/SusyEvent* .


echo "Chaining Files from: filelist_$JOB_NUMBER" 

while read file
do
  	sed -i '19i chain.Add("'$file'");' ANALYZER
	echo "Chaining File: $file"
done < filelist_$JOB_NUMBER
ls filelist_outputdir.txt
output_dir=$(sed -n 1p filelist_outputdir.txt)
echo $output_dir

root -b -q -l ANALYZER

#touch result.root
echo "Moving files: root file, analyzer, and filelist"
cp *.root $output_dir/hist_analysis_$JOB_NUMBER.root
mv ANALYZER $WORK_DIR/ANALYZER_$JOB_NUMBER
mv filelist_$JOB_NUMBER $WORK_DIR/filelist_used_$JOB_NUMBER
cd $WORK_DIR
rm -rf CMSSW_5_3_8_patch3/
rm fileLists.tgz
rm *.root
#rm filelist_*
