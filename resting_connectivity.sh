#!/bin/bash

##########################################
# Resting Connectivity Script v3
# 10/8/09 
# Donald McLaren (dm10@medicine.wisc.edu)
# Wisconsin Alzheimer's Disease Center Imaging Core (www.brainmap.wisc.edu)
#
# USAGE:
#       resting_connectivity -S study_directoy -s subject -V regionfile -R regionname --files rest
#
# This script will compute the seed based connectivity for a specified brain region with every other voxel.
# It has the option to regress out the effects of cardiac/respiratory data. If cardiac/respiratory data is not available, 
# then it regresses out the effects of the ventricles and deep white matter. It will also output the results of
# regressing out the global signal.
#
# Additional Options:
#    -- Partial correlation (--partial regionlist.txt), this will regress out the effects of other brain regions
#       (see van den Heuval, et al. Microstructural Organization of the Cingulum Tract and the Level of Default 
#       Mode Functional Connectivity. J. Neuroscience 2008; 28(43):10844-10851.)
#    -- Regress out a task (--task convolvedtaskregressorfile.txt), this option is NOT recommended!!!! If you feel
#       to use task data, against recommendations, use this option.USER WARNING: This has not been fully tested!!!!
#       Use at your own risk.
# For complete required and optional input listing, type: resting_connectivity.sh --help
#
# Before running, one must customize the default options below for the current system configuration and data structure.
#
##########################################

###################
# Default Options 
###################
lowpass=.08                                        # Defualt option
highpass=.009                                      # Defualt option
motion=`echo "rp*${files}*.txt"`                   # Default option from SPM
ROIDIR=/Data/data1/ROI/mclaren/                    # Directory of white matter and ventricle masks
subject_dir=/fmri                                  # Where smoothed and normalized data live; must have / in front; files must be sw* and not zipped in any fashion
regonly=0                                          # Default to do all the processing steps
corronly=0                                         # Default to do all the processing steps
cardiacfile=dummy.txt                              # Place holder for optional cardiac data
partialcorr=0;                                     # Default option, do not compute the partial correlations
task=0                                             # Default option, this is resting data -- there is no task effect to remove
RTB=/Data/data1/lab_scripts/RestingToolbox         # Location of the RestingToolbox directory
LPF=0                                              # Default option, do not process LPF with resting_correlation.sh
separate=0                                         # Default option, do not separately process each file with resting_correlation.sh
PATH=$PATH:${RTB}

#######################
# Function Setup
#######################
help_text () {
echo " "
echo "This script will compute the seed based connectivity for a specified brain region with every other voxel."
echo "It has the option to regress out the effects of cardiac/respiratory data. If cardiac/respiratory data is not available," 
echo "then it regresses out the effects of the ventricles and deep white matter. It will also output the results of"
echo "regressing out the global signal.Insert brief description"
echo " "
echo "Available Options:"
echo 
echo "      -S | --study    ) # Required! Sets Study Directory"
echo
echo "      -s | --subject  ) # Required! Sets Subject Number"
echo
echo "      -V | --VOI      ) # Required! Sets VOI file (directory necessary)"
echo
echo "      -R | --REGION   ) # Required! Adds REGION NAME to output files."
echo
echo "      --files         ) # Required! Sets the search term for files (e.g.rest)"
echo
echo "      --lowpass       ) # Optional. Changes the default low pass filter from .08"
echo
echo "      --highpass      ) # Optional. Changes the default high pass filter from .009"
echo
echo "      --motionfiles   ) # Optional. Changes the default string for motion files from rp*txt (e.g. motion*${files}*.txt)"
echo
echo "      --regonly       ) # Optional! Program will only do the regression computations. "
echo
echo "      --cardiac       ) # Optional! Use the cardiac file specified after the flag (e.g. cardiac*${files}*.txt)"
echo
echo "      --partial       ) # Optional! Compute the partial correlations. Use all VOI specified in the list. The list should be stored in a file."
echo
echo "      --separate      ) # Optional!!! In addition to processing the files jointly, compute maps for each run on its own within resting_correlation.sh"
echo
echo "      --LPF           ) # Optional!!! Process the LPF files with resting_correlation in addition to the BPF files. "
echo " "
echo "      --task          ) # Optional! NOT RECOMMENDED!!!!!! If you are inclined to use task data, then use to regress out the effect of task."
echo "                        After the flag, enter the name of the file that contains the convolved/FIR/basis function task regressors. USER WARNING: this"
echo "                        option has not been fully tested! Use at your own risk."
echo " "
echo "      -h | --help     ) # Get this Help message."
}

run_regress () {
  if [ $1 == LPF ] && [ `ls PSC_LPF${lowpass}_sw*.nii.gz | wc -l` -gt 0 ]; then
    FILES2=(`ls PSC_LPF${lowpass}_sw*.nii.gz`)
  elif [ $1 == BPF ] && [ `ls PSC_BPF${lowpass}${highpass}_sw*.nii.gz | wc -l` -gt 0 ]; then
    FILES2=(`ls PSC_BPF${lowpass}${highpass}_sw*.nii.gz`)
  else
    echo "Not a valid option OR .gz files do not exist, program will now exit"
    exit
  fi
  MOTION2=(`ls ../${motion}`) 
  if [ $cardiacfile != dummy.txt ]; then
     CARDIAC2=(`ls ../${cardiacfile}`)
  fi
  CNT=${#FILES[@]}
  i="0"
  while [ $i -lt $CNT ]
  do
    J=${FILES2[$i]%.nii.gz}
    3dROIstats -mask ${voi} ${FILES2[$i]} | awk 'NR != 1 {print $3}' > rest${J}_${region}.txt
    if [ $partialcorr == 1 ]; then
     partialdir=`basename ${partiallist} | awk -F. '{print $1}' `
     if ! [ -d ${partialdir} ]; then
          mkdir ${partialdir}
          pushd $partialdir
     else
          pushd $partialdir
     fi
     partialdata=`basename ${partiallist} | awk -F. '{print $1}' `
     partialdata=${partialdata}_rest${J}.txt
     if [ -e $partialdata ]; then
        echo "Partial Data already exists, rename your list"
        exit
     else
        touch $partialdata
     fi
     for j in `more ${partiallist} `
     do
         partialregion=`basename ${j} | awk -F. '{print $1}' `
         3dcalc -a ${j} -expr 'a' -prefix tmp.nii -datum short; mv tmp.nii ${j}
         3dROIstats -mask ${j} ../${FILES2[$i]} | awk 'NR != 1 {print $3}' > rest${J}_${partialregion}.txt
         paste $partialdata rest${J}_${partialregion}.txt > $partialdata
     done
     PARTIAL=${partialdir}/${partialdata}
     popd
    fi
    
    3dROIstats -mask ${subject_num}_WM.nii.gz ${FILES2[$i]} | awk 'NR != 1 {print $3}' > rest${J}_WM.txt
    3dROIstats -mask ${ROIDIR}/wm.nii ${FILES2[$i]} | awk 'NR != 1 {print $3}' > F_rest${J}_WM.txt
    3dROIstats -mask ${subject_num}_VENT.nii.gz ${FILES2[$i]} | awk 'NR != 1 {print $3}' > rest${J}_VENT.txt
    3dROIstats -mask ${ROIDIR}/ventricles.nii ${FILES2[$i]} | awk 'NR != 1 {print $3}' > F_rest${J}_VENT.txt
    3dROIstats -mask ${ROIDIR}/BRAIN.nii ${FILES2[$i]} | awk 'NR != 1 {print $3}' > rest${J}_BRAIN.txt

    VOI=rest${J}_${region}.txt
    WM=rest${J}_WM.txt
    VENT=rest${J}_VENT.txt
    BRAIN=rest${J}_BRAIN.txt    

    if [ $cardiacfile == dummy.txt ] && [ $partialcorr == 0 ]; then
      if [ $task == 1 ]; then
        matlab -nojvm -r "addpath ${RTB}; Connectivity_regress('${subject_num}','${region}','${VOI}','motion','${MOTION2[$i]}','whitematter','${WM}','ventricles','${VENT}','task','${taskfile}')"
        matlab -nojvm -r "addpath ${RTB}; Connectivity_regress('${subject_num}','${region}','${VOI}','motion','${MOTION2[$i]}','whitematter','${WM}','ventricles','${VENT}','brain','${BRAIN}','task','${taskfile}')"
      else
        matlab -nojvm -r "addpath ${RTB}; Connectivity_regress('${subject_num}','${region}','${VOI}','motion','${MOTION2[$i]}','whitematter','${WM}','ventricles','${VENT}')"
        matlab -nojvm -r "addpath ${RTB}; Connectivity_regress('${subject_num}','${region}','${VOI}','motion','${MOTION2[$i]}','whitematter','${WM}','ventricles','${VENT}','brain','${BRAIN}')"
      fi
    elif [ $cardiacfile == dummy.txt ]; then
      if [ $task == 1 ]; then
        matlab -nojvm -r "addpath ${RTB}; Connectivity_regress('${subject_num}','${region}','${VOI}','motion','${MOTION2[$i]}','whitematter','${WM}','ventricles','${VENT}','partials','${PARTIAL}','task','${taskfile}')"
        matlab -nojvm -r "addpath ${RTB}; Connectivity_regress('${subject_num}','${region}','${VOI}','motion','${MOTION2[$i]}','whitematter','${WM}','ventricles','${VENT}','partials','${PARTIAL}','brain','${BRAIN}','task','${taskfile}')"
      else
        matlab -nojvm -r "addpath ${RTB}; Connectivity_regress('${subject_num}','${region}','${VOI}','motion','${MOTION2[$i]}','whitematter','${WM}','ventricles','${VENT}','partials','${PARTIAL}')"
        matlab -nojvm -r "addpath ${RTB}; Connectivity_regress('${subject_num}','${region}','${VOI}','motion','${MOTION2[$i]}','whitematter','${WM}','ventricles','${VENT}','partials','${PARTIAL}','brain','${BRAIN}')"
      fi
    elif [ $partialcorr == 0 ]; then
      if [ $task == 1 ]; then
        matlab -nojvm -r "addpath ${RTB}; Connectivity_regress('${subject_num}','${region}','${VOI}','motion','${MOTION2[$i]}','cardiac','${CARDIAC2[$i]}','task','${taskfile}')"
        matlab -nojvm -r "addpath ${RTB}; Connectivity_regress('${subject_num}','${region}','${VOI}','motion','${MOTION2[$i]}','cardiac','${CARDIAC2[$i]}','brain','${BRAIN}','task','${taskfile}')"
      else
        matlab -nojvm -r "addpath ${RTB}; Connectivity_regress('${subject_num}','${region}','${VOI}','motion','${MOTION2[$i]}','cardiac','${CARDIAC2[$i]}')"
        matlab -nojvm -r "addpath ${RTB}; Connectivity_regress('${subject_num}','${region}','${VOI}','motion','${MOTION2[$i]}','cardiac','${CARDIAC2[$i]}','brain','${BRAIN}')"
      fi
    else
      if [ $task == 1 ]; then
        matlab -nojvm -r "addpath ${RTB}; Connectivity_regress('${subject_num}','${region}','${VOI}','motion','${MOTION2[$i]}','cardiac','${CARDIAC2[$i]}','partials','${PARTIAL}','task','${taskfile}')"
        matlab -nojvm -r "addpath ${RTB}; Connectivity_regress('${subject_num}','${region}','${VOI}','motion','${MOTION2[$i]}','cardiac','${CARDIAC2[$i]}','partials','${PARTIAL}','brain','${BRAIN}','task','${taskfile}')"
      else
        matlab -nojvm -r "addpath ${RTB}; Connectivity_regress('${subject_num}','${region}','${VOI}','motion','${MOTION2[$i]}','cardiac','${CARDIAC2[$i]}','partials','${PARTIAL}')"
        matlab -nojvm -r "addpath ${RTB}; Connectivity_regress('${subject_num}','${region}','${VOI}','motion','${MOTION2[$i]}','cardiac','${CARDIAC2[$i]}','partials','${PARTIAL}','brain','${BRAIN}')"
      fi
    fi
    i=$[$i+1]
  done
}

run_correlation () {
PARTIAL=`echo $PARTIAL | awk -F/ '{print $NF}' | awk -F"_rest" '{print $1}'` 
if [ $cardiacfile == dummy.txt ] && [ $partialcorr == 0 ]; then
      if [ $task == 1 ]; then
        cmd=`echo "resting_correlation.sh -s ${subject_num} -S ${DIR} -R ${region} --files ${1} -c motion_wm_vent_deriv_task"`
        cmd1=`echo "resting_correlation.sh -s ${subject_num} -S ${DIR} -R ${region} --files ${1} -c motion_wm_vent_brain_deriv_task"`
      else
        cmd=`echo "resting_correlation.sh -s ${subject_num} -S ${DIR} -R ${region} --files ${1} -c motion_wm_vent_deriv"`
        cmd1=`echo "resting_correlation.sh -s ${subject_num}  -S ${DIR} -R ${region} --files ${1} -c motion_wm_vent_brain_deriv"`
      fi
elif [ $cardiacfile == dummy.txt ]; then
      if [ $task == 1 ]; then
        cmd=`echo "resting_correlation.sh -s ${subject_num} -S ${DIR} -R ${region} --files ${1} -c motion_wm_vent_$PARTIAL_task" --partial $PARTIAL`
        cmd1=`echo "resting_correlation.sh -s ${subject_num} -S ${DIR} -R ${region} --files ${1} -c motion_wm_vent_brain_deriv_$PARTIAL_task" --partial $PARTIAL`
      else
        cmd=`echo "resting_correlation.sh -s ${subject_num} -S ${DIR} -R ${region} --files ${1} -c motion_wm_vent_deriv_$PARTIAL" --partial $PARTIAL`
        cmd1=`echo "resting_correlation.sh -s ${subject_num} -S ${DIR} -R ${region} --files ${1} -c motion_wm_vent_brain_deriv_$PARTIAL" --partial $PARTIAL`
      fi
elif [ $partialcorr == 0 ]; then
      if [ $task == 1 ]; then
        cmd=`echo "resting_correlation.sh -s ${subject_num} -S ${DIR} -R ${region} --files ${1} -c motion_deriv_cardiac_task"`
        cmd1=`echo "resting_correlation.sh -s ${subject_num} -S ${DIR} -R ${region} --files ${1} -c motion_brain_deriv_cardiac_task"` 
      else
        cmd=`echo "resting_correlation.sh -s ${subject_num} -S ${DIR} -R ${region} --files ${1} -c motion_deriv_cardiac"`
        cmd1=`echo "resting_correlation.sh -s ${subject_num} -S ${DIR} -R ${region} --files ${1} -c motion_brain_deriv_cardiac"`
      fi
else
      if [ $task == 1 ]; then
        cmd=`echo "resting_correlation.sh -s ${subject_num} -S ${DIR} -R ${region} --files ${1} -c motion_deriv_$PARTIAL_cardiac_task" --partial $PARTIAL` 
        cmd1=`echo "resting_correlation.sh -s ${subject_num} -S ${DIR} -R ${region} --files ${1} -c motion_brain_deriv_$PARTIAL_cardiac_task" --partial $PARTIAL`
      else
        cmd=`echo "resting_correlation.sh -s ${subject_num} -S ${DIR} -R ${region} --files ${1} -c motion_deriv_$PARTIAL_cardiac" --partial $PARTIAL`
        cmd1=`echo "resting_correlation.sh -s ${subject_num} -S ${DIR} -R ${region} --files ${1} -c motion_brain_deriv_$PARTIAL_cardiac" --partial $PARTIAL`
      fi
fi
if [ $separate == 1 ]; then
   cmd=`echo ${cmd} " --separate"`
   cmd1=`echo ${cmd1} " --separate"`
   ${cmd}
   ${cmd1}
else
   ${cmd}
   ${cmd1}
fi
}

#################################
# Parse Command-line arguements 
#################################
if [ $# == 0 ];then help_text;exit 1;fi

while [ $# -gt 0 ]; do
  case "$1" in

      -S | --study    ) # Required! Sets Study Directory
                        DIR=$2;shift;shift;;

      -s | --subject  ) # Required! Set Subject Number
                        subject_num=$2;shift;shift;;

      -V | --VOI      ) # Required! Sets VOI file
                        voi=$2;shift;shift;;

      -R | --REGION   ) # Required! Adds REGION NAME to output files.
	                region=$2;shift;shift;;

      --lowpass       ) # Optional. Changes the default low pass filter from .08
                        lowpass=$2;shift;shift;;

      --highpass      ) # Optional. Changes the default high pass filter from .009
                        highpass=$2;shift;shift;;

      --motionfiles   ) # Optional. Changes the default string for motion files from rp*txt
                        motion=`echo "${2}"`;shift;shift;;
   
      --files         ) # Requires. Sets search term for files (e.g. rest)
                        files=$2;shift;shift;;

      --regonly       ) # Optional! Program will only do the regression computations. 
                        regonly=1;shift;;

      --cardiac       ) # Optional! Use the cardiac file search string after the flag (e.g. cardiac*${files}*.txt)
                        cardiacfile=`echo "${2}"`;shift;shift;;

      --partial       ) # Optional! Compute the partial correlations. Use all VOI specified in the list.
                        partialcorr=1;partiallist=$2;shift;shift;;

      --LPF           ) # Optional!!! Process the LPF files with resting_correlation in addition to the BPF files.
                        LPF=1;shift;;

      --separate      ) # Optional!!! In addition to processing the files jointly, compute maps for each run on its own within resting_correlation.sh
                        separate=1;shift;;

      --task          ) # Optional! NOT RECOMMENDED!!!!!! If you are inclined to use task data, then use to regress out the effect of task.
                        #After the flag, enter the name of the file that contains the convolved/FIR/basis function task regressors.
                        task=1;taskfile=$2;shift;shift;;

      -h | --help     ) # Get this Help message.
                        help_text;exit 0;;
     
      *               ) # Undefined Parameter
                        echo "ERROR: Unrecognized option: $1";help_text;exit 1;;
      esac
done



##################
# Error Checking 
##################
if which fslroi > /dev/null; then
    echo "" > /dev/null
else
    echo "fslroi does not exist and is required"
    echo "Program will exit."
    exit -1
fi
if which fslstats > /dev/null; then
    echo "" > /dev/null
else
    echo "fslstats does not exist and is required"
    echo "Program will exit."
    exit -1
fi
if which afni > /dev/null; then
    echo "" > /dev/null
else
    echo "afni does not exist and is required"
    echo "Program will exit."
    exit -1
fi
if which matlab > /dev/null; then
    echo "" > /dev/null
else
    echo "matlab does not exist and is required"
    echo "Program will exit."
    exit -1
fi

if [ -z "$files" ]; then
   echo "Files not set "
   help_text
   exit 1
fi
if [ -z "$DIR" ]; then
   echo "Directory not set "
   help_text
   exit 1
fi
if [ -z "$subject_num" ]; then
   echo "Subject number not set "
   help_text
   exit 1
fi
if [ -z "$voi" ]; then
   echo "VOI file not set "
   help_text
   exit 1
fi
if [ -z "$region" ]; then
   echo "Region name not set"
   help_text
   exit 1
fi

if ! [ -d ${DIR}/${subject_num}${subject_dir} ]; then
    echo "Subject directory does not exist"
    help_text;exit 1
fi
if ! [ -e $voi ];then
  help_text;exit 1
fi
if ls ${DIR}/${subject_num}${subject_dir}/${motion} > /dev/null; then
   echo "Motion files found"
else
   echo "Motion files do not exist, check name format and location of files";
   exit
fi

if [ `ls ${DIR}/${subject_num}${subject_dir}/sw*${files}*.[nh][id][ir] | wc -l` -gt 0 ]; then
   echo "Imaging Files found"
else
   echo "Files doe not exist, check name format and location of files"
   exit
fi
if [ ${partialcorr} == 1 ]; then
   if [ -e $partiallist ]; then
      for i in `more ${partiallist}`
      do
          if ! [ -e $i ]; then
             echo "Partial Correlation ROI does not exist";
	     exit
	  fi
      done
   else
      echo "Partial Correlation list does not exist"
      exit
   fi
fi
if [ "$cardiacfile" != "dummy.txt" ]; then
   if ls ${DIR}/${subject_num}${subject_dir}/${cardiac}; then
      echo "Cardiac files found"
   else
      echo "Cardiac files do not exist, check name format and location of files";
      exit
   fi
fi

echo "Preprocessed Directory is:" ${DIR}${subject_num}${subject_dir}

##########################
# Create Rest regressors 
##########################
cd ${DIR}/${subject_num}${subject_dir}
FILES=(`ls sw*${files}*.[nh][id][ir]`)
if ! [ -d images_rest ]; then
   mkdir images_rest
   cd images_rest
else
   cd images_rest
fi
dirpath=`pwd`

if [ $regonly == 1 ]; then
   run_regress BPF
   run_correlation PSC_BPF${lowpass}${highpass}
   if [ $LPF == 1 ]; then
      run_regress LPF
      run_correlation PSC_LPF{lowpass}
   fi
   exit
fi

for i in `echo ${FILES[@]}`
do
  ln -s ../$i $i
done

CNT=${#FILES[@]}
i="0"
if ! [ -e tmpmask.nii ]; then
   fslroi ${FILES[0]} tmpmask.nii 0 1
fi
sf=`fslstats tmpmask.nii.gz -P 95`

if ! [ -e ${subject_num}_WM.nii.gz ] ; then
  3dcalc -a ${ROIDIR}/wm.nii -b ${FILES[$i]}[0] -expr "(1-astep(b/$sf,.95))*step(a)" -prefix ${subject_num}_WM.nii.gz -datum short
fi
if ! [ -e ${subject_num}_VENT.nii.gz ]; then
  3dcalc -a ${ROIDIR}/ventricles.nii -b ${FILES[$i]}[0] -expr "(1-astep(b/$sf,.75))*step(a)" -prefix ${subject_num}_VENT.nii.gz -datum short
fi


if [ $LPF == 1 ]; then
  while [ $i -lt $CNT ]
  do
    FILE=`echo ${FILES[$i]} | awk -F. '{print $1}' `
    if ! [ -e mean_${FILE}.nii ] && ! [ -e mean_${FILE}.nii.gz ]; then
     3dTstat -prefix mean_${FILE}.nii.gz -datum float ${FILES[$i]}
    fi
    if ! [ -e LPF${lowpass}_${FILE}.nii ] && ! [ -e LPF${lowpass}_${FILE}.nii.gz ]; then
      echo "Begin Processing: " ${FILES[$i]}
      3dFourier -lowpass ${lowpass} -prefix LPF${lowpass}_${FILE}.nii.gz ${FILES[$i]} -datum float
      3dcalc -a LPF${lowpass}_${FILE}.nii.gz -b mean_${FILE}.nii.gz -expr 'a/b*100' -prefix PSC_LPF${lowpass}_${FILE}.nii.gz -datum float
    fi
    i=$[$i+1]
  done
fi

i="0"
while [ $i -lt $CNT ]
do 
  FILE=`echo ${FILES[$i]} | awk -F. '{print $1}' `
  if ! [ -e mean_${FILE}.nii ] && ! [ -e mean_${FILE}.nii.gz ]; then
     3dTstat -prefix mean_${FILE}.nii.gz -datum float ${FILES[$i]}
  fi
  if ! [ -e BPF${lowpass}${highpass}_${FILE}.nii ] && ! [ -e BPF${lowpass}${highpass}_${FILE}.nii.gz ]; then
     3dFourier -lowpass ${lowpass} -highpass ${highpass} -prefix BPF${lowpass}${highpass}_${FILE}.nii.gz ${FILES[$i]} -datum float
     3dcalc -a BPF${lowpass}${highpass}_${FILE}.nii.gz -b mean_${FILE}.nii.gz -expr 'a/b*100' -prefix PSC_BPF${lowpass}${highpass}_${FILE}.nii.gz -datum float
  fi
  i=$[$i+1]
done

gzip *LPF${lowpass}_sw*.nii
gzip *BPF${lowpass}${highpass}_sw*.nii

run_regress BPF
run_correlation PSC_BPF${lowpass}${highpass}

if [ $LPF == 1 ]; then
   run_regress LPF
   run_correlation PSC_LPF${lowpass}
fi

exit