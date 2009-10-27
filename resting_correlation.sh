#!/bin/bash

###########################
# Resting Correlation Script v2
# 10/8/09 
# Donald McLaren (dm10@medicine.wisc.edu) & Kim Farbota (farbota@wisc.edu)
# Wisconsin Alzheimer's Disease Center Imaging Core (www.brainmap.wisc.edu)
#
# USAGE: 
#       resting_correlation.sh -s 2041 -R DMPFC -S /Data/vtrak1/epilepsy -f BPF.008.09_dsw -c motion_brain_deriv_partials_cardiac
#       resting_correlation.sh --help will list all possible options and an explanation
#
# This script will compute the residuals for the whole brain (in file 'resid*') and then it will compute the correlation between the residuals 
# in the seed region with the voxel-wise residuals (in file 'rsqrd*'). This produces an R^2 map. The correlation coefficent is computed by taking
# the square root of the R^2 map and multiplying by the sign of the slope (in file 'corr*'). Then the map is made to follow a Gaussian distribution 
# with the Fisher's z' Transformation (inverse hyperbolic tangent -- in file 'normed_corr*'). The slope of the line regression 
# is also output to a file called 'beta*'. NOTE: The 'beta' value is only interpretable if the data was already converted to percent signal change. 
#
# The output directory for these files is fmri_resting_stats/stats_${region}_${files}_${covariates} within the subject directory in the default options below. 
# If the directory does not exist, it will be created by the script. This output directory, if it exists, must be empty.
#
# One must configure the default options below for their system (e.g. ROI dir and subject_dir).
#
# This script is called from resting_connectivity.sh
#
###########################

##########################
# Default Options 
##########################
PATH=$PATH:/Data/data1/lab_scripts/RestingToolbox  # Add the RestingToolbox directory to the environment path
subject_dir=/fmri                                  # Directory for fmri data within each subject, must begin with  a /
ROIDIR=/Data/data1/ROI/mclaren                     # Directory containing the EPI mask -- rEPI.nii
subject_num=0 	                                   # Error Check Dummy
region=0  	                                   # Error Check Dummy
files=0     	                                   # Error Check Dummy
covariates=0                                       # Error Check Dummy
study=0		                                   # Error Check Dummy
separate=0                                         # Default option
partials=dummy                                     # Error Check Dummy

###########################
# Functions
###########################
help_text () {
echo "Available Options:"
echo " "
echo "      -s | --subject  ) # Required! Set Subject Number"
echo " "
echo "      -R | --REGION    ) # Required! Adds REGION NAME to output files."
echo " "
echo "      -S | --Study    ) # Required! Sets the study directory."
echo " "
echo "      -f | --files   ) # Required! Sets the regressor file to be the properly filtered data"
echo " "   
echo "      -c | --covariates   ) # Required! Sets the covariate search string."
echo " "
echo "      --partial    ) # Optional! Compute the partial correlations. Use all VOI specified in the list. Use the base of the filename (e.g. partial instead of partial.txt)"
echo " "
echo "      --separate      ) # Optional!!! In addition to processing the files jointly, compute maps for each run on its own. "
echo " "
echo "      -h | --help     ) # Get this Help message."
}

run_filecheck () {
if [ `ls ${study}/${subject_num}${subject_dir}/images_rest/${subject_num}_${covariates}*${files}*${region}*.txt | wc -l` -gt 0 ]; then
        if ! [ `ls ${study}/${subject_num}${subject_dir}/images_rest/${subject_num}_residuals_${covariates}*${files}*${region}*.txt | wc -l` -gt 0 ]; then
           echo " Seed Region Residual File does not exist."
           echo " Program will now exit."
           exit -1
	fi
else
        echo " Covariate file does not exist."
        echo " Program will now exit."
        exit -1
fi
}


#################################
# Parse Command-line arguements 
#################################
if [ $# == 0 ];then help_text;exit 1;fi
 
while [ $# -gt 0 ]; do
  case "$1" in

      -s | --subject  ) # Required! Set Subject Number
                        subject_num=$2;shift;shift;;

      -R | --REGION    ) # Required! Adds REGION NAME to output files.
	                region=$2;shift;shift;;

      -S | --Study    ) # Required! Sets the study name.
                        study=$2;shift;shift;;

      -f | --files   ) # Required! Sets the regressor file to be the properly filtered data
                        files=$2;shift;shift;;
     
      -c | --covariates ) # Required! Specifies which covariates to use.
			covariates=$2;shift;shift;;
      
      --partial       ) # Optional! Compute the partial correlations. Use all VOI specified in the list.
                       partials=$2;shift;shift;;

      --separate      ) # Optional!!! In addition to processing the files jointly, compute maps for each run on its own. 
                        separate=1;shift;;

      -h | --help     ) # Get this Help message.
                        help_text;exit 0;;

      *               ) # Undefined Parameter
                        echo "ERROR: Unrecognized option: $1";help_text;exit 1;;
      esac
done

##################
# Error Checking #
##################
echo "Processing Subject: " ${subject_num}

if [ $study == 0 ];then
  	echo "Study Directory Not Listed"
	help_text;exit 1
elif ! [ -d $study ]; then
        echo "Study directory:" $study " does not exist"
        help_text;exit -1   
fi
if [ $subject_num == 0 ];then
  	echo "Subject Number Missing"
	help_text;exit 1
elif ! [ -d ${study}${subject_num}${subject_dir}/images_rest ]; then
        echo "Subject directory:" ${study}/${subject_num}${subject_dir}/images_rest " does not exist"
        help_text;exit -1
fi
if [ $files == 0 ];then
  	echo "Files Not Specified"
	help_text;exit 1
elif ! [ `ls ${study}${subject_num}${subject_dir}/images_rest/*${files}*nii.gz | wc -l` -gt 0 ]; then
        echo "Files: " ${study}${subject_num}${subject_dir}/images_rest/*${files}*.nii.gz " do not exist."
        help_text;exit -1
fi
if [ $region == 0 ];then
  	echo "Region Not Specified"
	help_text;exit 1
fi

echo $covariates

if [ $covariates == 0 ];then
  	echo "Covariates Missing"
	help_text;exit 1
elif [ "$covariates" == "motion_wm_vent_brain_deriv_task" ]; then
    echo "continuing using covariates" 
elif [ "$covariates" == "motion_wm_vent_deriv_task" ]; then   
    echo "continuing using covariates"
elif [ "$covariates" == "motion_wm_vent_brain_deriv" ]; then
    echo "continuing using covariates"
elif [ "$covariates" == "motion_wm_vent_deriv" ]; then
    echo "continuing using covariates"
elif [ "$covariates" == "motion_wm_vent_brain_deriv_${partials}" ]; then
    echo "continuing using covariates"
elif [ "$covariates" == "motion_wm_vent_deriv_${partials}" ]; then
    echo "continuing using covariates"
elif [ "$covariates" == "motion_wm_vent_brain_deriv_${partials}_task" ]; then
    echo "continuing using covariates"
elif [ "$covariates" == "motion_wm_vent_deriv_${partials}_task" ]; then
    echo "continuing using covariates"
elif [ "$covariates" == "motion_brain_deriv_cardiac" ]; then
    echo "continuing using covariates"
elif [ "$covariates" == "motion_deriv_cardiac" ]; then
    echo "continuing using covariates"
elif [ "$covariates" == "motion_brain_deriv_cardiac_task" ]; then
    echo "continuing using covariates"
elif [ "$covariates" == "motion_deriv_cardiac_task" ]; then
    echo "continuing using covariates"
elif [ "$covariates" == "motion_brain_deriv_${partials}_cardiac" ]; then
    echo "continuing using covariates"
elif [ "$covariates" == "motion_deriv_${partials}_cardiac" ]; then
    echo "continuing using covariates"
elif [ "$covariates" == "motion_brain_deriv_${partials}_cardiac_task" ]; then
    echo "continuing using covariates"
elif [ "$covariates" == "motion_deriv_${partials}_cardiac_task" ]; then
    echo "continuing using covariates"
else
	echo "Covariates are not correct"
        echo "Covariates are not correct" > ERROR.log
	help_text; exit 1
fi

run_filecheck

######################
# Program begins here 
######################

if ! [ -d ${study}${subject_num}/fmri_resting_stats/stats_${files}_${region}_${covariates} ]; then
   mkdir -p ${study}${subject_num}/fmri_resting_stats/stats_${files}_${region}_${covariates}
fi   
cd ${study}/${subject_num}/fmri_resting_stats/stats_${files}_${region}_${covariates}
if [ "$(ls -A ${study}/${subject_num}/fmri_resting_stats/stats_${files}_${region}_${covariates})" ]; then
  a=`pwd`
  echo "Current directory: " $a
  echo "Do you want to remove all the files in this directory?"
  read response
  if [ $response == "y"  ] | [ $response == "yes" ]; then
    rm *
  else
    echo " ERROR: Directory not empty" > ERROR.log
    echo " The directory is not empty. Please move the files and rerun."
    echo " Program will now exit."
    exit -1
  fi
fi

pushd ${study}/${subject_num}/${subject_dir}/images_rest 
  FILES=(`lspath "${files}*"`)
  FILES_COV=(`lspath "${subject_num}_${covariates}*${files}*${region}.txt"`)
  FILES_RESID=(`lspath "${subject_num}_residuals_${covariates}_rest*${files}*${region}.txt"`)
popd


#####################
# Compute Voxel-wise Residuals
#####################
CNT=${#FILES[@]}
i="0"
stim=`awk 'END{print NF}' ${FILES_COV[$i]}`
while [ $i -lt $CNT ]
do 
   if ! [ -e resid_${subject_num}_${files}_${covariates}.${i}.nii.gz ]; then
      	echo "3dDeconvolve -input "${FILES[$i]}" -polort 1 -num_stimts " ${stim} "\\" > 3dDecon.${i}.sh
      	j="0"
       	while [ $j -lt ${stim} ]
        do
		echo "-stim_file " $[$j+1] " " ${FILES_COV[$i]}[$j] "\\" >> 3dDecon.${i}.sh
		j=$[$j+1]
	done
       	echo "-mask "${ROIDIR}"/rEPI.nii -errts resid_"${subject_num}"_"${files}"_"${covariates}"."${i}".nii.gz -nobucket" >> 3dDecon.${i}.sh
   fi
   chmod 770 3dDecon.${i}.sh
   ./3dDecon.${i}.sh 
   i=$[$i+1]
done

######################
# Remove Old Files
######################
if [ `ls corr_*.nii* | wc -l ` -gt 0 ]; then
   rm corr_*.nii*
fi
if [ `ls Z_*.nii* | wc -l ` -gt 0 ]; then
   rm Z_*.nii*
fi
if [ `ls normed_corr_*.nii* | wc -l ` -gt 0 ]; then
   rm normed_corr_*.nii*
fi
if [ `ls rsqrd_* | wc -l ` -gt 0 ]; then
   rm rsqrd_*
fi

######################
# Compute Voxel-wise Correlations and Related Values
######################
if [ -e ${subject_num}_${files}_${covariates}_stim_file.1D ]; then
   rm ${subject_num}_${files}_${covariates}_stim_file.1D 
fi
cat ${FILES_RESID[*]} >> ${subject_num}_residuals_${files}_${covariates}_stim_file.1D

if [ $separate == 1 ]; then
   i="0"
   while [ $i -lt $CNT ]
   do 
      3dDeconvolve -input resid_${subject_num}_${files}_${covariates}.${i}.nii.gz -polort 0 -num_stimts 1 \
      -stim_file 1 ${FILES_RESID[$i]}\
      -mask ${ROIDIR}/rEPI.nii -rout -tout -bucket rsqrd_${subject_num}_${files}_${region}_${covariates}.${i}.nii.gz 
      
     3dcalc -a rsqrd_${subject_num}_${files}_${region}_${covariates}.${i}.nii.gz[0] -b rsqrd_${subject_num}_${files}_${region}_${covariates}.${i}.nii.gz[2] -expr '(sqrt(a)*ispositive(b))-(sqrt(a)*isnegative(b))' -prefix corr_${subject_num}_${files}_${region}_${covariates}.${i}.nii.gz
     3dcalc -a corr_${subject_num}_${files}_${region}_${covariates}.${i}.nii.gz -expr 'atanh(a)' -prefix normed_corr_${subject_num}_${files}_${region}_${covariates}.${i}.nii
     3dcalc -a rsqrd_${subject_num}_${files}_${region}_${covariates}.${i}.nii.gz[2] -expr 'a' -prefix beta_${subject_num}_${files}_${region}_${covariates}.${i}.nii
     
     # Set the variance term appropriately for computing the Z at the single subject level. Since df is constant within a study, the correction doesn't effect group analyses. 
     #F=`awk 'END{print NR}' ${FILES_RESID[$i]}`
     #df=`echo "sqrt(1/(("$F"/2.34)-3))" | bc -l `
     #3dcalc -a normed_corr_${subject_num}_${files}_${region}_${covariates}.${i}.nii -expr a/${df} -prefix Z_correlation_${subject_num}_${files}_${region}_${covariates}.${i}.nii
     i=$[$i+1]
   done
fi

3dDeconvolve -input resid_${subject_num}_${files}_${covariates}.*.nii.gz -polort 0 -num_stimts 1 \
-stim_file 1 ${subject_num}_residuals_${files}_${covariates}_stim_file.1D \
-mask ${ROIDIR}/rEPI.nii -rout -tout -bucket rsqrd_${subject_num}_${files}_${region}_${covariates}.nii.gz 
      
3dcalc -a rsqrd_${subject_num}_${files}_${region}_${covariates}.nii.gz[0] -b rsqrd_${subject_num}_${files}_${region}_${covariates}.nii.gz[2] -expr '(sqrt(a)*ispositive(b))-(sqrt(a)*isnegative(b))' -prefix corr_${subject_num}_${files}_${region}_${covariates}.nii.gz
3dcalc -a corr_${subject_num}_${files}_${region}_${covariates}.nii.gz -expr 'atanh(a)' -prefix normed_corr_${subject_num}_${files}_${region}_${covariates}.nii
3dcalc -a rsqrd_${subject_num}_${files}_${region}_${covariates}.nii.gz[2] -expr 'a' -prefix beta_${subject_num}_${files}_${region}_${covariates}.nii

# Set the variance term appropriately for computing the Z at the single subject level. Since df is constant within a study, the correction doesn't effect group analyses. 
#F=`awk 'END{print NR}' ${subject_num}_residuals_${covariates}_${files}_stim_file.1D`
#df=`echo "sqrt(1/(("$F"/2.34)-3))" | bc -l `
#3dcalc -a normed_corr_${subject_num}_${files}_${region}_${covariates}.nii -expr a/${df} -prefix Z_correlation_${subject_num}_${files}_${region}_${covariates}.nii

#Remove intermediate files
rm resid*

exit