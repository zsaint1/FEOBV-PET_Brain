#!/bin/bash

## directories and mask
project_dir="/BICNAS2/tuominen/FEOBV-pet"
gr_analyses_dir="${project_dir}/PET_group_analyses"
vw_mask="${project_dir}/atlases/MNI152-2mm__full-voxelwise_MASK_finality.nii.gz"

cd ${gr_analyses_dir}

## GLM-related variables
######################################################################################
GLM="18patients-33controls--FINALITY-nLin_precise_SUVR-in-MNI-2mm__contrast_g2v0"
######################################################################################

## start time
start_time=$(date +%s)

echo -e "\n\e[32;1m ******* Commencing Permutations ******* \e[0m \n"

## define GLM and Permutations directories
GLM_dir="${gr_analyses_dir}/GLMfit__${GLM}_sm8_output"
GLM_perm_dir="${GLM_dir}___permutations_10000-abs--vwthr1.3"

## Exits if uncorrected voxelwise dir does not exists
if [ ! -d "${GLM_dir}" ]; then
	echo -e "\e[31;1m !?! ERROR: Uncorrected voxelwise model does not exist: ${GLM_dir} \e[0m"
		exit
fi

## archive if already exists
if [ -d "${GLM_perm_dir}" ]; then
	mod_time=$(stat -c %y "${GLM_perm_dir}" | cut -d'.' -f1 | sed 's/ /_/; s/:/h/g' | cut -c1-15)
	mv "${GLM_perm_dir}" "${GLM_perm_dir}__${mod_time}"
fi

## copy GLM dir for new permutations
cp -rv ${GLM_dir} ${GLM_perm_dir}
cd ${gr_analyses_dir}

### run Permutations
mri_glmfit-sim \
	--overwrite \
	--glmdir ${GLM_perm_dir} \
	--mask /BICNAS2/tuominen/FEOBV-pet/atlas-for-ROIs/MNI152-2mm__full-voxelwise_MASK_finality.nii.gz \
	--perm 10000 1.3 abs \
	--cwp .05 \
	--centroid \
	--a2009s \
	--bg 45	## adjust to # of cores available


echo ""

if [ -e ${GLM_perm_dir}/contrast_g2v0/sig.nii.gz ] ; then
	echo "=========================================================================================================================================================================="
	echo -e "\e[32;1m *** Permutations results saved here:"
	echo -e "${GLM_perm_dir}/contrast_g2v0/ \e[0m"
	echo "=========================================================================================================================================================================="
fi

## time elapsed
end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
elapsed_formatted=$(printf '%02dh:%02dm:%02ds' $((elapsed/3600)) $(( (elapsed%3600)/60 )) $((elapsed%60)))
echo ""
echo "#########################################################"
echo -e "############# \e[32;1m  Time Elapsed: ${elapsed_formatted} \e[0m  #############"
echo -e "############# \e[32;1m All permutations completed. \e[0m #############"
echo "#########################################################"
