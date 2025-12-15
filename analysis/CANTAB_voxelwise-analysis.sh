#!/bin/bash

## directories and mask
project_dir="/BICNAS2/tuominen/FEOBV-pet"
gr_analyses_dir="${project_dir}/PET_group_analyses"
PET_data_dir="${gr_analyses_dir}/final_SUVR-in-MNI"
vw_mask="${project_dir}/atlases/MNI152-2mm__full-voxelwise_MASK_finality.nii.gz"

## GLM-related variables
####################################################################
sample="18patients-CANTAB"
contrast="contrast_g1v1_slope"
GLM="18patients-CANTAB--FINALITY-nLin_precise_SUVR-in-MNI-2mm__contrast_g1v1_slope"
SUVR_prefix="nLin_precise_SUVR-in-MNI-2mm"
FSGD="${gr_analyses_dir}/18patients-CANTAB--FINALITY-nLin_precise_SUVR-in-MNI-2mm__contrast_g1v1_slope.fsgd"
CONTRAST="${gr_analyses_dir}/contrast_g1v1_slope.mtx"
####################################################################

## Define subjects' list
readarray -t subjects < ${gr_analyses_dir}/${GLM}.subjectlist.dat

echo "============================================================================================="
echo "   Subjects:"
echo "   ${subjects[@]} "
echo "============================================================================================="

PET_data=()

for subj in "${subjects[@]}"; do
	PET_data+=("${subj}__sm8-${SUVR_prefix}.nii.gz")
done

concatenated_PET_data="concatenated__${GLM}__sm8-finality-${SUVR_prefix}.nii.gz"
#### concatenate PET data
cd ${PET_data_dir}
fslmerge -t ${concatenated_PET_data} "${PET_data[@]}"
	echo ""
	echo ""

### set GLM dir
cd ${gr_analyses_dir}
GLM_dir="${gr_analyses_dir}/GLMfit__${GLM}_sm8_output"

## archive if already exists
if [ -d "${GLM_dir}" ]; then
	mod_time=$(stat -c %y "${GLM_dir}" | cut -d'.' -f1 | sed 's/ /_/; s/:/h/g' | cut -c1-15)
	mv "${GLM_dir}" "${GLM_dir}__${mod_time}"
fi

# Run GLM on
echo ">>> Running GLM"
echo ""
mri_glmfit \
	--y ${PET_data_dir}/${concatenated_PET_data} \
	--fsgd ${gr_analyses_dir}/${GLM}.fsgd \
	--C ${gr_analyses_dir}/${contrast}.mtx \
	--glmdir ${GLM_dir} \
	--mask ${vw_mask} \
	--eres-save \
	--nii.gz
