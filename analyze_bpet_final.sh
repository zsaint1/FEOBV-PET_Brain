#!/bin/bash

## verify if $SUBJECTS_DIR is set
if [ -z "${SUBJECTS_DIR}" ]
	then
		echo "------------------------------------------------------------------------------"
		echo "   ERROR :::  \$SUBJECTS_DIR is not set. Please set it before proceeding."
		echo "------------------------------------------------------------------------------"
			exit 1
fi

## define general directories
project_dir="$(dirname "${SUBJECTS_DIR}")"
mc_qc_dir="${project_dir}/QC/motion"
ref_qc_dir="${project_dir}/QC/reference_region"
atlas_dir="${project_dir}/atlas-for-ROIs"
stats_files_dir="${project_dir}/results/ROI_stats_files"
PET_data_dir="${project_dir}/PET_group_analyses/final_SUVR-in-MNI"
# scripts_dir="${project_dir}/scripts"
scripts_dir="/group/tuominen/FEOBV/scripts"
pscript="/group/tuominen/anaconda3/bin/python"


## define outer function
analyze_pbet_function () {

	subj=${1}
	subj_=${subj}
	echo "
====================================================================================================================================
===================================================================================================================================="
	echo -e " \e[33;1m

 ${subj_}${subj_}${subj_}  ${subj_}${subj_}           ${subj_}            FEOBV${subj_}${subj_}    ${subj_}                  ${subj_}
 ${subj_}${subj_}${subj_}  ${subj_}${subj_}       ${subj_}${subj_}        ${subj_}${subj_}${subj_}  ${subj_}                ${subj_}
 ${subj_}${subj_}${subj_}  ${subj_}${subj_}     ${subj_}     ${subj_}     ${subj_}         ${subj_}  ${subj_}              ${subj_}
 ${subj_}                  ${subj_}            ${subj_}       ${subj_}    ${subj_}          ${subj_}  ${subj_}            ${subj_}
 ${subj_}                  ${subj_}           ${subj_}         ${subj_}   ${subj_}         ${subj_}    ${subj_}          ${subj_}
 ${subj_}${subj_}          ${subj_}${subj_}  ${subj_}           ${subj_}  ${subj_}${subj_}${subj_}      ${subj_}        ${subj_}
 ${subj_}${subj_}          ${subj_}${subj_}  ${subj_}           ${subj_}  ${subj_}         ${subj_}      ${subj_}      ${subj_}
 ${subj_}                  ${subj_}           ${subj_}         ${subj_}   ${subj_}          ${subj_}      ${subj_}    ${subj_}
 ${subj_}                  ${subj_}            ${subj_}       ${subj_}    ${subj_}          ${subj_}       ${subj_}  ${subj_}
 ${subj_}                  ${subj_}${subj_}     ${subj_}     ${subj_}     ${subj_}         ${subj_}         ${subj_}${subj_}
 ${subj_}                  ${subj_}${subj_}       ${subj_}${subj_}        ${subj_}${subj_}${subj_}           ${subj_}FEOBV
 ${subj_}                  ${subj_}${subj_}           ${subj_}            FEOBV${subj_}${subj_}                FEOBVFEOBV

\e[0m"
	echo "
====================================================================================================================================
===================================================================================================================================="

	sleep 1

	## define relevant subdirectories for the subject
	PET_dir="${project_dir}/subjects/${subj}/pet/001"
	ref_dir="${PET_dir}/reference_region"
	recon_dir="/BICNAS2/tuominen/FEOBV-brain-PET-recons/2min-frames-no-filter/BRAIN-FEOBV-${subj:5}-Dyn-3i21ss-pseudoCT-0mmGausss-z1-172IM-2minframes-OP-DICOM"
	mri_dir="${SUBJECTS_DIR}/${subj}_pet/mri"
	T1_dir="${project_dir}/subjects/${subj}/anat/pet_T1"
	mc_dir="${PET_dir}/motion_correction"
	dynamic_qc_dir="${PET_dir}/QC_dynamic_PET"
	reg_mask_dir="${PET_dir}/registration_mask"


	## convert PET data to nifti (if not present in PET_dir)
	if [ ! -e ${PET_dir}/PET.nii ] ; then
		echo "=============================================================================================="
		echo -e " \e[33;1m ${subj} --> Converting PET dicoms to nifti (PET.nii). \e[0m"
		echo "=============================================================================================="
		echo ""
			${scripts_dir}/dependencies/mri_convert_for_petdata -i ${recon_dir}/*0001.ima -o ${PET_dir}/PET.nii
	else
		echo -e "\e[32;1m || ${subj} --> PET data were already converted ('PET.nii'). \e[0m"
	fi


	## recon-all  ### can modify to check for the recon-all_w_FS_8.0.0--DONE
	if ! grep -q "finished without error" "${SUBJECTS_DIR}/${subj}_pet/scripts/recon-all-status.log" 2>/dev/null; then
		if [ -e ${mri_dir}/orig.mgz ]; then
			echo "=============================================================================================="
			echo -e " \e[33;1m ${subj} --> RECON-ALL was previously started but did not complete without error. \e[0m"
			echo "      Proceeding with RECON-ALL."
			echo "=============================================================================================="
			echo ""
				recon-all -s ${subj}_pet -all -parallel -no-isrunning -notify ${T1_dir}/recon-all_w_FS_8.0.0--DONE
		else
			echo "=============================================================================================="
			echo -e " \e[33;1m ${subj} --> RECON-ALL was never done. Proceeding... \e[0m"
			echo "=============================================================================================="
			echo ""
				recon-all -i ${T1_dir}/T1.nii -s ${subj}_pet -all -parallel -notify ${T1_dir}/recon-all_w_FS_8.0.0--DONE
		fi
	else
		echo -e "\e[32;1m || ${subj} --> RECON-ALL was already completed. \e[0m"
	fi


	if [[ ! -f "${mri_dir}/sclimbic.mgz" ]]; then
		echo "=============================================================================================="
		echo -e " \e[33;1m ${subj} --> Running ScLimbic (RECON-ALL) module. \e[0m"
		echo "=============================================================================================="
		echo ""
		mri_sclimbic_seg -s "${subj}_pet"
	else
		echo -e "\e[32;1m || ${subj} --> mri_sclimbic_seg already performed. \e[0m"
	fi


	if [[ ! -f "${reg_mask_dir}/synthseg_anat.nii.gz" ]]; then
		echo "=============================================================================================="
		echo -e " \e[33;1m ${subj} --> Performing SynthSeg. \e[0m"
		echo "=============================================================================================="
		echo ""
			mri_synthseg \
				--i ${T1_dir}/T1.nii \
				--o ${reg_mask_dir}/synthseg.nii.gz \
				--vol ${T1_dir}/synthseg.volumes.csv \
				--threads 5 \
				--cpu

			############# can transfer these steps to final registration_segmentation section?
			## convert aparc+aseg to nifti
			mri_convert \
				${mri_dir}/aparc.a2009s+aseg.mgz \
				${reg_mask_dir}/aparc.a2009s+aseg.nii.gz \
			## reslice in the typical anatomical space/format
			mri_convert \
				${reg_mask_dir}/synthseg.nii.gz \
				${reg_mask_dir}/synthseg_anat.nii.gz \
				--reslice_like ${reg_mask_dir}/aparc.a2009s+aseg.nii.gz \
				--resample_type nearest
			############# can transfer these steps to final registration_segmentation section?
	else
		echo -e "\e[32;1m || ${subj} --> mri_synthseg already performed ('synthseg_anat.nii.gz'). \e[0m"
	fi


	## motion analysis -> raw (before motion correction) PET
	raw_motion_outliers="${subj}_raw_PET_motion_outliers.dat"
	raw_fd_metric="${subj}_raw_PET_fd_metric.dat"
	raw_fd_plot="${subj}_raw_PET_fd_plot.png"

	if [[ ! -f "${PET_dir}/${raw_fd_metric}" ]]; then
		echo "=================================================================================================================="
		echo -e " \e[33;1m ${subj} --> Motion analysis of uncorrected PET data, retrieving raw motion frame-wise displacements. \e[0m"
		echo "=================================================================================================================="
			fsl_motion_outliers -i ${PET_dir}/PET.nii \
								-o ${PET_dir}/${raw_motion_outliers} \
								-s ${PET_dir}/${raw_fd_metric} \
								-p ${PET_dir}/${raw_fd_plot} \
								--fd -v
				## save to QC dir
				cp -v ${PET_dir}/${raw_fd_metric} ${mc_qc_dir}/.
				cp -v ${PET_dir}/${raw_fd_plot} ${mc_qc_dir}/.
	else
		echo -e "\e[32;1m || ${subj} --> Motion analysis of raw PET data was already completed. \e[0m"
	fi


	## Motion correction
	if [[ ! -f "${PET_dir}/mc_PET.nii.gz" ]]; then
		echo "=============================================================================================="
		echo -e " \e[33;1m ${subj} --> Conducting motion correction of PET data. \e[0m"
		echo "=============================================================================================="
			motion_correction_function ${subj}
	else
		echo -e "\e[32;1m || ${subj} --> Motion correction of raw PET data was already completed. \e[0m"
	fi


	## Summing motion corrected PET data
	if [[ ! -f "${PET_dir}/sum_mc_PET.nii.gz" ]]; then
		echo "=============================================================================================="
		echo -e " \e[33;1m ${subj} --> Combining motion-corrected PET data. \e[0m"
		echo "=============================================================================================="
			mri_concat ${PET_dir}/mc_PET.nii.gz --sum --o ${PET_dir}/sum_mc_PET.nii.gz
	else
		echo -e "\e[32;1m || ${subj} --> Summing of motion corrected PET data was already done ('sum_mc_PET.nii.gz'). \e[0m"
	fi


	## Co-register PET and Anatomical T1
	if [[ ! -f "${PET_dir}/NATIVE-to-ANAT.reg.lta" ]]; then
		echo "=============================================================================================="
		echo -e " \e[33;1m ${subj} --> Registering native PET to anatomical T1 space. \e[0m"
		echo "=============================================================================================="
			## calculate co-registration (native to anat transforms)
			mri_coreg \
				--s ${subj}_pet \
				--mov ${PET_dir}/sum_mc_PET.nii.gz \
				--reg ${PET_dir}/NATIVE-to-ANAT.reg.lta
	else
		echo -e "\e[32;1m || ${subj} --> Native PET to anatomical space Registration already exists (NATIVE-to-ANAT.reg.lta). \e[0m"
	fi


	## Move summed PET to Anatomical space
	if [[ ! -f "${PET_dir}/sum_mc_PET-in-ANAT.nii.gz" ]]; then
		echo "=============================================================================================="
		echo -e " \e[33;1m ${subj} --> Projecting summed PET in anatomical space. \e[0m"
		echo "=============================================================================================="
			## move summed (3D) PET to anatomical space
			mri_vol2vol \
				--reg ${PET_dir}/NATIVE-to-ANAT.reg.lta \
				--mov ${PET_dir}/sum_mc_PET.nii.gz \
				--fstarg \
				--o ${PET_dir}/sum_mc_PET-in-ANAT.nii.gz
	else
		echo -e "\e[32;1m || ${subj} --> Projection of summed native PET to anatomical space already done (sum_mc_PET-in-ANAT.nii.gz). \e[0m"
	fi


################################################################################################################


	## Make reference region mask
	if [[ ! -f "${PET_dir}/${subj}_clean_wm-ref-region-mask.nii.gz" ]]; then
		echo "=============================================================================================="
		echo -e " \e[33;1m ${subj} --> Creating reference region mask in anatomical space for ${subj}. \e[0m"
		echo "=============================================================================================="
		echo ""
			create_wm_ref_region_function ${subj}
	else
		echo -e "\e[32;1m || ${subj} --> Reference region mask was already created. \e[0m"
	fi


	## Calculate mean activity in reference region
	if [[ ! -s "${PET_dir}/wm_ref_region_activity" ]]; then
		echo "=============================================================================================="
		echo -e " \e[33;1m ${subj} --> Calculating and saving mean activity in reference region. \e[0m"
		echo "=============================================================================================="
		echo ""
		## Obtaining average activity from the white matter reference region
		fslstats \
			${PET_dir}/sum_mc_PET-in-ANAT.nii.gz \
			-k ${PET_dir}/${subj}_clean_wm-ref-region-mask.nii.gz \
			-M > ${PET_dir}/wm_ref_region_activity

		cp -v \
			${PET_dir}/wm_ref_region_activity \
			${ref_qc_dir}/${subj}_wm_ref_region_activity
	else
		echo -e "\e[32;1m || ${subj} --> Reference region mean activity was already calculated and saved. \e[0m"
	fi


	## Calculate Standardized Uptake Value Ratios (in ANAT)
	if [[ ! -f "${PET_dir}/SUVR-in-ANAT.nii.gz" ]]; then
		echo "=============================================================================================="
		echo -e " \e[33;1m ${subj} --> Calculating Standardized Uptake Value Ratios (SUVRs) in T1 anatomical space. \e[0m"
		echo "=============================================================================================="
		echo ""
		# Load average value of reference region in a variable
		R="$(cat ${PET_dir}/wm_ref_region_activity)"
		## Divide brain PET in ANAT by white matter reference
		fslmaths \
			${PET_dir}/sum_mc_PET-in-ANAT.nii.gz \
			-div $R \
			${PET_dir}/SUVR-in-ANAT.nii.gz
	else
		echo -e "\e[32;1m || ${subj} --> SUVR calculation already performed. \e[0m"
	fi


	## Create precise brainmask and strip T1
	if [[ ! -f "${PET_dir}/ANAT_T1.nii.gz" ]]; then
		echo "======================================================================================================="
		echo -e " \e[33;1m ${subj} --> Creating brainmask and stripping. \e[0m"
		echo "======================================================================================================="
		echo ""
			creating_brainmask_and_stripping__function ${subj}
	else
		echo -e "\e[32;1m || ${subj} --> Final anatomical brainmask already created ('ANAT_brainmask.nii.gz'). \e[0m"
	fi

	## Create precise brainmask and strip T1
	if [[ ! -f "${reg_mask_dir}/${subj}__registration_parcellation.nii.gz" ]]; then
		echo "======================================================================================================="
		echo -e " \e[33;1m ${subj} --> Creating Registration parcellation. \e[0m"
		echo "======================================================================================================="
		echo ""
			${pscript} ${scripts_dir}/create_registration_parcellation_final.py "${subj}"
	else
		echo -e "\e[32;1m || ${subj} --> Registration parcellation already created ('${subj}__registration_parcellation.nii.gz'). \e[0m"
	fi


################## Define registration files

	### Images for Affine
	step0_telencephale_ANAT="${PET_dir}/telencephale_seg.nii.gz"
	step0_telencephale_MNI="${atlas_dir}/MNI152_telencephale_seg.nii.gz"

	MNI_seg_Final="${atlas_dir}/MNI152__registration_parcellation.nii.gz"
	MNI_T1="${atlas_dir}/MNI152_T1_1mm_brain.nii.gz"
	reg_MNI_T1="${atlas_dir}/MNI152_T1_1mm_Registration_brain.nii.gz"

	ANAT_seg_Final="${reg_mask_dir}/${subj}__registration_parcellation.nii.gz"
	ANAT_T1="${PET_dir}/ANAT_T1.nii.gz"

	SUVR_in_ANAT="${PET_dir}/SUVR-in-ANAT.nii.gz"
	SUVR_in_MNI_nLin="${PET_dir}/SUVR-in-MNI_nLin.nii.gz"

	Lin__prefix_Final="${PET_dir}/Lin_telencephale-only__ANAT-to-MNI-1mm__"
	Lin__FWD_transforms_Final="${Lin__prefix_Final}0GenericAffine.mat"
	Lin__ANAT_seg_Final="${Lin__prefix_Final}__Final.nii.gz"
	Lin__ANAT_T1="${Lin__prefix_Final}__ANAT_T1.nii.gz"

################ REGISTRATIONS + PROJECTION

	if [[ ! -f "${Lin__FWD_transforms_Final}" ]]; then
		############# STEP 0:  Affine registration
		echo -e "\n\n\e[33;1m *******  ${subj}  -->  STAGE 0:  LINEAR (Affine) REGISTRATION \e[0m"
		antsRegistration \
			--dimensionality 3 \
			--initial-moving-transform [${step0_telencephale_MNI},${step0_telencephale_ANAT},1] \
			--metric Mattes[${step0_telencephale_MNI},${step0_telencephale_ANAT},1,50,regular,0.2] \
			--transform Affine[0.25] \
			--convergence 2100x1200x1200x10 \
			--smoothing-sigmas 3x2x1x0 \
			--shrink-factors 6x4x2x1 \
			--use-histogram-matching 1 \
			--collapse-output-transforms 1 \
			--output ${Lin__prefix_Final} \
			--masks [NA,NA] \
			--float 1 \
			--write-composite-transform 0 \
			--verbose 1
	fi

	if [[ ! -f "${Lin__ANAT_seg_Final}" ]]; then
		####################### Project ${ANAT_seg_Final}
		echo -e "\n\n ******* \e[33;1m ${subj}  -->   Applying LINEAR transforms  (from STAGE 0) \e[0m"
		antsApplyTransforms \
			--dimensionality 3 \
			--input ${ANAT_seg_Final} \
			--reference-image ${MNI_seg_Final} \
			--transform ${Lin__FWD_transforms_Final} \
			--interpolation NearestNeighbor \
			--output ${Lin__ANAT_seg_Final} \
			--verbose 1
		verify_projection ${Lin__ANAT_seg_Final}
	fi

	if [[ ! -f "${Lin__ANAT_T1}" ]]; then
		##################### Project ${ANAT_T1}
		echo -e "\n\n ******* \e[33;1m ${subj}  -->   Applying LINEAR transforms  (from STAGE 0) \e[0m"
		antsApplyTransforms \
			--dimensionality 3 \
			--input ${ANAT_T1} \
			--reference-image ${reg_MNI_T1} \
			--transform ${Lin__FWD_transforms_Final} \
			--interpolation NearestNeighbor \
			--output ${Lin__ANAT_T1} \
			--verbose 1
		verify_projection ${Lin__ANAT_T1}
	fi


	##############################   STAGE 1
	metric='CC'
	regularization=12.0
	radius=4
	gradient_step=0.9
	suffix="CC-rad${radius}-grad${gradient_step}-rglr${regularization}"
	nLin_prefix_Final="${PET_dir}/nLin_precise-doublemetric-SyN-${suffix}__inMNI-1mm__"
	nLin_FWD_transforms_Final="${nLin_prefix_Final}1Warp.nii.gz"
	nLin_ANAT_seg_Final="${nLin_prefix_Final}__Final.nii.gz"
	nLin__ANAT_T1="${nLin_prefix_Final}__ANAT_T1.nii.gz"

	if [[ ! -f "${nLin_FWD_transforms_Final}" ]]; then
		echo -e "\n\n============================================================================="
		echo -e "\e[33;1m *******  ${subj}  -->  STAGE 1:  NON-LINEAR (SyN) REGISTRATION \e[0m"
		antsRegistration \
			--dimensionality 3 \
			--initial-moving-transform [${MNI_seg_Final},${Lin__ANAT_seg_Final},0] \
			--metric ${metric}[${MNI_seg_Final},${Lin__ANAT_seg_Final},0.5,${radius}] \
			--metric Mattes[${reg_MNI_T1},${Lin__ANAT_T1},0.5,128] \
			--transform SyN[${gradient_step},${regularization},0.0000000001] \
			--convergence [1000x700x400x100x40,1e-10,10] \
			--smoothing-sigmas 4x3x2x1x0 \
			--shrink-factors 7x5x3x2x1 \
			--use-histogram-matching 1 \
			--collapse-output-transforms 1 \
			--output ${nLin_prefix_Final} \
			--masks [NA,NA] \
			--float 1 \
			--write-composite-transform 0 \
			--verbose 1
		verify_projection ${nLin_FWD_transforms_Final}
	fi

	### This projection is purely for QC...
	if [[ ! -f "${nLin_ANAT_seg_Final}" ]]; then
		echo -e "\n ******* \e[33;1m ${subj}  -->  Applying NON-LINEAR transforms  (from STAGE 1: ${suffix}) \e[0m"
		antsApplyTransforms \
			--dimensionality 3 \
			--input ${Lin__ANAT_seg_Final} \
			--reference-image ${MNI_seg_Final} \
			--transform ${nLin_FWD_transforms_Final} \
			--interpolation NearestNeighbor \
			--output ${nLin_ANAT_seg_Final} \
			--verbose 1
		verify_projection ${nLin_ANAT_seg_Final}
	fi

	### This projection is purely for QC...
	if [[ ! -f "${nLin__ANAT_T1}" ]]; then
		echo -e "\n ******* \e[33;1m ${subj}  -->  Applying LINEAR + NON-LINEAR transforms  (from STAGE 0-1: ${suffix}) \e[0m"
		antsApplyTransforms \
			--dimensionality 3 \
			--input ${ANAT_T1} \
			--reference-image ${reg_MNI_T1} \
			--transform ${nLin_FWD_transforms_Final} \
			--transform ${Lin__FWD_transforms_Final} \
			--interpolation Linear \
			--output ${nLin__ANAT_T1} \
			--verbose 1
		verify_projection ${nLin__ANAT_T1}
	fi

	### This is the PET data in MNI (1mm) space.
	if [[ ! -f "${SUVR_in_MNI_nLin}" ]]; then
		echo -e "\n ******* \e[33;1m ${subj}  -->  Applying LINEAR + NON-LINEAR transforms to SUVR PET data \e[0m"
		antsApplyTransforms \
			--dimensionality 3 \
			--input ${SUVR_in_ANAT} \
			--reference-image ${MNI_T1} \
			--transform ${nLin_FWD_transforms_Final} \
			--transform ${Lin__FWD_transforms_Final} \
			--interpolation Linear \
			--output ${SUVR_in_MNI_nLin} \
			--verbose 1
		verify_projection ${SUVR_in_MNI_nLin}
	else
		echo -e "\e[32;1m || ${subj} --> Final SUVR in MNI already projected. \e[0m"
	fi


################################################################################

	SUVR_in_MNI_2mm_nLin="${PET_dir}/nLin_precise_SUVR-in-MNI-2mm.nii.gz"
	#######  Resampling to 2mm
	if [[ ! -f "${SUVR_in_MNI_2mm_nLin}" ]]; then
		echo "============================================================================================================================="
		echo -e " \e[33;1m ${subj} --> Resampling to 2mm (MNI152 template). \e[0m"
		echo "============================================================================================================================="
		echo ""
		mri_vol2vol \
			--mov ${SUVR_in_MNI_nLin} \
			--targ ${atlas_dir}/MNI152_T1_2mm_brain.nii.gz \
			--regheader \
			--cubic \
			--o ${SUVR_in_MNI_2mm_nLin}
	else
		echo -e "\e[32;1m || ${subj} --> Precise SUVRs in MNI were already resampled. \e[0m"
	fi


	############### Extraction segmentations in MNI152
	###################  SUVR <---  precise Registration

	# Associative array defining atlas names and corresponding paths
	declare -A atlases
	atlases=(
		[Glasser_2mm]="${atlas_dir}/glasser_MNI152NLin6Asym_labels_p20_2mm.nii.gz"
		[Tian_S1_2mm]="${atlas_dir}/Tian_Subcortex_S1_3T.nii.gz"
		[Tian_S2_2mm]="${atlas_dir}/Tian_Subcortex_S2_3T.nii.gz"
	)


	# Loop over each atlas using the associative array
	for atlas_name in "${!atlases[@]}"; do
		atlas="${atlases[$atlas_name]}"
		stat_file="${PET_dir}/${subj}__${atlas_name}__precise-SUVR-in-MNI.ROI.stats.dat"
			# Extract ROI values
		if [[ ! -f "${stat_file}" ]]; then
			echo "==============================================================================================================="
			echo -e " \e[33;1m ${subj} --> Extracting ${atlas} ROI segmentations in MNI space (from ${SUVR_in_MNI_2mm_nLin}). \e[0m"
			echo "==============================================================================================================="
			echo ""
				mri_segstats \
					--i "${SUVR_in_MNI_2mm_nLin}" \
					--seg "${atlas}" \
					--sum "${stat_file}"
				cp -v \
					"${stat_file}" \
					"${stats_files_dir}/"
		else
			echo -e "\e[32;1m || ${subj} --> Extraction of ${atlas} ROI values was already done. \e[0m"
		fi
	done


	####### Process Smoothing and masking in MNI  (Precise SUVR in MNI)
	MNI152_2mm_stripped__mask="${atlas_dir}/MNI152_T1_2mm_stripped-3mmborders-no-csf_mask.nii.gz"

	sm=("8")
	smoothed_nLin_precise_SUVR_in_MNI="${subj}__sm${sm}-nLin_precise_SUVR-in-MNI-2mm.nii.gz"

	if [[ ! -f "${PET_dir}/${smoothed_nLin_precise_SUVR_in_MNI}" ]]; then
		echo "=============================================================================================="
		echo -e " \e[33;1m ${subj} --> Smoothing (${sm}) brain PET data in MNI space (2mm). \e[0m"
		echo "=============================================================================================="
		echo ""
		mri_fwhm --smooth-only \
				--i "${SUVR_in_MNI_2mm_nLin}" \
				--fwhm ${sm} \
				--mask "${MNI152_2mm_stripped__mask}" \
				--o "${PET_dir}/${smoothed_nLin_precise_SUVR_in_MNI}"
		cp -v \
			"${PET_dir}/${smoothed_nLin_precise_SUVR_in_MNI}" \
			${PET_data_dir}/
	else
		echo -e "\e[32;1m || ${subj} --> ${smoothed_nLin_precise_SUVR_in_MNI} exists already. \e[0m"
	fi
}


######## Motion correction function
## MC is done in 2 pass:
## First pass to obtain second pass template;
## Second pass is performed on raw PET data.
## Motion correction with lowest final frame-wise displacement is kept as final motion-correction and summed into a static image.

motion_correction_function () {

	local subj="$1"
	FRAMES=("0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14")

	for refvol in "${FRAMES[@]}"; do

			## mc naming
			motion_correction_prefix="${subj}__mc_sm10_ref${refvol}_leastsquares"
			## first-pass files naming
			motion_corrected_PET_1stpass="${motion_correction_prefix}_PET_1stpass.nii.gz"
			# second pass files naming
			template_2ndpass="${motion_correction_prefix}_PET_mean.nii.gz"
			motion_corrected_PET_2ndpass="${motion_correction_prefix}_PET_2ndpass.nii.gz"
			## motion analysis
			motion_outliers_2ndpass="${motion_correction_prefix}_PET_2ndpass_motion_outliers.dat"
			fd_metric_2ndpass="${motion_correction_prefix}_PET_2ndpass_fd_metric.dat"
			fd_plot_2ndpass="${motion_correction_prefix}_PET_2ndpass_fd_plot.png"

			#1#1#1#1#1#1#1#1#1#1#1#1#1#1  FIRST PASS
			echo "========================================================================================================="
			echo -e " \e[33;1m ${subj} --> First-pass motion correction of PET data (to obtain the Second-pass template). \e[0m"
			echo "========================================================================================================="
			echo ""
			## motion correction of PET data
			if [[ ! -f "${mc_dir}/${motion_corrected_PET_1stpass}" ]]; then
				mcflirt -in ${PET_dir}/PET.nii -smooth 10 -cost leastsquares \
						-refvol ${refvol} \
						-out ${mc_dir}/${motion_corrected_PET_1stpass} \
						-plots -stats -report
			fi

			echo "======================================================================================================="
			echo -e " \e[33;1m ${subj} --> Averaging First-pass motion corrected frames to obtain Second-pass template. \e[0m"
			echo "======================================================================================================="
			echo ""
			## average first pass mc data
			if [[ ! -f "${mc_dir}/${template_2ndpass}" ]]; then
				mri_concat ${mc_dir}/${motion_corrected_PET_1stpass} --mean --o ${mc_dir}/${template_2ndpass}
			fi


			##2##2##2##2##2##2##2##2##2 SECOND PASS
			echo "==============================================================================================================="
			echo -e " \e[33;1m ${subj} --> Second-pass motion correction of PET data (to obtain final motion corrected frames). \e[0m"
			echo "==============================================================================================================="
			echo ""
			## motion correction of PET data
			if [[ ! -f "${mc_dir}/${motion_corrected_PET_2ndpass}" ]]; then
				mcflirt -in ${PET_dir}/PET.nii -smooth 10 -cost leastsquares \
								-reffile ${mc_dir}/${template_2ndpass} \
								-out ${mc_dir}/${motion_corrected_PET_2ndpass} \
								-plots -stats -report
			fi

			echo "======================================================================================================================="
			echo -e " \e[33;1m ${subj} --> Residual motion analysis of Second-pass (final) motion correction (frame-wise displacements). \e[0m"
			echo "======================================================================================================================="
			echo ""
			if [[ ! -f "${mc_dir}/${fd_metric_2ndpass}" ]]; then
			## motion analysis
			fsl_motion_outliers -i ${mc_dir}/${motion_corrected_PET_2ndpass} \
													-o ${mc_dir}/${motion_outliers_2ndpass} \
													-s ${mc_dir}/${fd_metric_2ndpass} \
													-p ${mc_dir}/${fd_plot_2ndpass} \
													--fd -v
			fi
	done

	## initialize variables
	best_fd="100"
	best_refvol=""
	## log the best frame for motion correction
	mc_frame_ref_file="${mc_dir}/mc_frame_ref.dat"
	mc_error_log="${mc_dir}/mc_errors.log"

	for refvol in "${FRAMES[@]}"; do

			fd_metric_file="${mc_dir}/${subj}__mc_sm10_ref${refvol}_leastsquares_PET_2ndpass_fd_metric.dat"

			if [[ -f "${fd_metric_file}" ]]; then
					# Compute average FD from fd metric file. Skips first row and Assumes one FD value per row.
					avg_fd=$(awk 'NR > 1 { sum += $1; count++ } END { if (count>0) print sum/count; else print 100 }' "${fd_metric_file}")

					# Compare using bc for floating-point comparison
					cmp=$(echo "$avg_fd < $best_fd" | bc -l)
					if [[ ${cmp} -eq 1 ]]; then
							best_fd=${avg_fd}
							best_refvol=${refvol}
					fi
			else
					echo "[ $(date '+%Y-%m-%d  %H:%M') ]  -->  Warning  :::  FD metric file not found for refvol/frame ${refvol}" >> "${mc_error_log}"
			fi
	done

	if [[ -n "${best_refvol}" ]]; then
			echo "${best_refvol}" > "${mc_frame_ref_file}"
			echo " *******  ${subj}:: Best reference volume is ${best_refvol} (average FD = ${best_fd})"
			echo "Saving best motion corrected PET data to PET processing directory with final name."
					## Save with final name in ${PET_dir}
					cp -v ${mc_dir}/${subj}__mc_sm10_ref${best_refvol}_leastsquares_PET_2ndpass.nii.gz \
							${PET_dir}/mc_PET.nii.gz
	else
			echo "Error: please check ${mc_error_log} !" > "${mc_frame_ref_file}"
			echo "Error:  ${subj}  -->  No FD metric files found to determine best reference volume."
				return
	fi

}



######## Function to create the REFERENCE REGION MASK
## Uses deep white matter (non-lobar) based on MRI anatomical T1, removes voxels too close to striatum (splill-over), and removes small isolated clusters
create_wm_ref_region_function () {

	local subj="$1"
	# convert wmparc
	mri_convert ${mri_dir}/wmparc.mgz ${ref_dir}/wmparc.nii.gz
	# Binarize region 5001-5002:  Undefined white matter (non-lobar)
	fslmaths ${ref_dir}/wmparc.nii.gz -thr 5001 -uthr 5002 -bin ${ref_dir}/non-lobar_wm_mask.nii.gz

	# convert segmentation
	mri_convert ${mri_dir}/aparc+aseg.mgz ${ref_dir}/aparc+aseg.nii.gz
	# Binarize region 11: Left-Caudate
	fslmaths ${ref_dir}/aparc+aseg.nii.gz -thr 11 -uthr 11 -bin ${ref_dir}/seg11.nii.gz
			# dilate
			fslmaths ${ref_dir}/seg11.nii.gz -dilM -dilM -dilM -dilM ${ref_dir}/seg11+4vx.nii.gz
	# Binarize region 50: Right-Caudate
	fslmaths ${ref_dir}/aparc+aseg.nii.gz -thr 50 -uthr 50 -bin ${ref_dir}/seg50.nii.gz
			# dilate
			fslmaths ${ref_dir}/seg50.nii.gz -dilM -dilM -dilM -dilM ${ref_dir}/seg50+4vx.nii.gz
	# Binarize region 12: Left-Putamen
	fslmaths ${ref_dir}/aparc+aseg.nii.gz -thr 12 -uthr 12 -bin ${ref_dir}/seg12.nii.gz
			# dilate
			fslmaths ${ref_dir}/seg12.nii.gz -dilM -dilM -dilM -dilM ${ref_dir}/seg12+4vx.nii.gz
	# Binarize region 51: Right-Putamen
	fslmaths ${ref_dir}/aparc+aseg.nii.gz -thr 51 -uthr 51 -bin ${ref_dir}/seg51.nii.gz
			# dilate
			fslmaths ${ref_dir}/seg51.nii.gz -dilM -dilM -dilM -dilM ${ref_dir}/seg51+4vx.nii.gz
	# Binarize region 26: Left-Accumbens-area
	fslmaths ${ref_dir}/aparc+aseg.nii.gz -thr 26 -uthr 26 -bin ${ref_dir}/seg26.nii.gz
			# dilate
			fslmaths ${ref_dir}/seg26.nii.gz -dilM -dilM -dilM -dilM ${ref_dir}/seg26+4vx.nii.gz
	# Binarize region 58: Right-Accumbens-area
	fslmaths ${ref_dir}/aparc+aseg.nii.gz -thr 58 -uthr 58 -bin ${ref_dir}/seg58.nii.gz
			# dilate
			fslmaths ${ref_dir}/seg58.nii.gz -dilM -dilM -dilM -dilM ${ref_dir}/seg58+4vx.nii.gz
	# Binarize region 10: Left-Thalamus-Proper
	fslmaths ${ref_dir}/aparc+aseg.nii.gz -thr 10 -uthr 10 -bin ${ref_dir}/seg10.nii.gz
			# dilate
			fslmaths ${ref_dir}/seg10.nii.gz -dilM -dilM -dilM -dilM ${ref_dir}/seg10+4vx.nii.gz
	# Binarize region 49: Right-Thalamus-Proper
	fslmaths ${ref_dir}/aparc+aseg.nii.gz -thr 49 -uthr 49 -bin ${ref_dir}/seg49.nii.gz
			# dilate
			fslmaths ${ref_dir}/seg49.nii.gz -dilM -dilM -dilM -dilM ${ref_dir}/seg49+4vx.nii.gz
	# Binarize region 16: Brainstem
	fslmaths ${ref_dir}/aparc+aseg.nii.gz -thr 16 -uthr 16 -bin ${ref_dir}/seg16.nii.gz
			# dilate
			fslmaths ${ref_dir}/seg16.nii.gz -dilM -dilM -dilM -dilM ${ref_dir}/seg16+4vx.nii.gz
	# Binarize region 28: left ventral DC
	fslmaths ${ref_dir}/aparc+aseg.nii.gz -thr 28 -uthr 28 -bin ${ref_dir}/seg28.nii.gz
			# dilate
			fslmaths ${ref_dir}/seg28.nii.gz -dilM -dilM -dilM -dilM ${ref_dir}/seg28+4vx.nii.gz
	# Binarize region 60: right ventral DC
	fslmaths ${ref_dir}/aparc+aseg.nii.gz -thr 60 -uthr 60 -bin ${ref_dir}/seg60.nii.gz
			# dilate
			fslmaths ${ref_dir}/seg60.nii.gz -dilM -dilM -dilM -dilM ${ref_dir}/seg60+4vx.nii.gz
	# Combine and binarize the striatal regions
	fslmaths ${ref_dir}/seg11+4vx.nii.gz \
			-add ${ref_dir}/seg50+4vx.nii.gz \
			-add ${ref_dir}/seg12+4vx.nii.gz \
			-add ${ref_dir}/seg51+4vx.nii.gz \
			-add ${ref_dir}/seg26+4vx.nii.gz \
			-add ${ref_dir}/seg58+4vx.nii.gz \
			-add ${ref_dir}/seg10+4vx.nii.gz \
			-add ${ref_dir}/seg49+4vx.nii.gz \
			-add ${ref_dir}/seg16+4vx.nii.gz \
			-add ${ref_dir}/seg28+4vx.nii.gz \
			-add ${ref_dir}/seg60+4vx.nii.gz \
			-bin ${ref_dir}/binary-mask_dilated_subcortex-brainstem.nii.gz

	# Invert binary mask
	fslmaths ${ref_dir}/binary-mask_dilated_subcortex-brainstem.nii.gz -binv ${ref_dir}/inverted_binary-mask_dilated_subcortex-brainstem.nii.gz

	# Remove dilated striatal voxels from wm ref region
	fslmaths ${ref_dir}/non-lobar_wm_mask.nii.gz \
			-mul ${ref_dir}/inverted_binary-mask_dilated_subcortex-brainstem.nii.gz \
			${ref_dir}/${subj}_wm-ref-region-mask_excl-dilated-subcortex-brainstem.nii.gz

	# remove small clusters (<500 voxels)
	${pscript} ${scripts_dir}/clean_wm-ref-region.py "${subj}"

}
################################################################################


################################################################################
## This function was used to confirm creation
verify_projection () {

	local projection="$1"
	if [[ -f ${projection} ]]; then
		echo -e "\e[31;1m Projection:\n ${projection}\e[0m \n\n"
	fi
}
################################################################################


######################### brainmask without external CSF, created using FreeSurfer recon-all denoised T1, stripped with mri_synthstrip, and adding aparc+aseg voxels to avoid overstripping.
creating_brainmask_and_stripping__function () {

	### Convert antsdn.brain.nii.gz
	if [[ ! -f "${PET_dir}/antsdn.brain.nii.gz" ]]; then
		echo "================================================================================================="
		echo -e " \e[33;1m ${subj} --> Preparing denoised T1 to create anatommical brainmask and stripped T1. \e[0m"
		echo "================================================================================================="
		echo ""
			## Preparing denoised T1
			mri_convert ${mri_dir}/antsdn.brain.mgz ${PET_dir}/antsdn.brain.nii.gz
	else
		echo -e "\e[32;1m || ${subj} --> Preparation of denoised T1 already done (antsdn.brain.nii.gz). \e[0m"
	fi

	### mri_synthstrip
	if [[ ! -f "${PET_dir}/synthstripped-brainmask_no-csf.nii.gz" ]]; then
		echo "=============================================================================================="
		echo -e " \e[33;1m ${subj} --> Creating stripped brainmask. \e[0m"
		echo "=============================================================================================="
		echo ""
			## synthstrip skull and CSF from denoised T1 to get a brainmask
			mri_synthstrip -i ${PET_dir}/antsdn.brain.nii.gz \
					-m ${PET_dir}/synthstripped-brainmask_no-csf.nii.gz \
					--no-csf
	else
		echo -e "\e[32;1m || ${subj} --> Creation of stripped brainmask was already completed. \e[0m"
	fi

	### Fill holes
	if [[ ! -f "${PET_dir}/synthstripped-brainmask_no-external-csf.nii.gz" ]]; then
		echo "=============================================================================================="
		echo -e " \e[33;1m ${subj} --> Filling brainmask ventricles. \e[0m"
		echo "=============================================================================================="
		echo ""
			fslmaths ${PET_dir}/synthstripped-brainmask_no-csf.nii.gz -fillh ${PET_dir}/synthstripped-brainmask_no-external-csf.nii.gz
	else
		echo -e "\e[32;1m || ${subj} --> Filling of stripped anatomical brainmask (ventricles) was already completed. \e[0m"
	fi

	### Include aparc+aseg.nii.gz in brainmask
	if [[ ! -f "${PET_dir}/ANAT_brainmask.nii.gz" ]]; then
		echo "======================================================================================================="
		echo -e " \e[33;1m ${subj} --> Creating final anatomical (T1) brainmask by adding brain recon-all segmentations. \e[0m"
		echo "======================================================================================================="
		echo ""
			fslmaths ${PET_dir}/synthstripped-brainmask_no-external-csf.nii.gz \
				-add ${ref_dir}/aparc+aseg.nii.gz \
				-fillh \
				-bin ${PET_dir}/ANAT_brainmask.nii.gz
	else
		echo -e "\e[32;1m || ${subj} --> Final anatomical brainmask already created ('ANAT_brainmask.nii.gz'). \e[0m"
	fi

	###
	if [[ ! -f "${PET_dir}/ANAT_T1.nii.gz" ]]; then
		echo "======================================================================================================="
		echo -e " \e[33;1m ${subj} --> Creating final anatomical (T1) brainmask by brain recon-all segmentations. \e[0m"
		echo "======================================================================================================="
		echo ""
			fslmaths ${PET_dir}/antsdn.brain.nii.gz \
				-mul ${PET_dir}/ANAT_brainmask.nii.gz \
				${PET_dir}/ANAT_T1.nii.gz
	else
		echo -e "\e[32;1m || ${subj} --> Anatomical T1 already stripped ('ANAT_T1.nii.gz'). \e[0m"
	fi
}

#########################
#########################


#####  Get subject list from /subjects/ directories
## FEOBV
subjects=($(find "${project_dir}/subjects" -maxdepth 1 -name 'FEOBV*' -exec basename {} \; | sort))

#####  Declare subject list overtly.
# declare -a subjects=("FEOBV002")

echo ""
echo ""
echo "#########################################################"
echo "###################   SUBJECTS LIST   ###################"
echo "#########################################################"
echo -e "\e[32;1m"
subj_cols=5
subj_count=0
for subj in "${subjects[@]}"; do
	printf "%-10s" "${subj}"
	((subj_count++))
	if ((subj_count % subj_cols == 0)); then
		echo ""
	fi
done
## limit of processes to reach
process_limit=10
## set process counter to 0
process_count=0

echo -e "\e[0m"
echo ""
echo "#############################################################"
echo "################     PARALLEL PROCESSING     ################"
echo "################   Concurrent processes: ${process_limit}  ################"
echo "#############################################################"
echo ""
echo ""

## wait momentarily to glance at the subjects list
sleep 1.5
## start time
start_time=$(date +%s)

## loop through subjects...
for subj in "${subjects[@]}"; do
	sleep 0.07
	analyze_pbet_function ${subj} &

	# Increment the counter
	((process_count++))

	# If we've reached <num> parallel processes, wait for them to finish
	if [ ${process_count} -eq ${process_limit} ]; then
		echo ""
		echo -e "\e[31;1m *!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*"
		echo " *!*!*  Reached limit of parallel processes (${process_limit}). Waiting for them to complete.  *!*!*!*"
		echo -e " *!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!* \e[0m"
		echo ""

		wait
		process_count=0
		echo ""
		echo -e "\e[31;1m *!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*"
		echo "*!*!*!*!*!*!*!*!* Completed batch. Resetting counter. *!*!*!*!*!*!*!*!*!*"
		echo -e "*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!* \e[0m"
		echo ""
	fi
done

# Wait for any remaining processes to finish
wait
## time elapsed
end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
elapsed_formatted=$(printf '%02dh:%02dm:%02ds\n' $((elapsed/3600)) $(( (elapsed%3600)/60 )) $((elapsed%60)))
echo ""
echo "#########################################################"
echo -e "############# \e[32;1m  Time Elapsed: ${elapsed_formatted} \e[0m  #############"
echo -e "############# * \e[32;1m All subjects processed. \e[0m * #############"
echo "#########################################################"
