import os
import sys
import glob
import time
import nibabel as nib
import numpy as np
import pandas as pd
from scipy.ndimage import distance_transform_edt
import argparse

## Parse subject ID argument
parser = argparse.ArgumentParser(description='Create registration parcellation.')
parser.add_argument('subj', type=str, help='Subject ID')
args = parser.parse_args()
subj = args.subj

### Directories
## Retrieve the SUBJECTS_DIR environment variable
FS_SUBJECTS_DIR = os.environ.get('SUBJECTS_DIR')
## Define project_dir
project_dir = os.path.dirname(FS_SUBJECTS_DIR)
## Define subject's PET_dir
PET_dir = f'{project_dir}/subjects/{subj}/pet/001'
## Define reference region dir
ref_dir = f'{PET_dir}/reference_region'
## Define registration mask dir
reg_mask_dir = f'{PET_dir}/registration_mask'

### File paths
aparc_aseg_path = f'{reg_mask_dir}/aparc.a2009s+aseg.nii.gz'
synthseg_path = f'{reg_mask_dir}/synthseg_anat.nii.gz'
wmparc_path = f'{ref_dir}/wmparc.nii.gz'
registration_parcellation_output_path = f'{reg_mask_dir}/{subj}__registration_parcellation.nii.gz'
telencephale_out_path = f"{PET_dir}/telencephale_seg.nii.gz"

print(f' ***  Making parcellation image for registration --> {subj}')
## Label Definitions
ctx_rh_IDs = list(range(12101, 12176))
ctx_lh_IDs = list(range(11101, 11176))
ctx_IDs = ctx_rh_IDs + ctx_lh_IDs
subctx_IDs = [10, 11, 12, 13, 17, 18, 26, 49, 50, 51, 52, 53, 54, 58]
subctx_l_IDs = [10, 11, 12, 13, 17, 18, 26]
subctx_r_IDs = [49, 50, 51, 52, 53, 54, 58]
l_cereb = [7, 46]
r_cereb = [8, 47]
mid_and_stem = [16, 28, 60]
csf_IDs = [4, 5, 14, 15, 43, 44, 24]  # Not remapped
wm_lh_IDs = [
	5001, 3001, 3002, 3003, 3005, 3006, 3007, 3008, 3009, 3010, 3011, 3012, 3013, 3014, 3015, 3016, 3017,
	3018, 3019, 3020, 3021, 3022, 3023, 3024, 3025, 3026, 3027, 3028, 3029, 3030, 3031, 3032, 3033,
	3034, 3035]
wm_rh_IDs = [
	5002, 4001, 4002, 4003, 4005, 4006, 4007, 4008, 4009, 4010, 4011, 4012, 4013, 4014, 4015,
	4016, 4017, 4018, 4019, 4020, 4021, 4022, 4023, 4024, 4025, 4026, 4027, 4028, 4029, 4030, 4031,
	4032, 4033, 4034, 4035]
corpuscallosum_IDs = [251, 252, 253, 254, 255]
all_wm_IDs = [
	5001, 3001, 3002, 3003, 3005, 3006, 3007, 3008, 3009, 3010, 3011, 3012, 3013, 3014, 3015, 3016, 3017,
	3018, 3019, 3020, 3021, 3022, 3023, 3024, 3025, 3026, 3027, 3028, 3029, 3030, 3031, 3032, 3033,
	3034, 3035,
	5002, 4001, 4002, 4003, 4005, 4006, 4007, 4008, 4009, 4010, 4011, 4012, 4013, 4014, 4015,
	4016, 4017, 4018, 4019, 4020, 4021, 4022, 4023, 4024, 4025, 4026, 4027, 4028, 4029, 4030, 4031,
	4032, 4033, 4034, 4035,
	251, 252, 253, 254, 255]
								####### now make masks... csf and ventricle excluded?
## Label Remapping Function
def spaced_labels(start, end, original_labels):
	interval = (end - start) // len(original_labels)
	return {old: start + i * interval for i, old in enumerate(original_labels)}

## Load Images
seg_img = nib.load(aparc_aseg_path)
seg_data = seg_img.get_fdata().astype(int)

synthseg_img = nib.load(synthseg_path)
synthseg_data = synthseg_img.get_fdata().astype(int)

wm_img = nib.load(wmparc_path)
wm_data = wm_img.get_fdata().astype(int)

## Start with synthseg as base
base_data = synthseg_data.copy()

## Replace cortical voxels in synthseg with values from aparc+aseg
aparc_ctx_mask = np.isin(seg_data, ctx_IDs)
synth_cortex_mask = (synthseg_data == 3) | (synthseg_data == 42)
replacement_mask = synth_cortex_mask & aparc_ctx_mask
base_data[replacement_mask] = seg_data[replacement_mask]

## Fill remaining synthseg cortical voxels using nearest aparc+aseg cortex label
ctx_gap_mask = synth_cortex_mask & (~aparc_ctx_mask)
distance, indices = distance_transform_edt(~aparc_ctx_mask, return_indices=True)
ctx_gap_indices = np.where(ctx_gap_mask)
nearest_labels = seg_data[tuple(idx[ctx_gap_indices] for idx in indices)]
base_data[ctx_gap_indices] = nearest_labels

## Overwrite white matter voxels from wmparc
wm_mask = np.isin(wm_data, all_wm_IDs)
base_data[wm_mask] = wm_data[wm_mask]
## Define synthseg white matter mask
synthseg_wm_mask = (synthseg_data == 2) | (synthseg_data == 41)
## Identify white matter voxels in synthseg that did NOT get a label from wmparc
wm_gap_mask = synthseg_wm_mask & (~wm_mask)
## Distance transform from valid wmparc labels
distance_wm, indices_wm = distance_transform_edt(~wm_mask, return_indices=True)
## Find gap voxel coordinates
wm_gap_indices = np.where(wm_gap_mask)
## Get nearest labels from wmparc for those voxels
nearest_wm_labels = wm_data[tuple(idx[wm_gap_indices] for idx in indices_wm)]
## Assign them to base_data
base_data[wm_gap_indices] = nearest_wm_labels

## Remap Labels
remap_dict = {
	**spaced_labels(14000, 17000, ctx_rh_IDs),
	**spaced_labels(14000, 17000, ctx_lh_IDs),
	**spaced_labels(7000, 10000, subctx_r_IDs),
	**spaced_labels(7000, 10000, subctx_l_IDs),
	**spaced_labels(2000, 3000, wm_rh_IDs),
	**spaced_labels(2000, 3000, wm_lh_IDs),
	**spaced_labels(3700, 3710, corpuscallosum_IDs),
	**spaced_labels(500, 700, r_cereb),
	**spaced_labels(500, 700, l_cereb),
	**spaced_labels(1100, 1400, mid_and_stem)
}

remapped_data = np.copy(base_data)
for old_val, new_val in remap_dict.items():
	remapped_data[base_data == old_val] = new_val

## Save Final Parcellation
final_img = nib.Nifti1Image(remapped_data.astype(np.int32), seg_img.affine, seg_img.header)
nib.save(final_img, registration_parcellation_output_path)
## Verification of existence
if os.path.exists(registration_parcellation_output_path):
	print()
	print(f"Final registration parcellation saved.  Download with:")
	print(f"BICget {registration_parcellation_output_path}")
else:
	print()
	print(f" %%%%% ERROR :: Failed to save final parcellation to:")
	print(f" {registration_parcellation_output_path}")
###############################################################################################


###################################
###############  telencephale only without CSF: this is used for intial linear registration to MNI template. Effectively, this favors cortical alignment.
## Create a copy of the remapped data
telencephale_data = np.copy(remapped_data)
## Set CSF to 0
telencephale_data[np.isin(telencephale_data, csf_IDs)] = 0
### Set excluded regions to 0
not_in_telencephale_IDs = mid_and_stem + r_cereb + l_cereb
telencephale_data[np.isin(telencephale_data, not_in_telencephale_IDs)] = 0
### Save
telencephale_img = nib.Nifti1Image(telencephale_data, seg_img.affine, seg_img.header)
nib.save(telencephale_img, telencephale_out_path)
## Verification of existence
if os.path.exists(telencephale_out_path):
	print()
	print(f"Telencephale segmentation for registration saved.  Download with:")
	print(f"BICget {telencephale_out_path}")
else:
	print()
	print(f" %%%%% ERROR :: Failed to create Telencephale parcellation... %%%%%")
