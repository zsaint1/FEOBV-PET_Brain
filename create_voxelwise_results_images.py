import os
import nibabel as nib
from nilearn import plotting
import matplotlib.pyplot as plt
from PIL import Image
from multiprocessing import Pool
from itertools import product
import numpy as np

## directories
project_dir = '/BICNAS2/tuominen/FEOBV-pet'
GLMfits_dir = os.path.join(project_dir, 'PET_group_analyses')
atlas_dir = os.path.join(project_dir, 'atlas-for-ROIs')
ROIs_dir = os.path.join(atlas_dir, 'ROI_masks')
results_dir = os.path.join(project_dir, 'results/voxelwise')
os.makedirs(results_dir, exist_ok=True)

## MNI & ROIs
MNI152_img = nib.load(os.path.join(atlas_dir, 'MNI152_T1_1mm_brain_stripped-for-voxelwise-results.nii.gz'))
BA9_contour_mask = os.path.join(ROIs_dir, 'BA9-glasser_2mm-ROI_mask.nii.gz')
hippocampus_contour_mask = os.path.join(ROIs_dir, 'hippocampus-Tian_in-MNI_2mm_ROI_mask.nii.gz')
BA46_contour_mask	= os.path.join(ROIs_dir, 'BA46.nii.gz')

##################################################
permutations_sig_masked = f"perm.th13.abs.sig.masked.nii.gz"

voxelwise_GLMs = {
	'Patients – Controls': {
		'folder': f'GLMfit__18patients-33controls--finality-nLin_precise_SUVR-in-MNI-2mm__contrast_g2v0_sm8_output___permutations_10000-abs--vwthr1.3/contrast_g2v0',
		'label': 'Group_diff',
		'vmin': -5,
		'vmax': -0.00001,
		'cmap': 'Blues_r',
		'ticks': [-5.0, -2.5, 0.0]
	},
	'CANTAB': {
		'folder': f'GLMfit__18patients-CANTAB--finality-nLin_precise_SUVR-in-MNI-2mm__contrast_g1v1_slope_sm8_output___permutations_10000-abs--vwthr1.3/contrast_g1v1_slope',
		'label': 'CANTAB',
		'vmin': 0.00001,
		'vmax': 5,
		'cmap': 'Reds',
		'ticks': [0, 2.5, 5.0]
	}
}

for info in voxelwise_GLMs.values():
	os.makedirs(os.path.join(results_dir, info['label']), exist_ok=True)


## Final Coordinate list for slices
slices_list = {
	'x': [-38, -24, -4, 4, 24, 38],
	'z': [-38, -22, -14, 26, 32, 46]
}


## Define processing function
def process_slice(args):
	axis, coord = args

	for title, info in voxelwise_GLMs.items():
		simple_slice_path = f"{results_dir}/{info['label']}/{info['label']}--{axis}_{coord}__no-cbar-no-title.png"

		glm_dir = os.path.join(GLMfits_dir, info['folder'])
		sig_path = f"{glm_dir}/{permutations_sig_masked}"
		if not os.path.isfile(sig_path):
			print(f"ERROR: Missing {sig_path}")
			continue

		## save a copy without cbar
		try:
			img = nib.load(sig_path)
			display = plotting.plot_stat_map(
				img,
				bg_img=MNI152_img,
				display_mode=axis,
				cut_coords=[coord],
				title=None,
				black_bg=False,
				alpha=0.9,
				vmax=info['vmax'],
				vmin=info['vmin'],
				cmap=info['cmap'],
				colorbar=False
			)

			display.savefig(simple_slice_path, dpi=700)
			display.close()
			plt.close('all')
			print(f"* Saved: {simple_slice_path}")

		except Exception as e:
			print(f"ERROR plotting {axis}={coord} for {title}: {e}")
			continue


## Parallellize processes
if __name__ == "__main__":
	all_coords = list(product(['x'], slices_list['x'])) + list(product(['z'], slices_list['z']))

	## randomize to balance load processing
	np.random.shuffle(all_coords)

	## adjust to available cores
	process_num = 10

	print(f"Launching {len(all_coords)} slice jobs with {process_num} cores...")
	with Pool(processes=process_num) as pool:
		pool.map(process_slice, all_coords)
