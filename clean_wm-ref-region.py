import numpy as np
import nibabel as nib
import sys
import os

## Define the 6 possible directions for neighboring voxels in a 3D matrix (adjacent voxels)
directions = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
## set minimum cluster size in voxels
min_cluster_size = 500

def in_bounds(x, y, z, matrix):
	""" Check if the given coordinates are within the bounds of the matrix. """
	return 0 <= x < matrix.shape[0] and 0 <= y < matrix.shape[1] and 0 <= z < matrix.shape[2]

def dfs(matrix, visited, x, y, z):
	""" Perform DFS to find all connected 1s and their size. """
	stack = [(x, y, z)]
	cluster = []
	while stack:
		cx, cy, cz = stack.pop()
		if not in_bounds(cx, cy, cz, matrix) or visited[cx, cy, cz] or matrix[cx, cy, cz] == 0:
			continue
		visited[cx, cy, cz] = True
		cluster.append((cx, cy, cz))
		## Check all 6 possible directions (neighboring voxels in the 3D matrix)
		for dx, dy, dz in directions:
			nx, ny, nz = cx + dx, cy + dy, cz + dz
			stack.append((nx, ny, nz))
	return cluster

def remove_small_clusters(matrix, min_size):
	visited = np.zeros(matrix.shape, dtype=bool)
	for x in range(matrix.shape[0]):
		for y in range(matrix.shape[1]):
			for z in range(matrix.shape[2]):
				if matrix[x, y, z] == 1 and not visited[x, y, z]:
					## Find the cluster using DFS
					cluster = dfs(matrix, visited, x, y, z)
					## Remove the cluster if it's smaller than the minimum size
					if len(cluster) < min_size:
						for cx, cy, cz in cluster:
							matrix[cx, cy, cz] = 0
	return matrix

def main(subj):
	### project dir
	## Retrieve the SUBJECTS_DIR environment variable
	FS_SUBJECTS_DIR = os.environ.get('SUBJECTS_DIR')
	## Set the parent directory of SUBJECTS_DIR as project_dir
	project_dir = os.path.dirname(FS_SUBJECTS_DIR)
	## set subj PET dir
	PET_dir = f'{project_dir}/subjects/{subj}/pet/001'
	## Define reference region dir
	ref_dir = f'{PET_dir}/reference_region'

	input_filepath = f'{ref_dir}/{subj}_wm-ref-region-mask_excl-dilated-subcortex-brainstem.nii.gz'
	output_filepath = f'{PET_dir}/{subj}_clean_wm-ref-region-mask.nii.gz'

	# Load the nifti file
	img = nib.load(input_filepath)
	data = img.get_fdata()

	## Binarize
	data = (data > 0).astype(int)

	## Remove small clusters
	cleaned_data = remove_small_clusters(data, min_cluster_size)

	## Save new reference region mask
	cleaned_img = nib.Nifti1Image(cleaned_data, img.affine, img.header)
	nib.save(cleaned_img, output_filepath)
	print(f"Saved cleaned brain mask to {output_filepath}")

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print("Usage: python clean_wm-ref-region.py <subj>")
	else:
		subj = sys.argv[1]
		main(subj)
