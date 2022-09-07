"""
Get synapse info: presynaptic cell and size
"""

import argparse
from sys import argv

import numpy as np
import pandas as pd
import tifffile as tif

from cloudvolume import CloudVolume


## Load volume
def load_volume(path, res):
	"""
	INPUT:
	path : Path of the volume
	res : (x, y, z) resolution in nm

	OUTPUT:
	vol : Image volume
	"""

	if path[:2] == "gs":

		cloud_vol = CloudVolume(path, mip=res, 
			parallel=True, progress=False)
		vol = cloud_vol[:,:,:][:,:,:,0]

	else:

		vol = tif.imread(path)

	return vol


if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument("--syn_seg", required=True, type=str,
		help="Synapse segmentation volume")
	parser.add_argument("--cell_seg", required=True, type=str,
		help="Cell segmentation volume")
	parser.ad_argument("--res", required=True, nargs=3, type=int,
		help="Resolution in nm")
	parser.add_argument("--outpath", required=True, type=str,
		default="synapse_info.csv",
		help="Output synapse info file path")

	opt = parser.parse_args()

	syn_seg_path = opt.syn_seg
	seg_path = opt.cell_seg

	res = opt.res

	# Load volumes
	syn_seg = load_volume(syn_seg_path, res)
	seg = load_volume(seg_path, res)

	pre_id_dict = {}
	size_dict = {}

	nsect = syn_seg.shape
	for i in range(nsect):

		syn_seg_vec = syn_seg[:,:,i].reshape((-1,))
		seg_vec = seg[:,:,i].reshape((-1,))

		syn_list_sect, idx_sect, syn_size_sect = np.unique(syn_seg_vec,
																											return_index=True,
																											return_counts=True)
		
		# Remove segment id 0
		syn_list_sect = syn_list_sect[1:]
		idx_sect = idx_sect[1:]
		syn_size_sect = syn_size_sect[1:]

		nlist = syn_list_sect.shape[0]
		for j in range(nlist):

			syn_id = syn_list_sect[j]
			pre_id_dict[syn_id] = seg_vec[idx_sect[j]]
			size_dict[syn_id] += syn_size_sect[j]*np.prod(res)

		if (i+1)<=10 or (i+1)%50==0:
			print(f"{i+1} / {nsect}")

	syn_id_list = list(pre_id_dict.keys())
	pre_id_list = list(pre_id_dict.values())
	syn_size_list = list(size_dict.values())

	syn_info = {"syn_id": syn_id_list,
							"pre_id": pre_id_list,
							"size": syn_size_list}

	out_df = pd.DataFrame(data=syn_info)
	out_df.to_csv(opt.outpath, indeex=False)
