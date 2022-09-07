"""
Get synapse info: presynaptic cell and size
"""

import argparse
from sys import argv

import numpy as np
import pandas as pd
import tifffile as tif

from cloudvolume import CloudVolume

from utils import load_volume

if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument("--syn_seg", required=True, type=str,
		help="Synapse segmentation volume")
	parser.add_argument("--cell_seg", required=True, type=str,
		help="Cell segmentation volume")
	parser.add_argument("--res", required=True, nargs=3, type=int,
		help="Resolution in nm")
	parser.add_argument("--outpath", type=str,
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
	pos_dict = {}

	nsect = syn_seg.shape[2]
	for i in range(nsect):

		syn_seg_sect = syn_seg[:,:,i]
		syn_seg_vec = syn_seg_sect.reshape((-1,))
		seg_sect = seg[:,:,i]
		seg_vec = seg_sect.reshape((-1,))

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

			syn_mask = syn_seg_sect==syn_id

			pre_cand, vx_count = np.unique(seg_sect[syn_mask], return_counts=True)
			if pre_cand[0] == 0 and pre_cand.shape[0] > 1:
				pre_cand = pre_cand[1:]
				vx_count = vx_count[1:]

			pre_id_dict[syn_id] = pre_cand[np.argmax(vx_count)]

			x, y = np.where(syn_mask)
			syn_coord = np.c_[x,y]
			syn_cent = np.mean(syn_coord, axis=0)
			min_idx = np.argmin(np.sum((syn_coord - syn_cent)**2)**0.5)
			syn_loc = np.array(list(syn_coord[min_idx,:]) + [i])

			keys = np.array(list(size_dict.keys()))
			if np.sum(keys==syn_id):
				size_dict[syn_id] += syn_size_sect[j]*np.prod(res)
				pos_dict[syn_id].append(syn_loc)
			else:
				size_dict[syn_id] = syn_size_sect[j]*np.prod(res)
				pos_dict[syn_id] = [syn_loc]

		if (i+1)<=10 or (i+1)%50==0:
			print(f"{i+1} / {nsect}")

	syn_id_list = list(pre_id_dict.keys())
	pre_id_list = list(pre_id_dict.values())
	syn_size_list = list(size_dict.values())

	nsyn = len(syn_id_list)
	syn_posx = []
	syn_posy = []
	syn_posz = []
	for i in range(nsyn):

		syn_pos_list = pos_dict[syn_id_list[i]]
		idx = len(syn_pos_list)//2
		syn_posx.append(syn_pos_list[idx][0]*res[0])
		syn_posy.append(syn_pos_list[idx][1]*res[1])
		syn_posz.append(syn_pos_list[idx][2]*res[2])

	syn_info = {"syn_id": syn_id_list,
							"pre_id": pre_id_list,
							"size": syn_size_list,
							"x_pos": syn_posx,
							"y_pos": syn_posy,
							"z_pos": syn_posz}

	out_df = pd.DataFrame(data=syn_info)
	out_df.to_csv(opt.outpath, index=False)
