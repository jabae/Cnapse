"""
Synapse partner assignment script
"""

import argparse
from sys import argv

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from cloudvolume import CloudVolume


## Partner assign
def assign(syn_seg, seg, res, r, mode="fixed"):
	"""
	INPUT:
	syn_seg : Synapse segmentation volume
	seg : Cell segmentation volume
	res : (x, y, z) resolution in nm
	r : Radius of area for partner assignment in nm
	mode : fixed or adaptive assignment
	
	OUTPUT:
	pre_list = Presynaptic cells
	post_partner_list = Postsynaptic cells for each synapse
	"""

	res = np.array(res)

	# Synapse segment list
	syn_seg_list = np.unique(syn_seg)
	syn_seg_list = syn_seg_list[1:] # Exclude 0
	nsyn = syn_seg_list.shape[0]

	vol_size = syn_seg.shape
	pre_list = np.zeros(N)
	post_partner_list = []

	for i in range(nsyn):
    
    syn_segid = syn_seg_list[i]
    syn_mask = syn_seg==syn_segid
    
    # Assign pre id
    segover_list, segover_size = np.unique(seg[syn_mask], return_counts=True) 
    if segover_list[0] == 0: # Remove 0
      segover_list = segover_list[1:]
      segover_size = segover_size[1:]
    
    # Assign 0 if there's no overlap
    if segover_list.shape[0] == 0:
      pre_list[i] = 0
    else:
      pre_list[i] = int(segover_list[np.argmax(segover_size)])
    
    ## Assign post ids
    x, y, z = np.where(syn_mask)
    syn_coord = np.c_[x,y,z]
    
    # Make bounding box for the synapse segment
    coord_max = np.max(syn_coord, axis=0)
    coord_min = np.min(syn_coord, axis=0)
    
    x, y, z = np.mgrid[np.max([coord_min[0]-r//res[0],0]):np.min([coord_max[0]+r//res[0],vol_size[0]]),
                       np.max([coord_min[1]-r//res[1],0]):np.min([coord_max[1]+r//res[1],vol_size[1]]),
                       np.max([coord_min[2]-r//res[2],0]):np.min([coord_max[2]+r//res[2],vol_size[2]])]
 
    # Search voxels within r nm from the synapse 
    pgrid = np.c_[x.ravel(),y.ravel(),z.ravel()]*res
    tree = cKDTree(pgrid)
    
    idx_list = tree.query_ball_point(syn_coord*res, r=r)
    idx_list = np.sum(idx_list)
    idx_list = np.unique(idx_list)
    
    # Extract post segment ids
    pval = (pgrid[idx_list,:]/res).astype("int")
    post_seg = seg[pval[:,0],pval[:,1],pval[:,2]]
    post_seg = np.unique(post_seg)
    
    # Remove 0 and autapse
    valid = ~np.isin(post_seg, [0, pre_list[i]])
    post_seg = post_seg[valid]
    post_partner_list.append(post_seg)
    
    if (i+1)% 100 == 0:
    	print(f"{i+1} / {nsyn} assigned.")

  return syn_seg_list, pre_list, post_partner_list

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

def save_synapse(syn_assign_info, outpath):

	n = syn_df.shape[0]

	syn_list = syn_assign_info[0]
	pre_list = syn_assign_info[1]
	post_list = syn_assign_info[2]

	syn_id_list = []
	pre_id_list = []
	post_id_list = []
	# x_pos_list = []
	# y_pos_list = []
	# z_pos_list = []
	# size_list = []
	for i in range(n):
	    
    syn_id = syn_list[i]
    pre_id = pre_list[i]
    post_id = post_list[i]
    
    for j in range(len(post_id)):
        
      post_id_list.append(post_id[j])
      syn_id_list.append(syn_df.iloc[i]["syn_id"])
      pre_id_list.append(syn_df.iloc[i]["pre_id"]))
        
	data = {"syn_id": syn_id_list,
	       "pre_id": pre_id_list,
	       "post_id": post_id_list}
	        
	syn_df = pd.DataFrame(data=data)
	syn_df.to_csv(outpath, index=False)


if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument("--syn_seg", required=True, type=str,
		help="Synapse segmentation volume")
	parser.add_argument("--cell_seg", required=True, type=str,
		help="Cell segmentation volume")
	parser.add_argument("--mip", required=True, type=int,
		help="Mip level number (2^[mip] nm)")
	parser.add_argument("--r", required=True, type=int,
		help="Radius of area for partner assignment (in nm)")
	parser.add_argument("--mode", required=False, type=str,
		default="fixed", help="Assignment mode (fixed or adaptive)")
	parser.add_argument("--outpath", required=False, type=str,
		default="synapse_list.csv", help="Output synapse file path")

	opt = parser.parse_args()

	syn_seg_path = opt.syn_seg
	seg_path = opt.cell_seg

	mip = opt.mip
	res = (2^mip, 2^mip, 50)

	outpath = opt.outpath

	# Load volumes
	syn_seg = load_volume(syn_seg_path, res)
	seg = load_volume(seg_path, res)

	# Assign partners
	syn_assign_info = assign(syn_seg, seg, res, r)

	# Save synapse list
	save_synapse(syn_assign_info, outpath)

