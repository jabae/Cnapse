"""
Neurotransmitter diffusion simulation
"""

import argparse
from sys import argv

import numpy as np
import pandas as pd
from skimage.morphology import dilation, ball

from cloudvolume import CloudVolume


## Fibonacci spiral model
def fibonacci_spiral_sphere(n_points):

 	vec = np.zeros((n_points, 3))
 	gr = (np.sqrt(5) + 1)/2 # Golden ratio = 1.6180...

 	p = np.arange(n_points)

 	u = 2*np.pi*gr*p
 	v = 1 - (2*p + 1)/n_points
 	sinTheta = np.sqrt(1 - v**2)

 	vec[:,0] = np.cos(u)*sinTheta
 	vec[:,1] = np.sin(u)*sinTheta
 	vec[:,2] = v

 	return vec


if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument("--syn_seg", required=True, type=str,
		help="Synapse segmentation volume")
	parser.add_argument("--cell_seg", required=True, type=str,
		help="Cell segmentation volume")
	parser.add_argument("--mip", required=True, type=int,
		help="Mip level number (2^[mip] nm)")
	parser.add_argument("--syn_info", required=True, type=str,
		help="Synapse info")

	opt = parser.parse_args()

	syn_seg_path = opt.syn_seg
	seg_path = opt.cell_seg

	mip = opt.mip
	res = (2^mip, 2^mip, 50)
	
	# Load volumes
	syn_seg = load_volume(syn_seg_path, res)
	seg = load_volume(seg_path, res)

	volsize = syn_seg.shape

	opt.diffusionsteplength_nm = 4
	opt.nrdiffusionvectors = 10000
	opt.nrofparticles = 1000
	opt.voxelsize = np.array([4,4,50])
	opt.maxiterations=10000
	opt.synexpandradiusxy_px = 6

	# Precompute diffusion vectors
	dvec = fibonacci_spiral_sphere(opt.nrdiffusionvectors)

	ds = opt.diffusionsteplength_nm/opt.voxelsize
	dvec = dvec*ds

	# Dilation params
	seg_se = ball(1)
	syn_se = np.zeros((3,3,3))
	syn_se[:,:,1] = 1

	# Synapse info
	syn_info_df = pd.read_csv(opt.syn_info)

	syn_id_list = np.array(syn_info_df["syn_id"])
	pre_id_list = np.array(syn_info_df["pre_id"])
	nsyn = syn_id_list.shape[0]

	errcount = 0

	for i in range(nsyn):

		syn_id = syn_id_list[i]
		pre_id = pre_id_list[i]

		syn_mask = (syn_seg==syn_id).astype("uint8")

		esyn = dilation(syn_mask, syn_se)
		pseg = (seg==pre_id).astype("uint8")
		epseg = dilation(pseg, seg_se)
		
		exc_valid = (seg!=0)*(esyn==0)
		epseg[exc_valid] = 0
		
		emittervoxels_loc = np.where(epseg)
		nemitvx = emittervoxels_loc[0].shape[0]

		if nemitvx==0:
		
			print(f"ERROR: Synapse with id {syn_id} has no emitter voxels!")
			errcount = errcount + 1

		else:

			sourcevoxnr = np.random.randint(nemitvx, size=opt.nrofparticles)
			x = emittervoxels_loc[sourcevoxnr]
			y = emittervoxels_loc[sourcevoxnr]
			z = emittervoxels_loc[sourcevoxnr]

			nt = np.zeros((nemitvx, 5))
			nt[:,0] = x
			nt[:,1] = y
			nt[:,2] = z
			nt[:,3] = np.zeros(nemitvx)
			nt[:,4] = np.ones(nemitvx)

			count = 0
			while (count<opt.maxiterations) and (np.sum(nt[:,4])>0):

				vnr = np.random.randint(opt.nrdiffusionvectors, size=opt.nrofparticles)
				v = dvec[vnr,:]*nt[:,4]
				nt[:,:3] = nt[:,:3] + v

				p = np.round(nt[:,:3])
				minp = np.min(p, axis=0)
				maxp = np.max(p, axis=0)

				if (np.min(minp)<0) or np.max(maxp-np.array(syn_vol.shape)):
					
					exited = (np.min(p, axis=1)<0) + (p[:,0]>=volsize[0]) + (p[:,1]>=volsize[1]) + (p[:,2]>=volsize[2])
					nt[exited,:3] = nt[exited,:3] - v[exited,:] # move back
					nt[exited,4] = 0
					p = np.round(nt[:,:3])

				vids = seg[p[:,0], p[:,1], p[:,2]]

				# Particles which moved to the presynaptic side
				vids_pre = vids==pre_id
				nt[vids_pre,:3] = nt[vids_pre,:3]-v[vids_pre,:]
				vids[vids_pre] = 0

				sp = np.where(vids>0)[0]
				nt[sp, 3] = vids[sp]
				nt[sp, 4] = 0

				count += 1