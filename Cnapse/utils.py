import numpy as np
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