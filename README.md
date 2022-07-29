# Cnapse
Semi-automated synapse detection in *C*. *elegans* electron microscopy (EM) images.

## Synapse detection pipeline
![](figures/synapse_detection.png)

*Cnapse* takes raw EM image as an input and use trained convolutional neural network (CNN) to detect the synapse candidates. Then, synaptic partners are assigned by running Monte Carlo simulation of neurotransmitters (Witvliet et al., 2021) and the size of the synaptic connections are determined by the proportion of neurotransmitters.

## Synpase assignment
```
assign.py --syn_seg [syn_seg.tif] --cell_seg [cell_seg.tif] --mip [mip] --r [r] --mode [mode] --outpath [outpath.tif]
```
- `syn_seg` : Synpase segmentation volume
- `cell_seg` : Cell segmentation volume
- `mip` : Mip level of the volumes (2^[mip] nm resolution)
- `r` : Radius in [nm] for partner assignment
- `mode` : Assignment mode (`fixed` or `adaptive`)
- `outpath` : Path to save result
