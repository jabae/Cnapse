# Cnapse
Automated synapse detection in C. elegans EM images.

## Synpase assignment
````
assign.py --syn_seg [syn_seg.tif] --cell_seg [cell_seg.tif] --mip [mip] --r [r] --mode [mode] --outpath [outpath.tif]

```
- `syn_seg` : Synpase segmentation volume
- `cell_seg` : Cell segmentation volume
- `mip` : Mip level of the volumes (2^[mip] nm resolution)
- `r` : Radius in [nm] for partner assignment
- `mode` : Assignment mode (`fixed` or `adaptive`)
- `outpath` : Path to save result