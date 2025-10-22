
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import scanpy as sc
import paste as pst
import anndata as ad
import pandas as pd
def load_slices(data_dir, slice_names=["slice1", "slice2"]):
    slices = []  
    for slice_name in slice_names:
        slice_i = sc.read_csv(data_dir + slice_name + ".csv")
        slice_i_coor = np.genfromtxt(data_dir + slice_name + "_coor.csv", delimiter = ',')
        slice_i.obsm['spatial'] = slice_i_coor
        slices.append(slice_i)
    return slices



def bulidSlice(slice_i,slice_i_coor):
    print("build slice......")
    slice = ad.AnnData(slice_i)
    slice.obsm['spatial'] = slice_i_coor
    print("done!")

    return slice

def getSlices(slice_pseudo, slice_real, slice_pseudo_coor, slice_real_coor):
    slices = []
    slices.append(bulidSlice(slice_pseudo, slice_pseudo_coor))
    slices.append(bulidSlice(slice_real, slice_real_coor))
    print("done!")
    return slices

def plotSlices(slices):
    print("plot the slices....................")
    slice_colors = ['#e41a1c','#377eb8']
    plt.figure(figsize=(7,7))
    for i in range(len(slices)):
        pst.plot_slice(slices[i],slice_colors[i],s=400)
    plt.legend(handles=[mpatches.Patch(color=slice_colors[0], label='1'),mpatches.Patch(color=slice_colors[1], label='2')])
    plt.gca().invert_yaxis()
    plt.axis('off')
    plt.show()




def getAnnMerge(slices):
    print("get ann begins............")
    merged_adata = ad.concat(slices, axis=0, join='outer') 
    num_pseudo_spots = slices[0].X.shape[0]
    num_real_spots = slices[1].X.shape[0]
    sample_types = pd.Series(index=merged_adata.obs_names, dtype=str)  
    sample_types.iloc[:num_pseudo_spots] = "S1"
    sample_types.iloc[num_real_spots:] = "S3"
    merged_adata.obs['data'] = sample_types
    return merged_adata





def pasteMy(slice_pseudo, slice_real, slice_pseudo_coor, slice_real_coor):
    print("get slices....................")
    slices = getSlices(slice_pseudo, slice_real, slice_pseudo_coor, slice_real_coor)
    print("get slices done!")
    # print("begin plot original slices.............")
    # plotSlices(slices)
    slice1, slice2 = slices
    initial_slice = slice1.copy()   
    lmbda = len(slices)*[1/len(slices)]
 
    initial_slice.obsm['spatial'] = np.array(initial_slice.obsm['spatial'])
    slices[0].obsm['spatial'] = np.array(slices[0].obsm['spatial'])
    slices[1].obsm['spatial'] = np.array(slices[1].obsm['spatial'])
    print("OK!center_align begins...................")
    center_slice, pis = pst.center_align(initial_slice, slices, lmbda, use_gpu = True)
    new_center, new_slices = pst.stack_slices_center(center_slice, slices, pis)
    # plotSlices(new_slices)
    return(getAnnMerge(new_slices))




def pasteGetSpatial(slice_pseudo, slice_real, slice_pseudo_coor, slice_real_coor):
    adata = pasteMy(slice_pseudo, slice_real, slice_pseudo_coor, slice_real_coor)
    spatial = pd.DataFrame(adata.obsm['spatial'])
    spatial.index = adata.obs_names
    spatial.columns = ['X','Y']
    return(spatial)



