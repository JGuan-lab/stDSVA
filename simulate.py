import pandas as pd
import numpy as np

def Simulated(spatial_rna, spatial_meta, spatial_loc, CoordinateXlable, CoordinateYlable, window,  CoordinateZlable = None):

  combined_spot = []
  combined_spot_loc = []
  window = window
  c = 0
  if CoordinateZlable != None:
    for z in np.arange(spatial_loc[CoordinateZlable].min(), spatial_loc[CoordinateZlable].max()):
      for x in np.arange((spatial_loc[CoordinateXlable].min()//window),spatial_loc[CoordinateXlable].max()//window+1):
        for y in np.arange((spatial_loc[CoordinateYlable].min()//window),spatial_loc[CoordinateYlable].max()//window+1):
          tmp_loc = spatial_loc[(x*window < spatial_loc[CoordinateXlable]) & (spatial_loc[CoordinateXlable] < (x+1)*window) & (y*window < spatial_loc[CoordinateYlable]) & (spatial_loc[CoordinateYlable] < (y+1)*window) & (spatial_loc[CoordinateZlable] == z)]
          if len(tmp_loc) > 0:
            c += 1
          combined_spot_loc.append([z,x,y])
          combined_spot.append(tmp_loc.index.to_list())
  else:
    for x in np.arange((spatial_loc[CoordinateXlable].min()//window),spatial_loc[CoordinateXlable].max()//window+1):
      for y in np.arange((spatial_loc[CoordinateYlable].min()//window),spatial_loc[CoordinateYlable].max()//window+1):
        tmp_loc = spatial_loc[(x*window < spatial_loc[CoordinateXlable]) & (spatial_loc[CoordinateXlable] < (x+1)*window) & (y*window < spatial_loc[CoordinateYlable]) & (spatial_loc[CoordinateYlable] < (y+1)*window)]
        if len(tmp_loc) > 0:
          c += 1
          combined_spot_loc.append([x,y])
          combined_spot.append(tmp_loc.index.to_list())


  combined_cell_counts = pd.DataFrame([len(s) for s in combined_spot],columns=['cell_count'])
  print ('The simulated spot has cells with ' + str(combined_cell_counts.min()[0]) + ' to ' + str(combined_cell_counts.max()[0]))
  combined_spot_loc = pd.DataFrame(combined_spot_loc, columns=['fov','x','y']) if CoordinateZlable != None else pd.DataFrame(combined_spot_loc, columns=['x','y'])


  combined_spot_exp = []
  for s in combined_spot:
    combined_spot_exp.append(spatial_rna.loc[s,:].sum(axis=0).values)
  combined_spot_exp = pd.DataFrame(combined_spot_exp, columns=spatial_rna.columns)

  combined_spot_clusters = pd.DataFrame(np.zeros((len(combined_spot_loc.index),len(np.unique(spatial_meta['cellType'])))),columns=np.unique(spatial_meta['cellType']))
  for i,c in enumerate(combined_spot):
    for clt in spatial_meta.loc[c,'cellType']:
      combined_spot_clusters.loc[i,clt] += 1
  print ('The simulated spot has size ' + str(combined_spot_clusters.shape[0]))
  print ('One simulated spot contains ' + str(combined_cell_counts.mean()))

  return combined_cell_counts, combined_spot_loc, combined_spot_exp, combined_spot_clusters

