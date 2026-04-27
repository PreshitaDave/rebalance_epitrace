import sys
import matplotlib.pyplot as plt
import scvelo as scv
import scanpy as sc
import cellrank as cr
from cellrank.external.kernels import WOTKernel
from cellrank.tl.kernels import ConnectivityKernel
from cellrank.tl.estimators import GPCCA
from cellrank.tl.kernels import VelocityKernel
from typing import Any
from copy import copy
from anndata import AnnData
import numpy as np
import scipy.sparse as sp
import petsc4py


# definition
class AgeKernel(cr.tl.kernels.Kernel):
  def __init__(self, adata: AnnData, obs_key: str = "EpiTraceAge_iterative", **kwargs:Any): 
    super().__init__(adata=adata, obs_key=obs_key, **kwargs)
  def _read_from_adata(self, obs_key: str, **kwargs:Any) -> None:
    super()._read_from_adata(**kwargs)
    print(f"Reading `adata.obs[{obs_key!r}]`")
    self.pseudotime = 1 -  self.adata.obs[obs_key].values
  def compute_transition_matrix(self, some_parameter: float = 0.5) -> "AgeKernel":
    if self._reuse_cache({"some_parameter": some_parameter}):
      print("Using cached values for parameters:", self.params)
      return self
    transition_matrix = sp.diags((some_parameter,) * len(self.adata), dtype=np.float64)
    self._compute_transition_matrix(transition_matrix, density_normalize=True)
    return self
  def copy(self) -> "AgeKernel":
    return copy(self)
      


# set verbosity levels
cr.settings.verbosity = 2
scv.settings.verbosity = 3

# figure settings
scv.settings.set_figure_params("scvelo", dpi_save=400, dpi=80, transparent=True, fontsize=20, color_map="viridis")

# read and set data
adata = cr.read('./Packed_Figure_Reproducibility/Source_Data_for_Figure/Figure6/SHARESEQ_brain_epitrace_obj_age_estimated_multiome.h5ad')
adata.layers['spliced'] = adata.X
seta = [0,1,2,3,4,5,6,7,]
sub_adata = adata[ adata.obs['Cluster.Name']!=8 ].copy()
adata = sub_adata
adata.obs['Real.Name'] = ''
adata.obs['Real.Name'][adata.obs['Cluster.Name'].isin([0])] = 'Radial Glia'
adata.obs['Real.Name'][adata.obs['Cluster.Name'].isin([1])] = 'Cyc. Prog'
adata.obs['Real.Name'][adata.obs['Cluster.Name'].isin([2])] = 'nIPC'
adata.obs['Real.Name'][adata.obs['Cluster.Name'].isin([3])] = 'ExN2'
adata.obs['Real.Name'][adata.obs['Cluster.Name'].isin([4])] = 'ExN3'
adata.obs['Real.Name'][adata.obs['Cluster.Name'].isin([5])] = 'ExN4'
adata.obs['Real.Name'][adata.obs['Cluster.Name'].isin([6])] = 'ExN5'
adata.obs['Real.Name'][adata.obs['Cluster.Name'].isin([7])] = 'mGPC/OPC'
adata.obs['Real.Name'][adata.obs['Cluster.Name'].isin([8])] = 'Subplate'
adata.obs['Real.Name'][adata.obs['Cluster.Name'].isin([9])] = 'InN3'
adata.obs['Real.Name'][adata.obs['Cluster.Name'].isin([10])] = 'InN2'
adata.obs['Real.Name'][adata.obs['Cluster.Name'].isin([11])] = 'InN1'
adata.obs['Real.Name'][adata.obs['Cluster.Name'].isin([12])] = 'Endo'

# preprocessing
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata, random_state=0)
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)




# kernels
ak = AgeKernel(adata).compute_transition_matrix()
vk = cr.kernels.VelocityKernel(adata).compute_transition_matrix()
ctk = cr.tl.kernels.CytoTRACEKernel(adata).compute_transition_matrix(threshold_scheme="soft", nu=0.5)


color_list_a = ["#FFA500","#8B0000","#98F5FF","#8EE5EE","#7AC5CD","#53868B","#A2CD5A","#BCEE68","#CAFF70","#FF0000","#68228B","#6495ED"]
color_list_macro = ["#8B0000","#9932CC","#68228B","#FF0000","#FFC0CB","#8EE5EE","#BCEE68","#CAFF70","#53868B","#FFA500"]
color_list_terminal = ["#8B0000","#9932CC","#68228B","#53868B","#FFA500"]



# project kernels 
scv.pl.velocity_embedding_stream(adata, color=['Real.Name'],palette=color_list_a, basis="umap", legend_loc="right",save='Figure5b_Brain_RNA_Velocity_only_streamline.png', title='RNA velocity', density=3,dpi=800,size=28,figsize=[8,8],show=False)

ctk.compute_projection(basis="umap")
scv.pl.velocity_embedding_stream(adata, color=['Real.Name'],palette=color_list_a, vkey="T_fwd", basis="umap", legend_loc="right",save='Figure5b_Brain_CytoTrace_only_streamline.png', title='CytoTRACE', density=3,dpi=800,size=28,figsize=[8,8],show=False)

ak.compute_projection(basis="umap")
scv.pl.velocity_embedding_stream(adata, color="EpiTraceAge_iterative", vkey="T_fwd", basis="umap", legend_loc="right",save='Figure5b_Brain_Epitrace_Age_only_streamline.png', title='EpiTrace', density=3,dpi=800,size=28,color_map='autumn',figsize=[8,8],show=False)




# combine kernels
combined_kernel = 0.6 * vk + 0.2 * ak + 0.2 * ctk
combined_kernel.compute_transition_matrix()
print(combined_kernel)
combined_kernel.compute_projection(basis="umap")
scv.pl.velocity_embedding_stream(adata, basis='umap',  vkey="T_fwd",color=['Real.Name'], save='Figure5a_Brain_Combine_Epitrace_Cytotrace_Velocity_streamline.png', title='Combined kernels',density=3,dpi=800,size=28,palette=color_list_a,figsize=[8,8],show=False, legend_loc='right')



