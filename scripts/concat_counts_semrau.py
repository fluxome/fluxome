import anndata
import pandas as pd
import numpy as np 
import scipy

from scipy.io import mmread

srrs = ['SRR3290185', 'SRR3290188','SRR3290190','SRR3290193','SRR3290196',
'SRR3290199','SRR3290202','SRR3290203','SRR3290204']
times = [2,6,12,24,36,48,60,72,96]
data_path = '/home/tchari/CellDynamicsData/harissa_timecourse/'

var_names = pd.read_csv(data_path+srrs[0]+'_outs/Solo.out/Velocyto/filtered/features.tsv',sep='\t',header=None)[0]

u_mats = []
s_mats = []
bars = []
all_ts = []
for i in range(len(srrs)):

	bs = pd.read_csv(data_path+srrs[i]+'_outs/Solo.out/Velocyto/filtered/barcodes.tsv',sep='\t',header=None)[0]
	bars += list(bs)
	all_ts += [times[i]]*len(bs)

	#Get U,S, then concatenate at end
	u_mats += [mmread(data_path+srrs[i]+'_outs/Solo.out/Velocyto/filtered/unspliced.mtx')] #gene x cell
	s_mats += [mmread(data_path+srrs[i]+'_outs/Solo.out/Velocyto/filtered/spliced.mtx')]


U = scipy.sparse.hstack(u_mats)
S = scipy.sparse.hstack(s_mats)	

U = np.asarray(U.T.todense())
S = np.asarray(S.T.todense())

adata = anndata.AnnData(X=S) #could also leave as sparse?
adata.var_names = var_names
adata.obs_names = bars
adata.layers['unspliced'] = U
adata.layers['spliced'] = S
adata.obs['time'] = all_ts
adata.obs['n_cells'] = [1]*len(adata)
print(adata)

adata.write_h5ad(data_path+'all_times_semrau.h5ad')
