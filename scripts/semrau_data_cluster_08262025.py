#!/usr/bin/env python
# coding: utf-8

# In[38]:


#extraembryonic endoderm (XEN): Sparc, Col4a1, Lama1, Dab2
#neuro-ectodermal: Prtg, Mdk, Fabp5, Cd24

#Test k-means clustering
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import seaborn as sns
import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

kmeans = KMeans(n_clusters=3, random_state=10)

data_path = '/home/tchari/CellDynamicsData/harissa_timecourse/'
adata = anndata.read_h5ad(data_path+'all_times_semrau.h5ad')
sub = adata[adata.X.sum(axis=1)>1e3,:] #[adata.obs['time'] == 96]
sub = sub[:,sub.X.sum(axis=0)>50] #> currently just spliced
gs = list(sub.var_names)

#From Semrau et al Supp Fig 3a
pc1_genes = [
    "Sparc", "Col4a1", "Lama1", "Lamb1", "Col4a2", "Lamc1", "Serpinh1", "Dab2", "Lrpap1",
    "Srgn", "Plod2", "Fst", "Pth1r", "Cryab", "Cubn", "Cyp26a1", "F3", "H19", "Sox17",
    "Gnas", "Krt8", "Nrk", "Flrt3", "Adam19", "Msl3"
]
pc2_genes = [
    "Prtg", "Ccnd2", "Mdk", "Crabp1", "Smoc2", "Mest", "Peg10", "Peg3", "Frem2", "Angptl1",
    "Peg3as", "Nepn", "Igf2", "Fabp5", "Ednrb", "Foxa2", "Tpm1", "Fat3", "Nr2f1", "Gnas",
    "Cd24a", "Malat1", "Id2", "Clu", "Flrt3"
]

print('sub adata: ',sub)


# In[ ]:





# In[ ]:





# In[39]:


print(sub.X.shape)
row_sums = sub.X.sum(axis=1, keepdims=True)
scaling_factors = np.median((np.sum(sub.X,axis=1))) / row_sums
normed = sub.X * scaling_factors

variances = np.sqrt(np.var(sub.X, axis=0))/np.mean(sub.X, axis=0) #do this pre-normalization
sorted_indices = np.argsort(variances)[::-1]
top_1k = np.array(gs)[sorted_indices[:1000]]



# In[40]:


plt.scatter(np.mean(sub.X, axis=0),variances,alpha=0.6)
plt.scatter(np.mean(sub.X, axis=0)[sorted_indices[:1000]],variances[sorted_indices[:1000]],color='red')
plt.loglog()


# In[41]:


#Add genes from first harissa paper
genes = ["Ascl1", "Cd24a", "Cdh2", "Cer1", "Col4a2", "Dab2", "Dnmt3a", "Dnmt3b",
    "Dppa2", "Dppa3", "Dppa4", "Esrrb", "Fgf4", "Fgf5", "Fst", "Gata4", "Gata6",
    "Hoxa1", "Hoxb2", "Jarid2", "Klf2", "Klf4", "Krt8", "Lamb1", "Meis1", "Nes",
    "Nr0b1", "Pax6", "Pdgfra", "Pou3f1", "Pou5f1", "Sox1", "Sox11", "Sox17",
    "Sox2", "Sox9", "Sparc", "Wnt8a", "Zfp42", "Zic2", "Zic3"] 



gs_for_filt = list(top_1k)+genes+pc1_genes+pc2_genes #pc1_genes+pc2_genes+genes  #genes+list(top_1k)
final_X = normed[:,pd.Categorical(gs).isin(gs_for_filt)]
final_gs = list(np.array(gs)[pd.Categorical(gs).isin(gs_for_filt)])

scaler = StandardScaler()
scaled = scaler.fit_transform(np.log1p(final_X))


pca = PCA(n_components=3)
X_pca = pca.fit_transform(np.log1p(final_X)) #np.log1p(final_X)

#Kmeans on pc genes from Semrau + variable genes + harissa genes
kmeans.fit(np.log1p(final_X)) #X_pca, np.log1p(final_X), scaled
assignments = kmeans.labels_


# In[ ]:





# In[43]:


#do the harissa/semrau genes meet expression threshold?
uniq = np.unique(pc1_genes+pc2_genes+genes)
frac = np.sum([i in gs for i in uniq])/len(uniq)
print(f'{frac*100:.2f}'+'% of PC/semrau genes in expr filter')


# In[52]:


#Plots
fig, axs = plt.subplots(4, 2,figsize=(6.5, 8))
ax = axs.flat

g=sns.scatterplot(x=X_pca[:, 0], y=X_pca[:, 1],hue=assignments,alpha=0.5,ax=ax[0])
g.legend(loc='right', bbox_to_anchor=(1.25, 0.5), ncol=1)
ax[0].set_xlabel('PC 1')
ax[0].set_ylabel('PC 2')
ax[0].set_title('Kmeans on PC genes')
ax[0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


g=sns.scatterplot(x=X_pca[:, 0], y=X_pca[:, 1],hue=pd.Categorical(list(sub.obs['time'])),alpha=0.5,ax=ax[1])
g.legend(loc='right', bbox_to_anchor=(1.25, 0.5), ncol=1)
ax[1].set_xlabel('PC 1')
ax[1].set_ylabel('PC 2')
ax[1].set_title('True times')
ax[1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


g=sns.scatterplot(x=X_pca[:, 0], y=X_pca[:, 1],hue=np.log1p(final_X[:,final_gs.index('Sparc')]),alpha=0.5,ax=ax[2])
g.legend(loc='right', bbox_to_anchor=(1.25, 0.5), ncol=1)
ax[2].set_xlabel('PC 1')
ax[2].set_ylabel('PC 2')
ax[2].set_title('Sparc-Xen log')
ax[2].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)



g=sns.scatterplot(x=X_pca[:, 0], y=X_pca[:, 1],hue=np.log1p(final_X[:,final_gs.index('Col4a1')]),alpha=0.5,ax=ax[3])
g.legend(loc='right', bbox_to_anchor=(1.25, 0.5), ncol=1)
ax[3].set_xlabel('PC 1')
ax[3].set_ylabel('PC 2')
ax[3].set_title('Col4a1-Xen log')
ax[3].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


g=sns.scatterplot(x=X_pca[:, 0], y=X_pca[:, 1],hue=np.log1p(final_X[:,final_gs.index('Prtg')]),alpha=0.5,ax=ax[4])
g.legend(loc='right', bbox_to_anchor=(1.25, 0.5), ncol=1)
ax[4].set_xlabel('PC 1')
ax[4].set_ylabel('PC 2')
ax[4].set_title('Prtg-Ect log')
ax[4].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


g=sns.scatterplot(x=X_pca[:, 0], y=X_pca[:, 1],hue=np.log1p(final_X[:,final_gs.index('Ccnd2')]),alpha=0.5,ax=ax[5])
g.legend(loc='right', bbox_to_anchor=(1.25, 0.5), ncol=1)
ax[5].set_xlabel('PC 1')
ax[5].set_ylabel('PC 2')
ax[5].set_title('Ccnd2-Ect log')
ax[5].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


g=sns.scatterplot(x=X_pca[:, 0], y=X_pca[:, 1],hue=np.log1p(final_X[:,final_gs.index('Fgf4')]),alpha=0.5,ax=ax[6])
g.legend(loc='right', bbox_to_anchor=(1.25, 0.5), ncol=1)
ax[6].set_xlabel('PC 1')
ax[6].set_ylabel('PC 2')
ax[6].set_title('Fgf4-Plur log')
ax[6].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


g=sns.scatterplot(x=X_pca[:, 0], y=X_pca[:, 1],hue=np.log1p(final_X[:,final_gs.index('Nr0b1')]),alpha=0.5,ax=ax[7])
g.legend(loc='right', bbox_to_anchor=(1.25, 0.5), ncol=1)
ax[7].set_xlabel('PC 1')
ax[7].set_ylabel('PC 2')
ax[7].set_title('Nr0b1-Plur log')
ax[7].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


plt.tight_layout()
plt.savefig('kmeans_expr_semrau_data.png')


#Fig 2c Semrau et al recreation plots

#mesc - Fgf4, Nr0b1, Klf4, Klf2, Jarid2, Zfp42, Dppa3, Esrrb, Sox2, Pou5f1
#ectoderm - Prtg, Ccnd2 
#Xen -  Sparc, Col4a1


# In[54]:


np.unique(assignments)


# In[71]:


np.unique(pd.Categorical(list(sub.obs['time'])))


# In[58]:


import scipy
import scipy.stats as st
from matplotlib.patches import Patch


# In[122]:


#Make unsupervised hierarhical clustering clustermap --> get high-level assignments for cells


gs_for_filt = list(top_1k)+genes+pc1_genes+pc2_genes#pc1_genes+pc2_genes+genes  #genes+list(top_1k)
final_X = normed[:,pd.Categorical(gs).isin(gs_for_filt)]
final_gs = list(np.array(gs)[pd.Categorical(gs).isin(gs_for_filt)])

plt.figure()

c_types = assignments
lut = dict(zip(pd.Categorical(np.unique(assignments)), "rbgmyc"))
k_colors = pd.Categorical(c_types).map(lut)

log_counts=np.log1p(final_X)
centered=log_counts-log_counts.mean(0)-log_counts.mean(1)[:,None]
normed_values=centered

row_colors = pd.DataFrame({'kmeans':np.array(k_colors)})

sns.clustermap(pd.DataFrame(normed_values),figsize = (14,4),cmap='bwr',row_colors=row_colors,
              row_cluster=True,col_cluster=True,metric='correlation',xticklabels=final_gs,yticklabels=[])	

plt.tick_params(rotation=90)

handles = [Patch(facecolor=lut[name]) for name in lut]
plt.legend(handles, lut.keys(), title='K Means',
   bbox_to_anchor=(1.1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')

plt.show()



# In[121]:


#Make unsupervised hierarhical clustering clustermap --> get high-level assignments for cells


gs_for_filt = genes+pc1_genes+pc2_genes#pc1_genes+pc2_genes+genes  #genes+list(top_1k)
final_X = normed[:,pd.Categorical(gs).isin(gs_for_filt)]
final_gs = list(np.array(gs)[pd.Categorical(gs).isin(gs_for_filt)])

plt.figure()

c_types = assignments
lut = dict(zip(pd.Categorical(np.unique(assignments)), "rbgmyc"))
k_colors = pd.Categorical(c_types).map(lut)

log_counts=np.log1p(final_X)
centered=log_counts-log_counts.mean(0)-log_counts.mean(1)[:,None]
normed_values=centered

row_colors = pd.DataFrame({'kmeans':np.array(k_colors)})

sns.clustermap(pd.DataFrame(normed_values),figsize = (14,4),cmap='bwr',row_colors=row_colors,
              row_cluster=True,col_cluster=True,metric='correlation',xticklabels=final_gs,yticklabels=[])	

plt.tick_params(rotation=90)

handles = [Patch(facecolor=lut[name]) for name in lut]
plt.legend(handles, lut.keys(), title='K Means',
   bbox_to_anchor=(1.1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')

plt.show()



# In[88]:


#Make unsupervised hierarhical clustering clustermap --> get high-level assignments for cells


gs_for_filt = genes+pc1_genes+pc2_genes#pc1_genes+pc2_genes+genes  #genes+list(top_1k)
final_X = normed[:,pd.Categorical(gs).isin(gs_for_filt)]
final_gs = list(np.array(gs)[pd.Categorical(gs).isin(gs_for_filt)])

plt.figure()

c_types = assignments
lut = dict(zip(pd.Categorical(np.unique(assignments)), "rbgmyc"))
k_colors = pd.Categorical(c_types).map(lut)

c_types = list(sub.obs['time'])
lut = dict(zip(pd.Categorical(np.unique(list(sub.obs['time']))), ['#deebf7', '#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b', '#041f40']))
t_colors = pd.Categorical(c_types).map(lut)


log_counts=np.log1p(final_X)
centered=log_counts-log_counts.mean(0)-log_counts.mean(1)[:,None]
normed_values=centered

row_colors = pd.DataFrame({'kmeans':np.array(k_colors),
                          'time':np.array(t_colors)})

sns.clustermap(pd.DataFrame(normed_values),figsize = (14,4),cmap='bwr',row_colors=row_colors,
              row_cluster=True,col_cluster=True,metric='correlation',xticklabels=final_gs,yticklabels=[])	

plt.tick_params(rotation=90)

handles = [Patch(facecolor=lut[name]) for name in lut]
plt.legend(handles, lut.keys(), title='K Means',
   bbox_to_anchor=(1.3, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')

plt.show()



# In[117]:


#Make unsupervised hierarhical clustering clustermap --> get high-level assignments for cells

#new genes from Fig 3d semrau

new_plur = ["Pou5f1", "Sox2", "Esrrb", "Dppa3", "Zfp42", "Jarid2", "Klf2", "Klf4"]
new_postimp = ["Nr0b1", "Fgf4", "Zic1", "Zic3", "Dppa2", "Dppa4", "Dnm1l", "Dnmt3a", "Dnmt3b", 
    "Fgf5", "Fgf1", "Pou3f1", "Cer1", "Wnt8a", "Cd24a", "Sox1", "Hoxb2", "Meis1"]
new_ect = ["Nes", "Klf1", "Ascl1", "Pax6", "Sox11", "Sox4", "Hoxa1", "Hoxc9", "Pdgfra",]
new_end = ["Fst", "Sparc", "Gata4", "Gata6", "Sox17", "Lrrb1", "Lamb1", "Dab2"]


gs_for_filt = genes+pc1_genes+pc2_genes+new_plur+new_postimp+new_ect+new_end #pc1_genes+pc2_genes+genes  #genes+list(top_1k)
final_X = normed[:,pd.Categorical(gs).isin(gs_for_filt)]
final_gs = list(np.array(gs)[pd.Categorical(gs).isin(gs_for_filt)])

plt.figure()

c_types = assignments
lut1 = dict(zip(pd.Categorical(np.unique(assignments)), "rbgmyc"))
k_colors = pd.Categorical(c_types).map(lut1)

c_types = list(sub.obs['time'])
lut2 = dict(zip(pd.Categorical(np.unique(list(sub.obs['time']))), ['#deebf7', '#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b', '#041f40']))
t_colors = pd.Categorical(c_types).map(lut2)


c_types = ['na']*len(genes+pc1_genes+pc2_genes)+['plur']*len(new_plur)+['post-impl']*len(new_postimp)+['neur ecto']*len(new_ect)+['endoderm']*len(new_end)
just_fig3 = ['plur']*len(new_plur)+['post-impl']*len(new_postimp)+['neur ecto']*len(new_ect)+['endoderm']*len(new_end)
just_fig3_gs = new_plur+new_postimp+new_ect+new_end 

lut3 = dict(zip(pd.Categorical(np.unique(c_types)), "bgmyc"))
g_colors = pd.Categorical([just_fig3[just_fig3_gs.index(i)] if i in just_fig3_gs else 'na' for i in final_gs]).map(lut3)


log_counts=np.log1p(final_X)
centered=log_counts-log_counts.mean(0)-log_counts.mean(1)[:,None]
normed_values=centered

row_colors = pd.DataFrame({'kmeans':np.array(k_colors),
                          'time':np.array(t_colors)})

col_colors = pd.DataFrame({'g_type':np.array(g_colors)})

g=sns.clustermap(pd.DataFrame(normed_values),figsize = (14,4),cmap='bwr',row_colors=row_colors,col_colors =col_colors,
              row_cluster=True,col_cluster=True,metric='correlation',xticklabels=final_gs,yticklabels=[])	

plt.tick_params(rotation=90)

handles1 = [Patch(facecolor=lut1[name]) for name in lut1]
g.ax_row_dendrogram.legend(
    handles1, lut1.keys(), title='K Means',
    bbox_to_anchor=(0.5, 1), loc='upper right'
)

# Add legend for 'time'
handles2 = [Patch(facecolor=lut2[name]) for name in lut2]
g.ax_row_dendrogram.legend(
    handles2, lut2.keys(), title='Time',
    bbox_to_anchor=(5.2, 1), loc='upper left'
)

# Add legend for 'g_type' (column annotation)
handles3 = [Patch(facecolor=lut3[name]) for name in lut3]
g.ax_col_dendrogram.legend(
    handles3, lut3.keys(), title='Gene Type',
    bbox_to_anchor=(0.5, 1.1), loc='lower center'
)

plt.show()



# In[113]:


#Make unsupervised hierarhical clustering clustermap --> get high-level assignments for cells

#new genes from Fig 3d semrau

new_plur = ["Pou5f1", "Sox2", "Esrrb", "Dppa3", "Zfp42", "Jarid2", "Klf2", "Klf4"]
new_postimp = ["Nr0b1", "Fgf4", "Zic1", "Zic3", "Dppa2", "Dppa4", "Dnm1l", "Dnmt3a", "Dnmt3b", 
    "Fgf5", "Fgf1", "Pou3f1", "Cer1", "Wnt8a", "Cd24a", "Sox1", "Hoxb2", "Meis1"]
new_ect = ["Nes", "Klf1", "Ascl1", "Pax6", "Sox11", "Sox4", "Hoxa1", "Hoxc9", "Pdgfra",]
new_end = ["Fst", "Sparc", "Gata4", "Gata6", "Sox17", "Lrrb1", "Lamb1", "Dab2"]


gs_for_filt = new_plur+new_postimp+new_ect+new_end #pc1_genes+pc2_genes+genes  #genes+list(top_1k)
final_X = normed[:,pd.Categorical(gs).isin(gs_for_filt)]
final_gs = list(np.array(gs)[pd.Categorical(gs).isin(gs_for_filt)])

plt.figure()

c_types = assignments
lut1 = dict(zip(pd.Categorical(np.unique(assignments)), "rbgmyc"))
k_colors = pd.Categorical(c_types).map(lut1)

c_types = list(sub.obs['time'])
lut2 = dict(zip(pd.Categorical(np.unique(list(sub.obs['time']))), ['#deebf7', '#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b', '#041f40']))
t_colors = pd.Categorical(c_types).map(lut2)


c_types = ['plur']*len(new_plur)+['post-impl']*len(new_postimp)+['neur ecto']*len(new_ect)+['endoderm']*len(new_end)
lut3 = dict(zip(pd.Categorical(np.unique(c_types)), "bgmyc"))
g_colors = pd.Categorical([c_types[gs_for_filt.index(i)] for i in final_gs]).map(lut3)


log_counts=np.log1p(final_X)
centered=log_counts-log_counts.mean(0)-log_counts.mean(1)[:,None]
normed_values=centered

row_colors = pd.DataFrame({'kmeans':np.array(k_colors),
                          'time':np.array(t_colors)})

col_colors = pd.DataFrame({'g_type':np.array(g_colors)})

g=sns.clustermap(pd.DataFrame(normed_values),figsize = (14,4),cmap='bwr',row_colors=row_colors,col_colors =col_colors,
              row_cluster=True,col_cluster=True,metric='correlation',xticklabels=final_gs,yticklabels=[])	

plt.tick_params(rotation=90)

handles1 = [Patch(facecolor=lut1[name]) for name in lut1]
g.ax_row_dendrogram.legend(
    handles1, lut1.keys(), title='K Means',
    bbox_to_anchor=(0.5, 1), loc='upper right'
)

# Add legend for 'time'
handles2 = [Patch(facecolor=lut2[name]) for name in lut2]
g.ax_row_dendrogram.legend(
    handles2, lut2.keys(), title='Time',
    bbox_to_anchor=(5.2, 1), loc='upper left'
)

# Add legend for 'g_type' (column annotation)
handles3 = [Patch(facecolor=lut3[name]) for name in lut3]
g.ax_col_dendrogram.legend(
    handles3, lut3.keys(), title='Gene Type',
    bbox_to_anchor=(0.5, 1.1), loc='lower center'
)

plt.show()



# In[ ]:





# In[144]:


assign_dict = dict(zip([0,1,2],['XEN','mESC','Ect']))
palette = {'XEN':'b','mESC':'darkorange','Ect':'r'}

fig, axs = plt.subplots(2, 5,figsize=(15,3))
ax = axs.flat
for i in range(len(np.unique(sub.obs['time']))):
    t = np.unique(sub.obs['time'])[i]

    assign_sub = assignments[list(sub.obs['time']==t)]
    names = [assign_dict[i] for i in assign_sub]
    ax[i].set_title(str(t)+'hr')
    g=sns.histplot(x=names,hue=names,palette=palette,ax=ax[i])
    sns.move_legend(ax[i], "upper left", bbox_to_anchor=(1, 1))

plt.tight_layout()


# In[45]:


# #Kmeans on variable genes like Semrau describe
# var_X = normed[:,pd.Categorical(gs).isin(list(top_1k)+genes+pc1_genes+pc2_genes)]
# #kmeans.fit(np.log1p(var_X)) #X_pca
# kmeans.fit(np.log1p(var_X))
# assign_var = kmeans.labels_


# #Check if k-means clustering changes much if using variable genes + PCA-Semrau genes
# plt.figure(figsize=(5,3))
# g=sns.scatterplot(x=X_pca[:, 0], y=X_pca[:, 1],hue=assign_var,alpha=0.6)
# g.legend(loc='right', bbox_to_anchor=(1.25, 0.5), ncol=1)
# plt.xlabel('PC 1')
# plt.ylabel('PC 2')
# plt.tight_layout()
# plt.title('Kmeans on Var+PC genes')
# plt.savefig('pca_kmeans_var_0826.png')


# **Save filtered and annotated data**

# In[145]:


sub


# In[148]:


np.unique(assignments)


# In[149]:


k_means_dict = {0:'Xen-like',2:'Ect-like',1:'mESC-like'}


# In[150]:


map_assignments = [k_means_dict[i] for i in assignments]


# In[152]:


sub.obs['k_means'] = map_assignments
sub


# In[155]:


sub.write_h5ad('semrau_filtered_kmeans.h5ad')


# In[ ]:





# In[ ]:




