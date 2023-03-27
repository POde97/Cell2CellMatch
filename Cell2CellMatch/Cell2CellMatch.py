from CellIDpy import ProteinCoding
import scanpy as sc
import networkx as nx
from CellIDpy import CellID
from CellIDpy import HyperG
import scipy
import pandas as pd 
import numpy as np
import igraph as ig
import leidenalg as la
import tqdm
from Cell2CellMatch import percolation 
from Cell2CellMatch import plot
import torch
from scipy import sparse
from scipy.stats import fisher_exact
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

#Cell2CellMatch build functonal network of cells leveraging CellID multiple correspondence analysis based representation.
#Arguments
#adata: AnnData object (expample scRna-seq)
#batch_key = key in adata metadata that identify different patients
#resolution = resolution parameter for community detection 
#n_clusters = if int (example: 9) perform community detection and find best resolution parameter to obtain 9 communities
#ref: Reference gene list: if we want to label our data as in CellID (See https://github.com/POde97/Cell-IDpy)
#gpu_cuda: if True use gpu for speed up calculation
#strategy_T: ["conservative","percolation"] 
             # if conservative -> increase tresholds until pruned edge do not disconnect only sigleton partion, use the last
             #treshold best treshold 
             # if percolation-> find best treshold using a percolation based methods. Best treshold is selected right before 
             #percolation treshold
#custom_T: float with custom treshold 
#strategy_G: ["local-topology","K-partite"] 
            #if local topology: build a network in which all cells were test against all cells through hypergeometric test 
            #if K-partite: Hypergeometrically test cells that belong to different dataset in a pairwise manner and build a K-partite network
#mply_p: max_chunks_per_worker of mapply parallelization https://pypi.org/project/mapply/


class Cell2CellMatch():

  def __init__(self,adata,batch_key="batch",resolution =1.0,n_clusters =None,ref = None,gpu_cuda =None,strategy_T ="conservative",custom_T=None,strategy_G ="local-topology",mply_p = 16):


    self.batch_key = batch_key
    self.n_clusters = n_clusters

    self.adata = adata
    adlist = [self.adata[self.adata.obs[batch_key]==i].copy() for i in list(self.adata.obs[batch_key].unique())]
    adlist = [self.preproc(x=x) for x in adlist]
    
    
    #run CELLID per batch + cell labelling if reference available
    #Cell-ID
    print("Performing Cell-ID per batch")
    self.adlist = [CellID.CellID(ad=x,verbose =False).ad for x in tqdm.tqdm(adlist)]
    #Labelling
    if ref is not None:
      print("Per-batch labelling with reference")
      for x in tqdm.tqdm(self.adlist): 
        HyperG.Hypergeom(x,ref,verbose = False,mply_p=mply_p)
      self.adata.obs = pd.concat([self.adata.obs,pd.concat([xf.obs["gs_prediction"] for xf in self.adlist])],axis=1)
    signature_list=[pd.DataFrame(self.adlist[i].obs["signature"]) for i in range(len(self.adlist))]
    signature_df =pd.concat(signature_list,axis=0) 
    self.adata.obs = pd.concat([self.adata.obs,signature_df],axis=1)
    

    
    #Hyper New -> Not Kpartite
    if strategy_G == "local-topology":
      print("Hypergeometric Test")
      H = HyperG.Hypergeom(self.adata,self.adata,prediction = False,verbose = True,gpu_cuda = gpu_cuda,mply_p=mply_p)
    

      df = H.Intersection.copy()
  

      maxx = df.values.flatten()[df.values.flatten()!=np.inf].max() + 1 
      df.values[df.values == np.inf] = maxx

    #G = nx.from_numpy_array(df.to_numpy())(Networkx 3.0)
      #G = nx.from_numpy_matrix(df.to_numpy())
      Adj = sparse.csr_matrix(df.values)
      G = nx.from_scipy_sparse_matrix(Adj)
      to_remove = [(a,b) for a, b, attrs in G.edges(data=True) if attrs["weight"] < 2.0]
      G.remove_edges_from(to_remove)
      names = list(df.index)
      mapping = dict(zip(G, names))
      G = nx.relabel_nodes(G, mapping)
      T=2.0

    #Strategy for build the network :
      print("Build network + sparsification") 
      if strategy_T == "conservative":
        t = [len(c) for c in sorted(nx.connected_components(G),key = len,reverse = True)]
        while  len([i for i in t if i>1]) == 1:
          T = T + 1.0
                        
          to_remove = [(a,b) for a, b, attrs in G.edges(data=True) if attrs["weight"] < T ]
          G.remove_edges_from(to_remove)
          t = [len(c) for c in sorted(nx.connected_components(G),key = len,reverse = True)]
    #
      if strategy_T == "percolation":
        T = percolation.find_threshold_bfs(np.array(nx.adjacency_matrix(G).todense()))

      if custom_T is not None:
        T = custom_T         

    
             
      T = T-2.0
      if T < 2.0:
        T = 2.0
      #print(T)
      self.T = T

      df = H.Intersection.copy()
      maxx = df.values.flatten()[df.values.flatten()!=np.inf].max() + 1 
      df.values[df.values == np.inf] = maxx
      self.df_ceck = df
      #G = nx.from_numpy_array(df.to_numpy())
      Adj = sparse.csr_matrix(df.values)
      G = nx.from_scipy_sparse_matrix(Adj)
      to_remove = [(a,b) for a, b, attrs in G.edges(data=True) if attrs["weight"] < T]
      G.remove_edges_from(to_remove)
      names = list(df.index)
      mapping = dict(zip(G, names))
      G = nx.relabel_nodes(G, mapping)
    
      #Select biggest connected community
      Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
      G0 = G.subgraph(Gcc[0])
      G = G0.copy()
      del G0

      self.G = G
    
    if strategy_G =="K-partite":
      #Hyper Old -> K-partite
      #Hyper for all couples (Old)
      self.PvalMS = self.RunHGT()
      G = self.KpG(OptT=True)
      self.G = G[0]
      



    #print("ok4")("Old")
    #self.G = BuildG[0]
    #self.T = BuildG[1]
    
    #Properly install natworkx with correct scipy dep to avoid error
    #Da Rimettere!!
    #self.Adj = nx.adjacency_matrix(self.G)

    #print("ok1")
    #Community Detection on G 
    
    #A) Modularity communities
    print("Find communities")
    if n_clusters == None:
      self.community = self.Leiden(res = resolution) 
      adata.obs = pd.concat([adata.obs,self.community],axis=1)
      self.CS = self.XClusterSignature()
    
    #B) Find proper Resolution for K_communities
    else:
      max_steps = 20 
      this_step = 0
      this_min = float(0)
      this_max = float(3)
      while this_step < max_steps:
        this_resolution = this_min + ((this_max-this_min)/2) 
        
        
        self.community = self.Leiden(res = this_resolution)
        
        n_C = self.community['C2C-cl'].nunique()
        
        

        if n_C > n_clusters:
          this_max = this_resolution
        elif n_C < n_clusters:
          this_min = this_resolution
      
        this_step +=1

      adata.obs = pd.concat([adata.obs,self.community],axis=1)
      #self.adata.obs["C2C-cl"] = self.adata.obs["C2C-cl"].fillna("unassigned")
      self.CS = self.XClusterSignature()
      
    

    
    

  def RunHGT(self):#PvalueMatrixs
    print("Run Hypergeometric Tests")
    listSignature1 = self.adlist.copy()
    listSignature2 = self.adlist.copy()

    CombinationAdata1=[]
    CombinationAdata2=[]
    self.PvalMS = []
    for x1 in tqdm.tqdm(listSignature1):
      listSignature2 = listSignature2[1:]
      #print(len(listSignature2))
      if(len(listSignature2)>0):
        for x2 in listSignature2:
    #CombinationAdata1.append(x1)
    #CombinationAdata2.append(x2)
          self.PvalMS.append(HyperG.Hypergeom(x1,x2,prediction = False,verbose = False).Intersection)

    return self.PvalMS


  

  def KpG(self,OptT):

    def Optimazed_Build_G_with_Treshold(matrix1,Treshold):
      G=nx.algorithms.bipartite.matrix.from_biadjacency_matrix(scipy.sparse.csr_matrix(matrix1))
      to_remove = [(a,b) for a, b, attrs in G.edges(data=True) if attrs["weight"] < Treshold]
      G.remove_edges_from(to_remove)
    
      return G

    def merge_two_dicts(x, y):
      z = x.copy()   # start with keys and values of x
      z.update(y)    # modifies z with keys and values of y
      return z


    list_of_sample_Seurat1 = self.PvalMS.copy()
    list_of_sample_Seurat2 = self.PvalMS.copy()
    
    T = 2.0
    G = nx.Graph()
    idx1=[i for i in range(len(list_of_sample_Seurat1))]
    idx2 = [i for i in range(len(list_of_sample_Seurat1))]

    for i in range(len(list_of_sample_Seurat1)):#
        
      matrix1 = list_of_sample_Seurat1[i]
      maxx = matrix1.values.flatten()[matrix1.values.flatten()!=np.inf].max() + 1 
      #print(maxx)
      matrix1.values[matrix1.values == np.inf] = maxx
      #matrix1 = matrix1.T
      #Rename
      T = 2.0
      G_temp =Optimazed_Build_G_with_Treshold(matrix1,T)
      sample2_name = list(matrix1.index)
      sample1_name= list(matrix1.columns)
      names = sample2_name+ sample1_name
      #print(len(names))
      mapping = dict(zip(G_temp, names))
      G_temp = nx.relabel_nodes(G_temp, mapping)
                
                
                
      if OptT ==True:
        t = [len(c) for c in sorted(nx.connected_components(G_temp),key = len,reverse = True)]
        while  len([i for i in t if i>1]) == 1:
          T = T + 1.0
                        
          to_remove = [(a,b) for a, b, attrs in G_temp.edges(data=True) if attrs["weight"] < T ]
          G_temp.remove_edges_from(to_remove)
          t = [len(c) for c in sorted(nx.connected_components(G_temp),key = len,reverse = True)]
                
                
      T = T-1.0
      if T < 2.0:
        T = 2.0
      G_temp = Optimazed_Build_G_with_Treshold(matrix1,T)
      G_temp = nx.relabel_nodes(G_temp, mapping)
      largest_cc = max(nx.connected_components(G_temp),key = len)
      Go = G_temp.subgraph(largest_cc)
      G_temp = Go
      G.add_nodes_from(G_temp.nodes(data=True))
        
      G.add_edges_from(G_temp.edges(data = True))

    return G,T


  def Leiden(self,res):
    h = ig.Graph.from_networkx(self.G)
    
    if self.n_clusters == None:
      part = la.find_partition(h, la.ModularityVertexPartition,weights="weight",n_iterations=-1)
    else:

      part = la.RBConfigurationVertexPartition(h,weights="weight")
      optimiser = la.Optimiser()
      diff = optimiser.optimise_partition(part,n_iterations =-1)
    
    d = {'node': h.vs['_nx_name'],'C2C-cl' : part.membership}
    
    return pd.DataFrame.from_dict(d).set_index("node")                 


  def XClusterSignature(self):
    df = pd.DataFrame()
    for i in self.adlist:
      df_temp = pd.DataFrame(pd.DataFrame(i.obs["signature"].tolist(),index=list(i.obs.index)))
      df = pd.concat([df,df_temp])
    geneCl = pd.concat([self.community,df],axis=1).dropna()

    perC = []
    dftot = pd.DataFrame()
    for clu in list(geneCl["C2C-cl"].dropna().unique()):
  
      geneSingC = pd.value_counts(geneCl[geneCl["C2C-cl"]==clu].drop("C2C-cl",axis=1).to_numpy().flatten()).iloc[:200]
      dftemp = pd.DataFrame(geneSingC,columns = ["cluster-"+str(clu)]).T.reset_index()
      dftot = pd.concat([dftot,dftemp])
        
    dftot = dftot.reset_index().drop("level_0",axis=1)
    signature = pd.DataFrame(dftot.apply(lambda x: list(x.dropna().index)[1:],axis=1),columns=["signature"])
    signature.index = dftot["index"]
    return signature
  
  def XclusterTotFreq(self):
    df = pd.DataFrame()
    for i in self.adlist:
      df_temp = pd.DataFrame(pd.DataFrame(i.obs["signature"].tolist(),index=list(i.obs.index)))
      df = pd.concat([df,df_temp])
    geneCl = pd.concat([self.community,df],axis=1).dropna()

    #Probability of a gene in cluster
    perC = []
    dftot = pd.DataFrame()
    for clu in list(geneCl["C2C-cl"].dropna().unique()):
  
      geneSingC = pd.value_counts(geneCl[geneCl["C2C-cl"]==clu].drop("C2C-cl",axis=1).to_numpy().flatten())
      dftemp = pd.DataFrame(geneSingC,columns = ["cluster-"+str(clu)]).T.reset_index()
      dftot = pd.concat([dftot,dftemp])
    dftot = dftot.reset_index().drop("level_0",axis=1)
    dftot = dftot.fillna(0.)
    return dftot

  def preproc(self,x):
    x = x[:,x.var_names.isin(ProteinCoding.HgProteinCoding())].copy()
    #x.raw = x
    sc.pp.normalize_total(x, target_sum=1e4)
    sc.pp.log1p(x)
    #sc.pp.scale(x)

    return x


  def CommonGenes(self):
    df = pd.DataFrame(self.G.edges(),columns=["i","j"])
    df["signature-i"] = df["i"]
    df["signature-j"] = df["j"]
    m_dict = self.adata.obs["signature"].to_dict()
    df["signature-j"] = df["signature-j"].map(m_dict)
    df["signature-i"] = df["signature-i"].map(m_dict)
    df["Intersection-genes"] = df.mapply(lambda x: intersection(x["signature-i"],x["signature-j"]),axis=1)

    return df
  
  def HyperGxClSignature(self,ref):
    df = self.XClusterSignature()
    adata_cl = sc.AnnData(np.zeros((len(df),len(self.adata.var_names))))
    adata_cl.var_names = list(self.adata.var_names)

    adata_cl.obs_names =  list(df.index)
    adata_cl.obs = pd.concat([adata_cl.obs,df],axis=1)
    HGT_cl= HyperG.Hypergeom(adata_cl,ref)
  
    return adata_cl.obs
    
    
  def fisher_exact_test(self,df):
    """
    Perform Fisher exact test on a dataframe and return the p-values.
    
    Args:
    - df: a pandas DataFrame where rows are genes and columns are cluster elements and ij is the counts of the genes in the cluster.
    
    Returns:
    - A pandas DataFrame with p-values for each gene and cluster.
    """
    
    #print(len(df.columns))
    p_values = pd.DataFrame(index=df.index, columns=df.columns)
    for col in tqdm.tqdm(df.columns):
      #print(col)
      n = df[col].sum()
      for gene in df.index:
        a = df.loc[gene, col]
        b = n - a
        c = df.loc[gene].sum() - a
        d = df.drop(col, axis=1).sum().sum() - c
        _, p = fisher_exact([[a, b], [c, d]], alternative='greater')
        p_values.loc[gene, col] = p
    return p_values

  def most_specific_genes(self, n):
    """
    Estimate the n most specific genes for each cluster based on the p-values from the Fisher exact test.
    
    Args:
    - df: a pandas DataFrame where rows are genes and columns are cluster elements and ij is the counts of the genes in the cluster.
    - n: the number of genes to return for each cluster.
    
    Returns:
    - A pandas DataFrame with the n most specific genes for each cluster, along with their p-values.
    """
    df = self.XclusterTotFreq()
    cl_namez = df["index"]
    df = df.drop("index",axis=1)
    df.index = list(cl_namez)
    df = df.T
    
    p_values = self.fisher_exact_test(df)
    
    results = pd.DataFrame(index=range(n), columns=df.columns)
    for col in df.columns:
      sorted_genes = p_values[col].sort_values().index[:n]
      results[col] = pd.Series(sorted_genes, index=range(n))
      results[col+'_pval'] = p_values.loc[sorted_genes, col].values
    return results



  #VISUALIZATION#######################
  def Clvis1(self,label,n=0,color=None):
    dftot = self.adata.obs.loc[~self.adata.obs["C2C-cl"].isna()]
    unique= dftot["C2C-cl"].unique().astype(int)
    CC = []
    for i in dftot[self.batch_key].unique():
      dftemp = dftot.loc[dftot[self.batch_key]==i]
      C_t= dftemp[[label,"C2C-cl"]].groupby([label,"C2C-cl"]).size().unstack(fill_value=0).T
      C_t = C_t/len(dftemp)

  #unique = list(set(list(df1.index) + list(df2.index)))
    #unique.remove('unasigned')
      C_t=pd.concat([pd.DataFrame(np.zeros(len(unique)),index = unique),C_t],axis=1).fillna(0).drop(0,axis=1)


      CC.append(C_t)

#from plot import*
    if n is not None:
      dg = dftot["C2C-cl"].value_counts()
      CC = [i[i.index.isin(list(dg[dg>n].index))] for i in CC]
    plt.figure(figsize=(24,8))
    plot.plot_clustered_stacked(CC,[1.01,0.1],[1.41,0.1],list(dftot[self.batch_key].unique()),title="Local-Topology",colours=color)
    
    
    

