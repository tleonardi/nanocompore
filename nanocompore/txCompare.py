import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
import numpy as np
import pickle
from collections import defaultdict
from collections import Counter
from scipy import stats
from statsmodels.stats.multitest import multipletests
from bedparse import bedline
import seaborn as sns
import pandas as pd
from collections import namedtuple


class txCompare(object):
    def __init__(self, data):

        (self.tx_name, rest) = data[0].split(sep="##")
        (self.tx_id, rest) = rest.split(sep="::")
        (self.chrom, rest) = rest.split(sep=":")
        self.start = rest.split(sep="-")[0]
        
        Events = namedtuple('Events', ['file', 'ref_name', 'ref_pos', 'ref_kmer', 'read_name', 'n_events', 'mean', 'median', 'std', 'var'])
        # Make a dictionary of Events using ref_pos as key
        self.data = defaultdict(list)
        for i in data[1][0]:
            self.data[int(i[1])].append(Events(*['F1']+i))
        for i in data[1][1]:
            self.data[int(i[1])].append(Events(*['F2']+i))
            
        # Save positions of interest
        file_one_A_positions = {k for k,v in self.data.items() for e in v if e.file=='F2' and e.ref_kmer[2]=='A'}
        file_two_A_positions = {k for k,v in self.data.items() for e in v if e.file=='F1' and e.ref_kmer[2]=='A'}
        self.positions = set.intersection( file_one_A_positions, file_two_A_positions)
        
        self.clusters = dict()
        pvals = dict()
        for p in self.positions:
            self.clusters[p] = self.cluster( [ (i.mean, i.n_events) for i in self.data[p] ] )
            file_labels = [ i.file for i in self.data[p] ]
            f1_counts = Counter(self.clusters[p][[i=="F1" for i in file_labels]])
            f2_counts = Counter(self.clusters[p][[i=="F2" for i in file_labels]])
            f_obs = np.array([[f1_counts[0],f1_counts[1]],[f2_counts[0],f2_counts[1]]], dtype="int64")
            if any([k<5 for i in f_obs for k in i ]):
                pvals[p]="1"
            else:
                try:
                    chi = stats.chi2_contingency(f_obs)
                    pvals[p] = chi[1]
                except:
                    pvals[p]="1"

        self.p = np.transpose( [ [k,p] for (k,p) in pvals.items() ] )
        adj = multipletests(np.array(self.p[1], dtype="float"), method="fdr_bh")[1]
        self.adj_pvals = { i[0]:float(i[1]) for i in np.transpose([self.p[0], adj]) }

    @staticmethod
    def cluster(x):
        X = StandardScaler().fit_transform(x)
        y_pred = KMeans(n_clusters=2, random_state=146).fit_predict(X)
        return y_pred
  

    def significant(self, thr=0.1):
        """ Return the list of significant regions. 
        If a bed file is provided also return the genome coordinates """
        results = []
        for k,v in self.adj_pvals.items():
            if v <= thr:
                results.append([self.tx_name, self.tx_id, self.chrom, self.start, int(k), v])
        return results


    def plotClust(self, pos):
        dat=np.transpose(self.data[pos])
        x=np.array(dat[6], dtype="float")
        y=np.array(dat[5], dtype="int")
        color_dict = { 'F1':'red', 'F2':'green', '0': "blue", '1': "red"}
        f, (ax1, ax2) = plt.subplots(1,2, figsize=(10,5))
        ax1.scatter(x,y, c=[ color_dict[i] for i in dat[0] ], alpha=0.5, s=0.5)
        ax1.set_xlabel("Mean intensity")
        ax1.set_ylabel("Events")
        ax2.scatter(x,y, c=[ color_dict[str(i)] for i in self.clusters[pos]], alpha=0.5, s=0.5)
        ax2.set_xlabel("Mean intensity")
        ax2.set_ylabel("Events")
        categories=["Cluster 1", "Cluster 2"]

    def plotKde(self, position=None, context=2, fig_height=5):
        if positions is None:
            return None
        interval=range(position-context+1, position+context)
        df=pd.DataFrame([i for p in interval for i in comp.data[p] ]).apply(pd.to_numeric, errors='ignore')
        df['n_events'] = np.log2(df['n_events']+0.1)
        sns.set_style("white")
        def mykdeplot(x, y, color, **kwargs):
            cmap = sns.light_palette(color, as_cmap=True)
            sns.kdeplot(x, y, alpha=1, linewidths=1, cmap=cmap, **kwargs)
        g = sns.FacetGrid(df, col="ref_pos", hue="file", height=fig_height)
        g = (g.map(mykdeplot, 'mean', 'n_events', shade_lowest=False).add_legend())
        g.set(xlabel='Mean intensity', ylabel='Dwell time')
        return g
