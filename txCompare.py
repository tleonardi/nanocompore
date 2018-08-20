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


class txCompare(object):
    def __init__(self, data):
        # data [0]: Tx name
        # data [1]: [0]: data_f1 [1]: data_f2
        self.significant=[]
        (self.tx_name, rest) = data[0].split(sep="##")
        (self.tx_id, rest) = rest.split(sep="::")
        (self.chrom, rest) = rest.split(sep=":")
        self.start = rest.split(sep="-")[0]

        data_f1 = defaultdict(list)
        data_f2 = defaultdict(list)
        # Make data dictionaries:
        # Position: [read, intensity, dwell]
        for i in data[1][0]:
            if i[2][2] == 'A':
                data_f1[i[1]].append([i[3], float(i[6]), float(i[8])])
         
        for i in data[1][1]:
            if i[2][2] == 'A':
                data_f2[i[1]].append([i[3], float(i[6]), float(i[8])])
        
        self.data_f1=data_f1
        self.data_f2=data_f2

        # Save all positions observed
        self.positions = set.intersection( set(data_f1.keys()), set(data_f2.keys()) )
        self.collapsed_data = dict()
        
        for k in self.positions:
            self.collapsed_data[k] = [ ["F1"] + i for i in self.collapseObs(data_f1[k]) ] + [ ["F2"] + i for i in self.collapseObs(data_f2[k]) ]
            
        self.clusters = dict()
        pvals = dict()
        for p in self.positions:
            self.clusters[p] = self.cluster( [ (i[1], i[2]) for i in self.collapsed_data[p] ] )
            file_labels =np.transpose(self.collapsed_data[p])[0]
            f1_counts = Counter(self.clusters[p][file_labels=="F1"])
            f2_counts = Counter(self.clusters[p][file_labels=="F2"])
            f_obs = np.array([[f1_counts[0],f1_counts[1]],[f2_counts[0],f2_counts[1]]], dtype="int64")
            try:
                chi = stats.chi2_contingency(f_obs)
                pvals[p] = chi[1]
            except:
                pvals[p]="1"

        self.p = np.transpose( [ [k,p] for (k,p) in pvals.items() ] )
        adj = multipletests(np.array(self.p[1], dtype="float"), method="fdr_bh")[1]
        self.adj_pvals = { i[0]:i[1] for i in np.transpose([self.p[0], adj]) }
        for (k,v) in self.adj_pvals.items():
            if float(v)<0.1:
                kmer=''.join(set([i[2] for i in data[1][0]  if i[1] == k ]))
                self.significant.append(["chr"+self.chrom, int(self.start)+int(k), int(self.start)+int(k)+5, self.tx_id, self.tx_name,kmer,k,v])


    @staticmethod
    def collapseObs(x):
        results=[]
        reads=set([i[0] for i in x])
        for r in reads:
            q = np.transpose([ (i[1], i[2]) for i in x if i[0]==r ])
            wm = np.average(q[0], weights=q[1])
            results.append([wm, sum(q[1])])
        return results
    
    @staticmethod
    def cluster(x):
        X = StandardScaler().fit_transform(x)
        y_pred = KMeans(n_clusters=2, random_state=146).fit_predict(X)
        return y_pred
    
    def plotData(self, pos):
        dat=np.transpose(self.collapsed_data[pos])
        x=np.array(dat[1], dtype="float")
        y=np.array(dat[2], dtype="float")
        color_dict = { 'F1':'red', 'F2':'green', '0': "blue", '1': "red"}
        f, (ax1, ax2) = plt.subplots(1,2, figsize=(10,5))
        ax1.scatter(x,y, c=[ color_dict[i] for i in dat[0] ], alpha=0.5, s=0.5)
        ax1.set_xlabel("Mean intensity")
        ax1.set_ylabel("Dwell time")
        ax2.scatter(x,y, c=[ color_dict[str(i)] for i in self.clusters[pos]], alpha=0.5, s=0.5)
        ax2.set_xlabel("Mean intensity")
        ax2.set_ylabel("Dwell time")
        categories=["Cluster 1", "Cluster 2"]
        #f1_counter=Counter(clusters[pos][d[0]=="F1"])
        #f2_counter=Counter(clusters[pos][d[0]=="F2"])
        #p1 = plt.bar(categories, f1_counter.values(), 0.55, color='#d62728')
        #p2 = plt.bar(categories, f2_counter.values(), 0.55, bottom=f1_counter.values())
        #plt.legend((p2[0], p1[0]), ('Male', 'Female'))
