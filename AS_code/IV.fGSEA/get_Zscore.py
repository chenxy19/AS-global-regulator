
###### For GSEA with shuffling ######


import argparse
import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np
import gc
import random
import time

def normal_distribution(x, mu, sigma):
    """
    return the value in a normal distribution.
    Params:
    - x
    - mu: mean of distribution.
    - sigma: standard deviation.
    """
    return np.exp( -1 * ( (x-mu) ** 2) / ( 2 * (sigma ** 2)) ) / (math.sqrt( 2 * np.pi ) * sigma)

def zc_sort(a,b):  
    """
    Insertion sort. Find rank of elements in a in a large sorted list b.
    
    Params:
    - a: smaller list.
    - b: larger sorted list.
    """
    
    ind_list = []
    a_sort = sorted(a)
    i = 0
    j = 0
    while i < len(a):
        temp = a_sort[i]
        while j < len(b):
            if temp < b[j]:
                ind_list.append(j-1)
                i = i + 1
                break
            elif temp >= b[-1]:
                ind_list.append(len(b)-2)
                i = len(a)
            j = j + 1

    return(ind_list)


class GetZScore(object):
    """
    A class representation to compute the Z score of correlation for each gene.
    """
    def __init__(self, with0, nbin, w, nbackground, nthres):
        """
        Params:
        - with0: 1 indicate the correlation is taken with 0s; 0 means without 0s.
        - nbin: number of bins for the genes.
        - w: control variation of background distribution. 
        - nbackground: number of background genes. 
        - nthres: minimal number of cells required for computing correlation. 
        """
        self.expr_psi = []
        self.with0 = with0
        self.corr = []
        self.expr_level = 0
        self.expr_drate = 0
        self.nbin = nbin
        self.bin_gap = 0
        self.w = w
        self.nthres = nthres
        self.background_prob = []
        self.nbackground = nbackground
        self.background_prob_norm = []
        self.gene_bin = []
        
    def load_file(self, prefix):
        """
        load the file containing expression level and mean psi (on last row) for all genes. Ncell * (N gene+1) Dataframe.
        
        Params:
        - filename: name of the matrix.
        """
        self.prefix = prefix
        print(f"{self.prefix}")
        self.expr_psi = pd.read_csv(self.prefix, sep=" ", index_col=0)
        
        self.expr_drate = [sum(self.expr_psi.iloc[:,i] > 0)/len(self.expr_psi.index) for i in range(len(self.expr_psi.columns)-1)]
        # only look at genes above detection threshold and take the last col as well.
        which_enough = np.sum(self.expr_psi!=0, axis = 0) >= self.nthres
        self.expr_psi = self.expr_psi.loc[:,which_enough]
         
        # update for new matrix and then binning
        self.expr_drate = [sum(self.expr_psi.iloc[:,i] > 0)/len(self.expr_psi.index) for i in range(len(self.expr_psi.columns)-1)]
        self.bin_gap = (max(self.expr_drate)-min(self.expr_drate)) / self.nbin
        self.expr_level = np.sum(self.expr_psi.iloc[:,:-1], axis=0)
        

    def corr_with0(self, gene, corr_method):
        """
        compute corr for all sample between a single gene and mean psi.
        Params:
        - gene: name of gene to compute.
        - corr_method: "pearson" or "spearman"

        """
        return self.expr_psi[gene].corr(self.expr_psi["bi_glbl_mean_psi"], method=corr_method)

    def corr_non0(self, gene, corr_method):
        """
        compute corr for non-0 sample between a single gene and mean psi.
        Params:
        - gene: name of gene to compute.
        - corr_method: "pearson" or "spearman"

        """

        return self.expr_psi[gene][self.expr_psi[gene]!=0].corr(self.expr_psi["bi_glbl_mean_psi"][self.expr_psi[gene]!=0], method=corr_method)

    def take_corr(self, corr_method):
        """
        compute the correlation between all genes and mean psi (self.corr) in 2 ways.
        
        with0: pearson/spearman corr for all samples.
        without0: Take correlation after deleting all cells with 0 expression level. 
        
        Only compute for cells above the minimal detection rate, otherwise the correlation is 0-padded.
        The minimal detection rate is determined by the bin size. (genes in the first bin will be discarded in the following calculation)

        Params:
        - corr_method: "pearson" or "spearman"
        """
        
        if self.with0 == 1:
            self.corr = [self.corr_with0(gene, corr_method) for gene in self.expr_psi.columns[:-1]]
        else: 
            self.corr = [self.corr_non0(gene, corr_method) for gene in self.expr_psi.columns[:-1]]
        self.corr=pd.DataFrame(self.corr)
        self.corr.index=self.expr_psi.columns[:-1]
        self.corr.columns=["observed"]
        self.generate_background("observed")
        
    def take_corr_rand(self, corr_method, it):
        """
        compute the correlation between all genes and random generated array.
        
        params:
        - corr_method: method of taking correlation
        - it: number of iteration
        """
        
        random.shuffle(self.expr_psi["bi_glbl_mean_psi"])
        
        if self.with0 == 1:
            self.corr[f"iter{it}"] = [self.corr_with0(gene, corr_method) for gene in self.expr_psi.columns[:-1]]
        else: 
            self.corr[f"iter{it}"] = [self.corr_non0(gene, corr_method) for gene in self.expr_psi.columns[:-1]]
        self.generate_background(f"iter{it}")
        
    def bin_genes(self):
        """
        1.Binning genes with 0/1 mapping matrix by evenly dividing their detection rate (gene_bin).
        2.Compute probability of sampling a gene from bin i for a gene in bin j (background_prob).
        3.Mapping bin-wise probability to gene-wise probability (background_prob_gene).
        4.Take CDF.
        5.Normalization to 1.

        """
        bin_start = [(min(self.expr_drate) + self.bin_gap * i) for i in range(self.nbin)]
        # to include the genes with max expr_drate
        bin_start.append(max(self.expr_drate) + self.bin_gap)
        self.gene_bin = [[int(j>=bin_start[i] and j<bin_start[i+1]) for i in range(self.nbin)] for j in self.expr_drate]
        
        # make sure every gene has a bin
        if sum([int(sum(self.gene_bin[i])==0) for i in range(len(self.gene_bin))])==0:
            print("Every gene has a bin")
        else:
            print("Error: Not every gene has a bin!")

        no_gene_per_bin = np.sum(self.gene_bin, axis=0)
        self.background_prob = np.array([[normal_distribution(i - j, 0, self.w)/no_gene_per_bin[ind] for ind,i in enumerate(bin_start[:-1])] for j in bin_start[:-1]])
        background_prob_gene = np.matmul(np.matmul(np.array(self.gene_bin) , np.array(self.background_prob)), np.array(self.gene_bin).T)
        
        # first make CDF
        background_prob_CDF = np.zeros(len(background_prob_gene))
        background_prob_CDF = [background_prob_CDF]
        non_print = [background_prob_CDF.append(background_prob_CDF[-1]+background_prob_gene[j]) for j in range(len(background_prob_gene))]
        background_prob_CDF = np.array(background_prob_CDF)
        del background_prob_gene, non_print
        
        # normalize per row
        scale_factor = background_prob_CDF[-1]
        background_prob_CDF = background_prob_CDF.T
        self.background_prob_norm = background_prob_CDF/scale_factor.reshape(-1,1)
        del background_prob_CDF
        
    def generate_background(self, col_name):
        """
        Generate background genes and compute Z score.
        
        Params:
        - colname: the column name for the correlation to convert. Overwrite the correlation with zscore.
        """
        # random number in [0,1)
        gene_rand = np.random.random_sample((len(self.background_prob_norm), self.nbackground))
        self.back_gene_ind = []
        [self.back_gene_ind.append(zc_sort(gene_rand[j],self.background_prob_norm[j])) for j in range(len(self.background_prob_norm))]
       
        corr2 = self.corr[col_name]
        corr_dict = dict(zip(range(len(corr2)), corr2))
        background_corr = [[corr_dict.get(j) for j in gen_list] for gen_list in self.back_gene_ind]
        background_mean = np.mean(background_corr, axis=1)
        background_std = np.std(background_corr, axis=1)
        del background_corr
        
        z_score = (corr2 - background_mean)/background_std
        self.corr[col_name] = z_score
        
    def save_file(self, output_file):
        self.corr.to_csv(output_file)

        
def main():
    my_parser = argparse.ArgumentParser(description='Specify the expression_psi file directory, and output file.')
    # mandatory variable
    my_parser.add_argument('expr_psi', action='store', help="expr_psi")
    my_parser.add_argument('output_file', action='store', help="output_file")
    my_parser.add_argument('n_iteration',  type = int, action='store', help="number of iteration")

    args = my_parser.parse_args()

    expr_psi = args.expr_psi
    output_file = args.output_file
    iterations = args.n_iteration
    
    #expr_psi_dir = "/Users/chenxinyi/Desktop/" 
    #output_file = "Zscore_met2.csv"
    #iterations=5

    # whether to take corr with zeros , number of bins for the genes, variation of background distribution (0.1 originally), number of background genes, number of non-zero samples required for taking correlation.
    get_zscore = GetZScore(0, 10, 0.1, 100, 6)  
    get_zscore.load_file(expr_psi)
    get_zscore.bin_genes()
    get_zscore.take_corr("spearman")

    start=time.time()
    for it in range(iterations):  # can be quite slow
        get_zscore.take_corr_rand("spearman", it)  # corr method 
    get_zscore.save_file(output_file)
    print(f"{iterations} iter time: {time.time()-start}")

if __name__ == '__main__':
    main()