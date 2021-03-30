import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import random
import re
import decimal
from decimal import Decimal, getcontext
import argparse

class Gene_Set(object):
    """
    A class representation of the KEGG gene set, Stores the name and gene contents of KEGG gene sets.
    """
    def __init__(self):
        """
        Initialize gene set class. 
        KEGG_annotation is [char],gene_list is [[char],[char],...], dictionary is {KEGG:gene_name}
        """
        self.KEGG_annotation=[] 
        self.gene_list=[]
        self.dictionary={}
        
    def add_gene_set(self,KEGG,genes):
        """
        add a pair of new KEGG annotation and gene list into the class.
        
        :param KEGG: a char starting with "KEGG_"
        :param genes: a list containing all genes in the set.
        """
        self.KEGG_annotation.append(KEGG)
        self.gene_list.append(genes)

    def make_dict(self):
        """
        make a dictionary mapping KEGG name to its gene list.
        """
        self.dictionary = dict(zip(self.KEGG_annotation, self.gene_list))
 
        
    def get_gene_set(self,KEGG):
        """
        Return all genes in the KEGG gene set specified.
        
        :param KEGG: a char starting with "KEGG_" specifying the gene set to retrieve.
        return: a list of genes in the set.
        """
        genes_to_report = self.dictionary.get(KEGG)
        return(genes_to_report)
    
    def get_mice_gene_set(self,KEGG):
        """
        Return all mice genes in the KEGG gene set specified.
        
        :param KEGG: a char starting with "KEGG_" specifying the gene set to retrieve.
        return: a list of genes in the set, converted to mice gene format.
        """
        genes_to_report = self.dictionary.get(KEGG)
        mice_genes = [gene[0].upper()+gene[1:].lower() for gene in genes_to_report]
        #print(f"mice genes:{mice_genes}")
        return(mice_genes)
    
    def drop_gene_set(self, if_save):
        """
        After certain numbers of iterations, check for gene sets that do not meet criteria, and omit them in further iterations.
        
        :param if_save: a list of True/False of same length with the number of KEGG. If false, drop from the list.
        """
        
        if_save = list(if_save[0])
        self.KEGG_annotation = [self.KEGG_annotation[i] for i in if_save]
        self.gene_list = [self.gene_list[i] for i in if_save]
        self.dictionary = dict(zip(self.KEGG_annotation, self.gene_list))
        

class GSEA(object):
    """
    A class to implement a given gene set enrichment analysis.
    """
    def __init__(self, species, iterations):
        """
        Initialize GSEA class. Specify the species and iterations to run. 
        
        """
        self.gene_sorted=[]  
        self.gene_set_list=[] # a Gene_Set instance
        self.species = species
        self.iterations = iterations
        
    def load_file(self, expfile, genesets,inverse):
        """
        Load data from the expression file, true label file and gene set file.
       
        :param expfile: the path to preprocessed expression matrix. row: genes, col: samples
        :param sampfile: the path to true labels of the samples. rownames are samples and 0/1 denotes classes.
        :param genesets: the path to the gene set data. 
        self.expression and self.sample_assignment are pandas DataFrames. self.gene_set_list is a Gene_Set instance.
        """
 
        rank=pd.read_table(expfile,sep=" ",index_col=0)
        self.gene_sorted = list(rank.iloc[:,0])
        if inverse==1:
            self.gene_sorted=self.gene_sorted[::-1]
            
        
        with open(genesets,"r") as f:
            KEGGs=f.readlines()
        gene_set = Gene_Set()
        for line in KEGGs:
            KEGG_gene = line.strip("\n").split("\t")
            desired_pathway = len(KEGG_gene[2:]) <= 200 and len(KEGG_gene[2:]) >= 30
            if desired_pathway:
                gene_set.add_gene_set(KEGG_gene[0],KEGG_gene[2:])
        self.gene_set_list = gene_set
        self.gene_set_list.make_dict()
        self.n_bonfer=len(self.gene_set_list.KEGG_annotation)
        print(f"KEGG set number:{len(self.gene_set_list.KEGG_annotation)}")
        
    def get_enrichment_score(self, geneset):
        """
        Get the enrichment score for a specified geneset under a specific gene rank.
       
        :param genesets: a given gene set like ‘KEGG_CITRATE_CYCLE_TCA_CYCLE’.
        return: the enrichment score, a float correct to two decimal places, for a given gene set.
        """
        #genes_sorted_FC = self.get_gene_rank_order() # genes_sorted_FC is all genes ranked by logFC.
        
        N_total = len(self.gene_sorted)  # total number of genes

        #Get ES for one KEGG
        if self.species=="human":
            gene_list=self.gene_set_list.get_gene_set(geneset)
            #print("using human gene sets")
        elif self.species=="mice":
            #print("using mice gene set")
            gene_list=self.gene_set_list.get_mice_gene_set(geneset)  # gene_list is a list of genes from the specified gene set.
        else:
            print("specify species for the GSEA class, either human or mice")

        N_gene_set= len(set(gene_list) & set(self.gene_sorted)) # numbers of genes both in the gene set and in the expression matrix
        try:
            P_hit = math.sqrt((N_total - N_gene_set) / N_gene_set)
            P_miss = (-1)*math.sqrt(N_gene_set / (N_total - N_gene_set))

            sum_score=[0]
            for gene in self.gene_sorted:
                if gene in gene_list:
                    new_score = sum_score[-1] + P_hit
                    sum_score.append(new_score)
                else:
                    new_score = sum_score[-1] + P_miss
                    sum_score.append(new_score)

            enrichment_score = max(sum_score)
        except:
            enrichment_score=0

            #plt.plot(range(len(sum_score)),sum_score)
            #plt.show()

            #with open("kegg_enrichment_scores.txt","a") as f:
            #    f.write(f"{geneset} {enrichment_score}\n")
        return(enrichment_score)
        
    def get_sig_sets(self,p,output_file):
        """
        Get the list of significant gene sets.
       
        :param p: a given threshold for p value (before Bonferonni correction).
        return: the list of significant gene sets (as strings), at a corrected threshold of p, by name. 
        """
        
        # a list of real ES scores without permutation
        real_ES_scores = np.array([self.get_enrichment_score(KEGG) for KEGG in self.gene_set_list.KEGG_annotation])

        #print(f"real_ES_score_computed:{time.time() - start}")
        
        # permutation
        #permutated_ES_scores = []
        p_real_number = np.zeros(len(self.gene_set_list.KEGG_annotation))
        total_iter = self.iterations
        
        #start=time.time()
        
        for n_iter in range(total_iter):  

            random.shuffle(self.gene_sorted)
            #permutated_ES_scores.append([self.get_enrichment_score(KEGG) for KEGG in self.gene_set_list.KEGG_annotation])
            permutated_ES_scores = [self.get_enrichment_score(KEGG) for KEGG in self.gene_set_list.KEGG_annotation]
            is_bigger = [int(permutated_ES_scores[i] >= real_ES_scores[i]) for i in range(len(self.gene_set_list.KEGG_annotation))]
            
            #print(f"one iter:{time.time() - start}")
            #start=time.time()
            
            p_real_number = p_real_number + is_bigger

            if n_iter in [10,20,50,100,200,500,1000,10000]:
                if_drop = np.where(p_real_number >= 5)
                if_save = np.where(p_real_number < 5)

                p_real_drop = p_real_number[if_drop]  
                real_ES_drop = real_ES_scores[if_drop]
                pd.DataFrame(p_real_drop/n_iter, index = np.array(self.gene_set_list.KEGG_annotation)[if_drop]).to_csv(f"{output_file}_last_drop.csv")

                self.gene_set_list.drop_gene_set(if_save)
                p_real_number = p_real_number[if_save]  
                real_ES_scores = real_ES_scores[if_save]
                
        try:
            #p_real = [sum(permutated_ES_scores.iloc[:,i] >= real_ES_scores[i])/total_iter for i in range(len(permutated_ES_scores.columns))]
            p_real = p_real_number / total_iter
            print(f"p_real:{p_real}")
            P_real_df = pd.DataFrame(p_real, index = self.gene_set_list.KEGG_annotation )        
            P_real_df.to_csv(f"{output_file}_{total_iter}.p.csv")

            significant_list=[]

            for i in range(len(p_real)):
                if p_real[i] < p/self.n_bonfer:
                    significant_list.append(self.gene_set_list.KEGG_annotation[i])

            with open(f"{output_file}_{total_iter}.txt","w") as f:
                f.write("\n".join(significant_list))

        except:
            print("None significant left")
    
    def implement_GSEA(self, expfile, genesets, out_sig_file, inverse=0):
        """
        Implement the GSEA pipeline.
       
        :param expfile: the path to preprocessed expression matrix. row: genes, col: samples
        :param sampfile: the path to true labels of the samples. rownames are samples and 0/1 denotes classes.
        :param genesets: the path to the gene set data. 
        return: the list of significant gene sets (as strings), at a corrected threshold of p, by name. 
        """
        significance_level = 0.05
        self.load_file(expfile, genesets, inverse) 
        self.get_sig_sets(significance_level,out_sig_file)


def main():
    my_parser = argparse.ArgumentParser(description='Specify the gene rank, gene set .gmt file, output file.')
    # mandatory variable
    my_parser.add_argument('gene_rank', action='store', help="gene_rank")
    my_parser.add_argument('gene_set', action='store', help="gene_set")
    my_parser.add_argument('output_file', action='store', help="output_file")
    my_parser.add_argument('n_iteration',  type = int, action='store', help="number of iteration")

    # Selective variable
    my_parser.add_argument('-inv', action='store', type = int, default = 0, help="1/0 take the inverse rank or not.")
    my_parser.add_argument('-sp', action='store', type = str, default = "human", help="species")
    args = my_parser.parse_args()

    gene_rank = args.gene_rank
    gene_set = args.gene_set
    output_file = args.output_file
    iterations = args.n_iteration
    inv = args.inv
    species = args.sp

    gsea = GSEA(species, iterations)
    gsea.implement_GSEA(gene_rank, gene_set , output_file, inverse=inv)



if __name__ == '__main__':
    main()
