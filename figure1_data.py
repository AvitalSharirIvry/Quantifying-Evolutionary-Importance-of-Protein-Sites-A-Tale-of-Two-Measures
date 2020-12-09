from functions import rank_consurf_file,separate_summary_consurfDB,pearsonr_slope_dist_cons_inPDB,read_propert_over_curated_file
from scipy.stats.stats import pearsonr,spearmanr,linregress
from scipy.stats import binned_statistic
import pandas as pd

def run(chosen_file,summary_consurf_file,out):
    chosen_file=open(chosen_file)
    out=open(out,'w')
    Cax_of_pdb,seq_of_pdb=read_propert_over_curated_file('Cax_over_curated_chains')
    consurfDB_file_lines=separate_summary_consurfDB(summary_consurf_file)
    overall_dataset_ranks,overall_dataset_pearsons,overall_dataset_slopes,overall_dataset_ranked_pearsons,overall_dataset_ranked_slopes=[],[],[],[],[]
    yeast_pdb=[]
    for line in chosen_file:
        yeast,pdb,chain=line.split()[0],line.split()[1].lower(),line.split()[2]
        print yeast,pdb,chain
        if pdb+chain in consurfDB_file_lines and pdb+chain in Cax_of_pdb and yeast+pdb not in yeast_pdb:
            yeast_pdb.append(yeast+pdb)
            lines=consurfDB_file_lines[pdb+chain]
            (zipped_sorted_consurfs,zipped_ranks,ranks,consurfs,most_conserved_res,res_of_interest_most_conserved)=rank_consurf_file(lines,[])
    #Pearson correlation (conservation,dist) for each residue
    # find for that protien, the correlation between pearson and consurf of all residue
            correlations,pval_correlations,slopes,ranked_pearsons,ranked_slopes=pearsonr_slope_dist_cons_inPDB(pdb,chain,lines)
            protein_corr_cor_consurf={}
            protein_corr_slope_consurf={}
            consurfs_list=[]
            correlations_list=[]
            slopes_list=[]
            ranks_list=[]
            ranked_correlations_list=[]
            ranked_slopes_list=[]
            for res in correlations:
                consurfs_list.append(consurfs[res])
                correlations_list.append(correlations[res])
                slopes_list.append(slopes[res])
                ranks_list.append(ranks[res])
                ranked_correlations_list.append(ranked_pearsons[res])
                ranked_slopes_list.append(ranked_slopes[res])
            try:
                protein_corr_cor_consurf[pdb+chain]=pearsonr(consurfs_list,correlations_list)[0]
                protein_corr_slope_consurf[pdb+chain]=pearsonr(consurfs_list,slopes_list)[0]
                out.write(yeast+' '+pdb+' '+chain+' '+str(protein_corr_cor_consurf[pdb+chain])+' '+str(protein_corr_slope_consurf[pdb+chain])+'\n')
                #print yeast,pdb,chain,protein_corr_cor_consurf[pdb+chain],protein_corr_slope_consurf[pdb+chain]
                overall_dataset_ranks=overall_dataset_ranks+ranks_list
                overall_dataset_pearsons=overall_dataset_pearsons+correlations_list
                overall_dataset_slopes=overall_dataset_slopes+slopes_list
                overall_dataset_ranked_pearsons=overall_dataset_ranked_pearsons+ranked_correlations_list
                overall_dataset_ranked_slopes=overall_dataset_ranked_slopes+ranked_slopes_list
            
            except:
                print 'problem', yeast,pdb,chain
                pass
       
    overall_dataset_correlaion_rank_slope=pearsonr(overall_dataset_ranks,overall_dataset_slopes)[0]
    print 'overall_dataset_correlaion_rank_slope',overall_dataset_correlaion_rank_slope
    overall_dataset_correlaion_rank_pearson=pearsonr(overall_dataset_ranks,overall_dataset_pearsons)[0]
    print 'overall_dataset_correlaion_rank_pearsons',overall_dataset_correlaion_rank_pearson
# overall dataset average pearson vs binned ranks
#    y=binned_statistic(overall_dataset_ranks,overall_dataset_pearsons,statistic='mean',bins=[0.  , 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1 , 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2 , 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3 , 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4 , 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5 , 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6 , 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7 , 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8 , 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9 , 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.  ])
# all residues!!
    binned_residues={}
#    for i in zip(overall_dataset_ranks,overall_dataset_pearsons):
#        print i[0],i[1]
    y=binned_statistic(overall_dataset_ranks,overall_dataset_pearsons,statistic='mean',bins=[0.  , 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50 , 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0  ])
    counter=-1
    for i in y[2]: # y[2] is the assignment of bin number to each residue
        counter=counter+1
        if i not in binned_residues:
            binned_residues[i]=[]
            binned_residues[i].append(overall_dataset_pearsons[counter])
        else:
            binned_residues[i].append(overall_dataset_pearsons[counter])
    df=pd.DataFrame.from_dict(binned_residues,orient='index')
    print df
    df=df.transpose()
    print df
    df.to_csv('figure1a.txt',sep=',')
    return


#run('chosen_BSnoncat_corrected','../ConsurfDB_grades_with_functional_sites/summary_all_BSnoncat_correct','BSnoncat_pearsons_between_pearsonORslopes_and_consurf')
run('chosen_structural_homologs','consurfDB_file','figure1b.txt')

