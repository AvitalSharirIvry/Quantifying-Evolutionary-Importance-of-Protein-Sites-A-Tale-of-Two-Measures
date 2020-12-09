from functions import rank_consurf_file,separate_summary_consurfDB,pearsonr_slope_dist_cons_inPDB,read_propert_over_curated_file
from scipy.stats.stats import pearsonr,spearmanr,linregress
from scipy.stats import binned_statistic
import pandas as pd

def read_funct_res(list_file):
    list_file=open(list_file)
    funct_residues={}
    for line in list_file:
        yeast,pdb,chain=line.split()[0],line.split()[1],line.split()[2]
        if (line.split()[3]).isalpha()==True:
            lig=line.split()[3]
            funct_ress=line.split()[4:]
        else:
            funct_ress=line.split()[3:]
        funct_residues[yeast+pdb+chain]=funct_ress
    return funct_residues

def run(chosen_file,summary_consurf_file,funct_res_file):
    chosen_file=open(chosen_file)
    Cax_of_pdb,seq_of_pdb=read_propert_over_curated_file('Cax_over_curated_chains')
    consurfDB_file_lines=separate_summary_consurfDB(summary_consurf_file)
    overall_dataset_ranks,overall_dataset_pearsons,overall_dataset_slopes,overall_dataset_ranked_pearsons,overall_dataset_ranked_slopes=[],[],[],[],[]
# read functional residues
    functional_ress=read_funct_res(funct_res_file)
# for each protein - yeast with model
    yeast_pdb_chain=[]
    for line in chosen_file:
        yeast,pdb,chain=line.split()[0],line.split()[1].lower(),line.split()[2]
        if pdb+chain in consurfDB_file_lines and pdb+chain in Cax_of_pdb and yeast+pdb+chain in functional_ress and yeast+pdb+chain not in yeast_pdb_chain:
            lines=consurfDB_file_lines[pdb+chain]
            yeast_pdb_chain.append(yeast+pdb+chain)
    # get the ranks and sorted scores of all the residues in the protein
            (zipped_sorted_consurfs,zipped_ranks,ranks,consurfs,most_conserved_res,res_of_interest_most_conserved)=rank_consurf_file(lines,[])
    #Pearson correlation (conservation,dist) for each residue
            correlations,pval_correlations,slopes,ranked_pearsons,ranked_slopes=pearsonr_slope_dist_cons_inPDB(pdb,chain,lines)       

            consurfs_list=[]    # lists of residues: consurfs,earsons,slopes,ranks
            correlations_list=[]
            slopes_list=[]
            ranks_list=[]
            ranked_correlations_list=[]
            ranked_slopes_list=[]
            if yeast+pdb+chain in functional_ress:
                print yeast,pdb,chain
                for res in functional_ress[yeast+pdb+chain]:    # populate the lists
                    if res in correlations:         # make sure it has a correlations (appears in consurf)
                        consurfs_list.append(consurfs[res])
                        correlations_list.append(correlations[res])
                        slopes_list.append(slopes[res])
                        ranks_list.append(ranks[res])
                        ranked_correlations_list.append(ranked_pearsons[res])
                        ranked_slopes_list.append(ranked_slopes[res])
                try:
                # concatenate the lists to an overall list
                    overall_dataset_ranks=overall_dataset_ranks+ranks_list
                    overall_dataset_pearsons=overall_dataset_pearsons+correlations_list
                    overall_dataset_slopes=overall_dataset_slopes+slopes_list
                    overall_dataset_ranked_pearsons=overall_dataset_ranked_pearsons+ranked_correlations_list
                    overall_dataset_ranked_slopes=overall_dataset_ranked_slopes+ranked_slopes_list

                except:
                    print 'problem', yeast,pdb,chain
                    pass
# calculate the overall plot of residue pearson vs rank     &   resiodue slope vs rank over ALL proteins
#    overall_dataset_correlaion_rank_slope=pearsonr(overall_dataset_ranks,overall_dataset_slopes)[0]
#    print 'overall_dataset_correlaion_rank_slope',overall_dataset_correlaion_rank_slope
#    overall_dataset_correlaion_rank_pearson=pearsonr(overall_dataset_ranks,(overall_dataset_pearsons)[0])
#    print 'overall_dataset_correlaion_rank_pearsons',overall_dataset_correlaion_rank_pearson
#    print 'all residues allos'
 #   for i in zip(overall_dataset_ranks,overall_dataset_pearsons):
 #       print i[0],i[1]
# overall dataset average pearson vs binned ranks 
    binned_residues={}

    y=binned_statistic(overall_dataset_ranks,overall_dataset_pearsons,statistic='mean',bins=[0.  , 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50 , 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0  ])
    counter=-1
    for i in y[2]:
        counter=counter+1
        if i not in binned_residues:
            binned_residues[i]=[]
            binned_residues[i].append(overall_dataset_pearsons[counter])
        else:
            binned_residues[i].append(overall_dataset_pearsons[counter])
    df=pd.DataFrame.from_dict(binned_residues,orient='index')
    df=df.transpose()
    df.to_csv('figure2_BSnoncat.txt',sep=',')

    return

#run('chosen_BSnoncat_corrected','../ConsurfDB_grades_with_functional_sites/summary_all_BSnoncat_correct','BSnoncat_pearsons_between_pearsonORslopes_and_consurf')
run('chosen_structural_homologs','consurfDB_file','../list_res_BSnoncat_no_allos')

