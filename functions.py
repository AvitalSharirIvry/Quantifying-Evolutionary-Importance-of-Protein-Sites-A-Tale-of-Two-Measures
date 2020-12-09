
import random
import numpy
import numpy 
import os 
import gzip 
from pattern_align_local import align 
from scipy.stats.stats import pearsonr,spearmanr,linregress
from scipy.optimize import curve_fit

def intersections_of_lists(list1,list2):
    intersection=[value for value in list1 if value in list2]
    return intersection

def read_HGMD():
	three_letters_to_1letters_aa_text=open('/home/avitaliv/projects/ctb-yxia/avitaliv//xialab/Data_files/3letters_to-1letters_aa')
	HGMD_text=open('/home/avitaliv/projects/ctb-yxia/avitaliv//xialab/Data_files/HGMD_database_2011_with_header.txt')
# read HGMD data
	HGMD_polymorphism={}
	HGMD_disease={}
	HGMD_genes=[]
	HGMD_diseaseAssociated={}
        refseq_ofGene={}
	for line in HGMD_text:
		if line.split()[0]=='Type':
			break
	for line in HGMD_text:
		gene=line.split()[8]
                refseq_accession=line.split()[7].split(':')[0]
                refseq_ofGene[gene]=refseq_accession
		if gene not in HGMD_genes:
			HGMD_genes.append(gene)
		if gene=='-':
			break
		variant_class=line.split()[1]
		if line.split()[7]!='null':
			AA_change=line.split()[7].split('.')[2]
		else:	
			AA_change=line.split()[7]
		if variant_class=='FP':
			if gene in HGMD_polymorphism:
				HGMD_polymorphism[gene].append(AA_change)
			else:
				HGMD_polymorphism[gene]=[]
				HGMD_polymorphism[gene].append(AA_change)
		elif variant_class=='DM':
			if gene in HGMD_disease:
				HGMD_disease[gene].append(AA_change)
			else:
				HGMD_disease[gene]=[]
				HGMD_disease[gene].append(AA_change)
		elif variant_class=='DFP' or variant_class=='DP':
			if gene in HGMD_diseaseAssociated:
				HGMD_diseaseAssociated[gene].append(AA_change)
			else:
				HGMD_diseaseAssociated[gene]=[]
				HGMD_diseaseAssociated[gene].append(AA_change)
# read 3 letter code to 1 letter code
	three_to_1={}
	for line in three_letters_to_1letters_aa_text:
		three_to_1[line.split()[0]]=line.split()[1]

	return (HGMD_polymorphism,HGMD_disease,refseq_ofGene)



def identify_cat_res_general():
    cat_res_file=open('/home/avitaliv/projects/ctb-yxia/avitaliv//xialab/Data_files/catalytic_res_atlas_wout_homs_2018.txt')
    cat_res_ofPDB,cat_res_nm_ofPDB={},{}
    lines_cat_res_file=cat_res_file.readlines()
    for line in lines_cat_res_file[1:]:
        pdb=line.split(',')[0]
        chain=line.split(',')[3]
        if pdb+chain in cat_res_ofPDB:
            cat_res_ofPDB[pdb+chain].append(line.split(',')[4])
            cat_res_nm_ofPDB[pdb+chain].append(line.split(',')[2].upper()+line.split(',')[4])
        else:
            cat_res_ofPDB[pdb+chain]=[]
            cat_res_nm_ofPDB[pdb+chain]=[]
            cat_res_ofPDB[pdb+chain].append(line.split(',')[4])
            cat_res_nm_ofPDB[pdb+chain].append(line.split(',')[2].upper()+line.split(',')[4])
    return cat_res_ofPDB,cat_res_nm_ofPDB
def identify_cat_res_general_w_added_homs():
    cat_res_file=open('/home/avitaliv/projects/ctb-yxia/avitaliv//xialab/Data_files/catalytic_res_atlas_w_homs_from_assemblies')
    cat_res_ofPDB={}
    lines_cat_res_file=cat_res_file.readlines()
    for line in lines_cat_res_file:
        pdb,chain,ress=line.split()[0],line.split()[1],line.split()[2:]
        cat_res_ofPDB[pdb+chain]=ress
    return cat_res_ofPDB

def identify_cat_res_general_oldV():
    cat_res_file=open('/home/avitaliv/projects/ctb-yxia/avitaliv//xialab/Data_files/catalytic_res_atlas.txt')
    cat_res_ofPDB,cat_res_nm_ofPDB={},{}
    lines_cat_res_file=cat_res_file.readlines()
    for line in lines_cat_res_file[1:]:
        pdb=line.split(',')[0]
        chain=line.split(',')[3]
        if pdb+chain in cat_res_ofPDB:
            cat_res_ofPDB[pdb+chain].append(line.split(',')[4])
            cat_res_nm_ofPDB[pdb+chain].append(line.split(',')[2].upper()+line.split(',')[4])
        else:
            cat_res_ofPDB[pdb+chain]=[]
            cat_res_nm_ofPDB[pdb+chain]=[]
            cat_res_ofPDB[pdb+chain].append(line.split(',')[4])
            cat_res_nm_ofPDB[pdb+chain].append(line.split(',')[2].upper()+line.split(',')[4])
    return cat_res_ofPDB,cat_res_nm_ofPDB

#pdb_chains_of_genes is a list of gene_pdbchain items
def identify_cat_res(gene_pdb_chain_list):
    cat_res_file=open('/home/avitaliv/projects/ctb-yxia/avitaliv//xialab/Data_files/catalytic_res_atlas.txt')
    cat_res_of_gene_pdb={}
    cat_res_ofPDB={}
    lines_cat_res_file=cat_res_file.readlines()
    for line in lines_cat_res_file[1:]:
        pdb=line.split(',')[0]
        chain=line.split(',')[3]
        if pdb+chain in cat_res_ofPDB:
            cat_res_ofPDB[pdb+chain].append(line.split(',')[4])
        else:
            cat_res_ofPDB[pdb+chain]=[]
            cat_res_ofPDB[pdb+chain].append(line.split(',')[4])

    for case in gene_pdb_chain_list:
        gene=case.split('_')[0]
        pdb_chain=case.split('_')[1]
        if pdb_chain in cat_res_ofPDB:
            cat_res_of_gene_pdb[case]=cat_res_ofPDB[pdb_chain]
    
    return cat_res_of_gene_pdb

def reduced_aligned_res_to_cat(aligned_res,out):
    aligned_file=open(aligned_res)
    out=open(out,'w')
    cat_res_file=open('/home/avitaliv/projects/ctb-yxia/avitaliv//xialab/Data_files/catalytic_res_atlas.txt')
    cat_res_ofPDB={}
    lines_cat_res_file=cat_res_file.readlines()
    for line in lines_cat_res_file[1:]:
        pdb=line.split(',')[0]
        chain=line.split(',')[3]
        if pdb+chain in cat_res_ofPDB:
            cat_res_ofPDB[pdb+chain].append(line.split(',')[4])
        else:
            cat_res_ofPDB[pdb+chain]=[]
            cat_res_ofPDB[pdb+chain].append(line.split(',')[4])

    for line in aligned_file:
        pdb_chain=line.split()[0].split('_')[1]
        if pdb_chain in cat_res_ofPDB:
            out.write(line)
    out.close()
    return

def reduce_aligned_res_to_bestpdb(aligned_res,out):
    aligned_file=open(aligned_res)
    out=open(out,'w')
    lines=aligned_file.readlines()
    pdb_seqs_of_gene={}
    for line in lines:
        gene_pdb=line.split()[0]
        gene=line.split()[0].split('_')[0]
        if gene not in pdb_seqs_of_gene:
            pdb_seqs_of_gene[gene]=[]
            pdb_seqs_of_gene[gene].append(line)
        else:
            pdb_seqs_of_gene[gene].append(line)

    for gene in pdb_seqs_of_gene:
      #  print gene
        best_line=pdb_seqs_of_gene[gene][0]
        gaps_num=best_line.count('-- ')
       # print best_line,gaps_num
        for line in pdb_seqs_of_gene[gene]:
           # print line.count('-- ')
            if line.count('-- ')<gaps_num:
                gaps_num=line.count('-- ')
                best_line=line
        out.write(best_line)
    out.close()
    return

# provides mut_gene_pdb[gene_pdbchain]=list_pdb_res_nums_which_are_muts
def identify_mut_in_pdb(type_mut,aligned_res):

    from functions import read_HGMD
    (hgmd_polymorphism,hgmd_disease,refseq_ofgene)=read_HGMD()

    mut_gene_pdb={}
    aligned_pdb_res_to_refseq_file=open(aligned_res)
    lines=aligned_pdb_res_to_refseq_file.readlines()
    for line in lines:
        gene_pdb=line.split()[0]
        gene=line.split('_')[0]
        refseq_seq=line.split()[1]
        hgmd_mut=hgmd_disease
        if gene in hgmd_mut:
            muts=hgmd_mut[gene]
            muts=[mut for mut in muts if mut!='null']
            for mut in muts :
                index_mut_on_refseq=int(mut[1:-1])
                mut_refseq_res=mut[0]
                if len(refseq_seq)>index_mut_on_refseq and refseq_seq[index_mut_on_refseq-1]==mut_refseq_res:
                    index_mut_on_pdb=str(line.split()[2:][index_mut_on_refseq-1])
                    if gene_pdb in mut_gene_pdb:
                        mut_gene_pdb[gene_pdb].append((index_mut_on_refseq,index_mut_on_pdb))
                    else:
                        mut_gene_pdb[gene_pdb]=[]
                        mut_gene_pdb[gene_pdb].append((index_mut_on_refseq,index_mut_on_pdb))
# mut_gene_pdb[gene_pdb]=[(res in refseq,res in pdb),(),...]
    return mut_gene_pdb
                
def wget_pdb(pdb):
    pdb_code=pdb+'.pdb'
    basic_path='http://www.rcsb.org/pdb/files/'
    path=basic_path+pdb_code
    os.system('wget --no-check-certificate -P pdb_files/ %s' %path )
    return

def calc_dist_from_res_list(pdb_chain,res_list):
    min_dist_from_res_list={}
    pdb_code=pdb_chain[:-1]+'.pdb1.gz'
    basic_path='/home/avitaliv/projects/ctb-yxia/avitaliv//singen_largerYdata2/output/assemblies_curated/'
    path=basic_path+pdb_code
    with gzip.open(path,'r') as f:
#		os.system('wget --no-check-certificate -P pdb_files/ %s' %path )	
        lines_f_pdb=f.readlines()
        res_list_x=[]
        res_list_y=[]
        res_list_z=[]
        for line in lines_f_pdb:
            if line.split()[0]=='ATOM' and 'CA' in line and line[21]==pdb_chain[-1] and line[22:26].strip() in res_list:
	        res_list_x.append(float(line[30:37]))
	        res_list_y.append(float(line[38:45]))
	        res_list_z.append(float(line[46:53]))
    res_num_used=[] # checker to make sure we will not take the same residue several times with different x,y,z (we take the 1st confor as default)
    for line in lines_f_pdb:
        if line.split()[0]=='ATOM' and line[13:15]=='CA' and line[21]==pdb_chain[-1]: 
	    x=float(line[30:37])
	    y=float(line[38:45])
	    z=float(line[46:53])
	    res_num=line[22:26].strip()
	    res_nm=line[17:20]
            dist_list=[]
	    counter=-1
# make sure that no duplicate residues are taken in cases where there are a few conformations in the pdb (1st is default) (CA hardly affected anyway)
	    if res_num in res_num_used:
                res_num_used.append(res_num)    # collecting all used residues (prevent repeat due to diff conf of residues etc)
            else:
		for X in res_list_x:
		    counter=counter+1
		    dist=numpy.sqrt((x-X)**2+(res_list_y[counter]-y)**2+(res_list_z[counter]-z)**2) 
		    dist_list.append(dist)
		if dist_list!=[]:
		    min_dist_from_cat=min(dist_list)
                    min_dist_from_res_list[res_num]=min_dist_from_cat

        else:
            continue


#        os.system('rm pdb_files/*')
    return  min_dist_from_res_list

def calc_dist_from_res_list_based_on_Ca_over_entire_pdb(pdb_chain):
    f_aligned_Cax=open('/home/avitaliv/projects/ctb-yxia/avitaliv//singen_largerYdata2/src_fix_chosen_hom/t_Cax_over_curated_chains')
    f_aligned_Cay=open('/home/avitaliv/projects/ctb-yxia/avitaliv//singen_largerYdata2/src_fix_chosen_hom/t_Cay_over_curated_chains')
    f_aligned_Caz=open('/home/avitaliv/projects/ctb-yxia/avitaliv//singen_largerYdata2/src_fix_chosen_hom/t_Caz_over_curated_chains')
    f_aligned_res_num=open('/home/avitaliv/projects/ctb-yxia/avitaliv//singen_largerYdata2/src_fix_chosen_hom/res_num_over_curated_chains')
    Cax={}
    Cay={}
    Caz={}
    res_num={}
    for line in f_aligned_Cax:
        if len(line.split())>3:
            Cax[line.split()[0]+line.split()[1]]=line.split()[3].split('|')[:-1]
    for line in f_aligned_Cay:
        if len(line.split())>3:
            Cay[line.split()[0]+line.split()[1]]=line.split()[3].split('|')[:-1]
    for line in f_aligned_Caz:
        if len(line.split())>3:
            Caz[line.split()[0]+line.split()[1]]=line.split()[3].split('|')[:-1]
    for line in f_aligned_res_num:
        if len(line.split())>3:
            res_num[line.split()[0]+line.split()[1]]=line.split()[3].split('|')[:-1]
   # print Cax['4kqxA']
    zipped_dicts=zip(Cax[pdb_chain],Cay[pdb_chain],Caz[pdb_chain],res_num[pdb_chain])
   # print zipped_dicts
    dict_of_dists_from_res_of_interest={}
    for cat_res in zipped_dicts: # this is the ref residue
        dists_from_cat_res={}
        for residue in zipped_dicts:
            d=numpy.sqrt((float(residue[0])-float(cat_res[0]))**2+(float(residue[1])-float(cat_res[1]))**2+(float(residue[2])-float(cat_res[2]))**2)
            dists_from_cat_res[residue[3]]=d
        dict_of_dists_from_res_of_interest[cat_res[3]]=dists_from_cat_res
    return dict_of_dists_from_res_of_interest # dists_from_res_of_interest[res]=min_dist_from_residues_of_interest


def find_closest_res_of_res(pdb_inp,chain_inp):
    pdb_code=pdb_inp+'.pdb'
    f_pdb=open('pdb_files/'+pdb_code)
    coordinates=[]
    closest_res={}
    for line in f_pdb:
        if line.split()[0]=='ATOM' and 'CA' in line and line[21] in chain_inp:
            x=float(line[30:37])
            y=float(line[38:45])
            z=float(line[46:53])
            res_num=line[22:26].strip()
            coordinates.append(res_num,x,y,z)
            closest_res[res_num]=[]
    for i in range(0,len(coordinates)+1):
        for res in coordinates[i+1:]:
            d=numpy.sqrt((coordinates[i][1]-res[1])**2+(coordinates[i][2]-res[2])**2+(coordinates[i][3]-res[3])**2)
            if d<7.0:
                closest_res[coordinates[i][0]].append(res[0])
                closest_res[res[0]].append(coordinates[i][0])
    print closest_res

def extract_list_pdb_res_from_aligned(aligned_res):
    aligned_file=open(aligned_res)
    pdb_residues={}
    for line in aligned_file:
        residues=[r for r in line.split()[2:] if '--' not in r]
        pdb_gene=line.split()[0]
        pdb_residues[pdb_gene]=residues
#pdb_residues[pdb_gene]= list of residues
    return pdb_residues


def read_dssp(pdb_chain):
    Max=read_maxsasa()
    ACC_ofPdb={}
    RSA_ofPdb={}
    pdb=pdb_chain[:-1]
    chain=pdb_chain[-1]
    file=open('/home/avitaliv/projects/ctb-yxia/avitaliv//xialab/Data_files/dssp_files/'+pdb+'.dssp')
    lines=file.readlines()
    counter_line=-1
    for line in lines:
         counter_line=counter_line+1
         if '#  RESIDUE AA STRUCTURE' in line:
             start_line=counter_line
    for line in lines[start_line:]:
        if line[11]==chain:
            AA,resNum,ACC=line[13],line[6:10].strip(),line[35:38].strip()
            if AA in Max:
                ACC_ofPdb[resNum]=ACC
                RSA_ofPdb[resNum]=str(float(ACC)/float(Max[AA]))
    return ACC_ofPdb,RSA_ofPdb # ACC_ofPdb[resNum]=ACC_value, RSA_ofPdb[resNum]=RSA_value

def read_enz_pdbs_Brenda():
    list_enz_pdbs=[]
    Brenda_file=open('/home/avitaliv/projects/ctb-yxia/avitaliv//xialab/Data_files/BRENDA_EC_PDB')
    for line in Brenda_file:
        list_enz_pdbs.append(line.split()[1]+line.split()[2])
    return list_enz_pdbs

def resLetter_of_resNum_from_aligned(aligned_res):
    f=open(aligned_res)
    Let_ofResNum={}
    for line in f:
        let_of_resNum={}
        gene_pdb=line.split()[0]
        seq=line.split()[1]
        res_nums=line.split()[2:]
        i=-1
        while i<len(seq)-1:
            i=i+1
            let_of_resNum[res_nums[i]]=seq[i]
        Let_ofResNum[gene_pdb]=let_of_resNum
    return Let_ofResNum

def read_maxsasa():
    f=open('/home/avitaliv/projects/ctb-yxia/avitaliv//singen_largerYdata2/data/maxsasa.dat')
    Max={}
    for line in f:
        Max[line.split()[0]]=float(line.split()[1])
    return Max

def bootstrap_array(array):
    counter=0
    shuffled_arrays=[]
    l=len(array)
    while counter<10:
        counter=counter+1
        shuffled_array=[random.choice(array) for _ in range(len(array))] # random.choice is sampling with replacement
        shuffled_arrays.append(shuffled_array)
    return shuffled_arrays

def read_column_from_file(input_file,column_num):
    array=[float(x.split()[column_num]) for x in open(input_file).readlines()]
    return array

def create_RSA_histogram(array):
    RSA_bins=[-1.0,0.00000001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.5]
    hist,edges=numpy.histogram(array,bins=RSA_bins)
    return hist, edges

def read_aligned_res_num(aligned_res_num_file):
    aligned_res_num_file=open(aligned_res_num_file)
    res_nums_of_y={}
    seq_of_y={}
    for line in aligned_res_num_file:
        yeast=line.split()[0]
        seq=line.split()[1]
        res_nums=line.split()[2:]
        res_nums_of_y[yeast]=res_nums
        seq_of_y[yeast]=seq
    return res_nums_of_y,seq_of_y

def three_to_one_amino():
    file=open('/home/avitaliv/projects/ctb-yxia/avitaliv//xialab/Data_files/3letters_to-1letters_aa')
    three_to_1={}
    for line in file:
        three_to_1[line.split()[0].upper()]=line.split()[1]
    return three_to_1

def read_pdb_gzipped(path,zipped_file_nm,pdb_chain):
    three_to_1=three_to_one_amino()
    seq=''
    res_nums=''
    with gzip.open(path+zipped_file_nm,'r') as f:
        for line in f:
            if line.split()[0]=='ATOM' and line.split()[2]=='CA' and line[21]==pdb_chain[-1] and line[16]!='B': # only A conformations in the sequence
                seq=seq+three_to_1[line[17:20]]
                res_nums=res_nums+line[22:26].strip()+'|'

    return seq,res_nums
def reduce_aligned1_according_to_aligned2(aligned1,aligned2,out_aligned):
    f1=open(aligned1)
    f2=open(aligned2)
    out=open(out_aligned,'w')
    list1=[]
    list2=[]
    for line in f2:
        list2.append(line.split()[0])
    for line in f1:
        if line.split()[0] in list2:
            out.write(line)

    out.close()
    return
def reduce_aligned1_according_to_not_aligned2(aligned1,aligned2,out_aligned):
    f1=open(aligned1)
    f2=open(aligned2)
    out=open(out_aligned,'w')
    list1=[]
    list2=[]
    for line in f2:
        list2.append(line.split()[0])
    for line in f1:
        if line.split()[0] not in list2:
            out.write(line)

    out.close()
    return

def rank_consurf_file(consurf_file_lines,list_res_of_interest):
    ranks,consurfs,list_ranks,list_res,list_consurf={},{},[],[],[]
    for line in consurf_file_lines[15:-4]:
        if line.split()[2]!='-':
            consurf,res=float(line.split()[3]),line.split()[2].split(':')[0][3:]
            list_consurf.append(consurf)
            list_res.append(res)
            consurfs[res]=consurf
    l=len(list_consurf)
    zipped_sorted_consurfs=zip(list_consurf,list_res)
    zipped_sorted_consurfs.sort()
    counter=-1
    sorted_list_res=[]
    for res in zipped_sorted_consurfs:
        counter=counter+1
        rank=float(counter)/float(l)
        list_ranks.append(rank)
        sorted_list_res.append(res[1])
        ranks[res[1]]=rank
    zipped_ranks=zip(list_ranks,sorted_list_res)
    most_conserved_res=zipped_sorted_consurfs[0][1]
    res_of_interest_most_conserved=''
    for i in zipped_sorted_consurfs:
        if i[1] in list_res_of_interest:
            res_of_interest_most_conserved=i[1]
            break
    return (zipped_sorted_consurfs,zipped_ranks,ranks,consurfs,most_conserved_res,res_of_interest_most_conserved)  # list_of_(cons,res)_sorted, most_conserved_res_of_all_pdb_file, most conseved_res out of the specific list of interest


def separate_summary_consurfDB(summary_consurfDB_file):
    summary_consurfDB_file=open(summary_consurfDB_file)
    consurfDB_file_lines={}
    for line in summary_consurfDB_file:
        if '_consurf_summary.txt' in line:
            pdb_chain=line[0:4].lower()+line[4]
            consurfDB_file_lines[pdb_chain]=[]
        else:
            consurfDB_file_lines[pdb_chain].append(line)
    return consurfDB_file_lines

def read_hits_file(hits_file):
    yeast_models={}
    f=open(hits_file)
    for line in f:
        if line.split()[0] not in yeast_models:
            yeast_models[line.split()[0]]=[]
            yeast_models[line.split()[0]].append(line.split()[1][:4]+line.split()[1][5])
        else:
            yeast_models[line.split()[0]].append(line.split()[1][:4]+line.split()[1][5])
    return yeast_models
def read_ASD_allosetric_sites_file():
    f=open('/home/avitaliv/projects/ctb-yxia/avitaliv//xialab/Data_files/ASD_Release_062015_XF/allosteric_list')
    allosteric_sites={}
    for line in f:
        temp_sites=[]
        for data in line.split()[1].split(';'):
            if 'site' in data:
                temp_sites.append(data.split(',')[1])
            if 'range' in data:
                upper=int(data.split(',')[2])
                lower=int(data.split(',')[1])
                for i in range(lower,upper+1):
                    temp_sites.append(str(i))
        allosteric_sites[line.split()[0]]=temp_sites
    return allosteric_sites # allosteric_sites[ASD_case]=list_of_allos_res

def read_ASD_pdb_per_ASD_file():
    f=open('/home/avitaliv/projects/ctb-yxia/avitaliv//xialab/Data_files/ASD_Release_062015_XF/pdb_id_list')
    allosteric_sites=read_ASD_allosetric_sites_file() # needs this function to be loaded!!
    allosteric_pdb={}
    for line in f:
        if line.split()[0] in allosteric_sites:
            if line.split()[0] not in allosteric_pdb:
                allosteric_pdb[line.split()[0]]=[]
            if ';' in line.split()[1]:
                temp_pdb_list=line.split()[1].split(';')
                for p in temp_pdb_list:
                    allosteric_pdb[line.split()[0]].append(p)
            else:
                allosteric_pdb[line.split()[0]].append(line.split()[1])
    return allosteric_pdb # allosteric_pdb[ASD_case]=list_of_pdbs
def read_propert_over_curated_file(property_over_curated_file): # important!! maybe parameters here are different for different files (RSA,Cax etc)
    property_over_curated_file=open(property_over_curated_file)
    property_of_pdb={}
    seq_of_pdb={}
    for line in property_over_curated_file:
        pdb_chain=line.split()[0]+line.split()[1]
        if len(line.split())>2:
            seq=line.split()[2] # !!! changable !!!
            propert=line.split()[3].replace('|',' ') #!!! changable!!), also the -1
      #      propert=line.split()[4].replace('|',' ')[:-1] #!!! changable!!), also the -1
            property_of_pdb[pdb_chain]=propert
            seq_of_pdb[pdb_chain]=seq

    return property_of_pdb,seq_of_pdb
def read_propert_RSA_over_curated_file(property_over_curated_file): # important!! only for RSA 
    property_over_curated_file=open(property_over_curated_file)
    RSA_of_pdb={}
    seq_of_pdb={}
    for line in property_over_curated_file:
        pdb_chain=line.split()[0]+line.split()[1]
        if len(line.split())>3:
            seq=line.split()[2] 
            propert=line.split()[3].replace('|',' ')[:-1]
            RSA_of_pdb[pdb_chain]=propert
            seq_of_pdb[pdb_chain]=seq
    return RSA_of_pdb,seq_of_pdb

def read_Scer():
    yeast_seq={}
    Scer_file=open('/home/avitaliv/projects/ctb-yxia/avitaliv//singen_largerYdata2/data/Scer.aa')
    for line in Scer_file:
        if '>' in line:
            yeast=line[1:-2]
            yeast_seq[yeast]=next(Scer_file)[:-2]
    return yeast_seq # yeast_seq[yeast]=seq

def read_chosen_homologs(chosen_homs_file):
    chosen_homs_file=open(chosen_homs_file)
    chosen_hom={}
    for line in chosen_homs_file:
        chosen_hom[line.split()[0]]=line.split()[1]+line.split()[2][0]
    return chosen_hom

def read_chosen_homologs_into_list(chosen_homs_file):
    chosen_homs_file=open(chosen_homs_file)
    chosen_hom=[]
    for line in chosen_homs_file:
        chosen_hom.append(line.split()[1]+line.split()[2][0])
    return chosen_hom

def align_property(chosen_homs_file,property_over_curated_file,out):
    out=open(out,'w')
#    property_of_pdb,seq_of_pdb=read_propert_over_curated_file(property_over_curated_file)
    property_of_pdb,seq_of_pdb=read_propert_RSA_over_curated_file(property_over_curated_file) 
    yeast_seq=read_Scer()
    chosen_hom=read_chosen_homologs(chosen_homs_file)
    for yeast in chosen_hom:
        pdb_chain=chosen_hom[yeast]
        Y_seq=yeast_seq[yeast]
        if pdb_chain in seq_of_pdb:
            PDB_seq=seq_of_pdb[pdb_chain]
            Property=property_of_pdb[pdb_chain]
#            print pdb_chain,Y_seq,PDB_seq,Property
            try:
	    #	print yeast,pdb_chain,Property
            	alignment_out=align(Y_seq,PDB_seq,Property)
            	out.write(yeast + ' ' + Y_seq+ ' '+' '.join(alignment_out)+ '\n')
            except:
                print 'Error'
                pass
    out.close()
    return

def reduce_dict_keys_according_to_list(old_dict,reducer_list):
    new_dict={k: old_dict[k] for k in reducer_list if k in old_dict}
    return new_dict

def pearsonr_slope_dist_cons_inPDB(pdb,chain,consurfDB_file_lines):
    consurf_scores={}
    res_nums=[]
    for line in consurfDB_file_lines:
        if line[0]!='-' and len(line.split())>0 and 'Below' not in line and 'estimated' not in line and '==' not in line and 'Amino' not in line and '3LATOM' not in line and 'normalized'  not in line and line.split()[2]!='-':
	    res_num=line.split()[2][3:].split(':')[0]
	    consurf_scores[res_num]=float(line.split()[3])
            res_nums.append(res_num)
    correlations={}
    slopes={}
    pval_correlations={}
    combined_cons_cor={}
    ranked_pearsons,ranked_slopes={},{}

    dict_of_min_dist_from_res_list=calc_dist_from_res_list_based_on_Ca_over_entire_pdb(pdb+chain)
    for res in res_nums:	# go over residues that are relevant . BUT - maybe some have consurf score but no dist..(will be taken care below)
        dists=[]
        consurf_scores_wout_problems={}
        consurf_scores_wout_problems_list=[]
        if res in dict_of_min_dist_from_res_list:
	    min_dist_from_res_list=dict_of_min_dist_from_res_list[res]
            for res2 in res_nums:
                if res2 in min_dist_from_res_list :
	            dists.append(float(min_dist_from_res_list[res2]))
	            consurf_scores_wout_problems[res2]=consurf_scores[res2] # only scores for res that have a pdb coordinates
                    consurf_scores_wout_problems_list.append(consurf_scores[res2])
            cor=pearsonr(consurf_scores_wout_problems_list,dists)
            slope=linregress(consurf_scores_wout_problems_list,dists)
        if dists!=[]:
	    correlations[res]=cor[0]
	    pval_correlations[res]=cor[1]
            slopes[res]=slope[0]
    tuple_slopes,tuple_correlations=[(slopes[res],res) for res in slopes],[(correlations[res],res) for res in correlations]
    tuple_slopes.sort()
    tuple_correlations.sort()
  #  sorted_zipped_slopes,sorted_zipped_correlations=tuple_slopes.sort(),tuple_correlations.sort()
    l=float(len(slopes))
    for residue in tuple_correlations:
        res,ranked_pearson=residue[1],float(tuple_correlations.index(residue))/l
        ranked_pearsons[res]=ranked_pearson
    for residue in tuple_slopes:
        res,ranked_slope=residue[1],float(tuple_slopes.index(residue))/l
        ranked_slopes[res]=ranked_slope
    return correlations,pval_correlations,slopes,ranked_pearsons,ranked_slopes #correaltions[res]=pearson_cor   pval_correlations[res]=pval of pearson 

def pearsonr30_slope_dist_cons_inPDB(pdb,chain,consurfDB_file_lines):
    consurf_scores={}
    res_nums=[]
    for line in consurfDB_file_lines:
        if line[0]!='-' and len(line.split())>0 and 'Below' not in line and 'estimated' not in line and '==' not in line and 'Amino' not in line and '3LATOM' not in line and 'normalized'  not in line and line.split()[2]!='-':
	    res_num=line.split()[2][3:].split(':')[0]
	    consurf_scores[res_num]=float(line.split()[3])
            res_nums.append(res_num)
    correlations={}
    slopes={}
    pval_correlations={}
    combined_cons_cor={}
    ranked_pearsons,ranked_slopes={},{}

    dict_of_min_dist_from_res_list=calc_dist_from_res_list_based_on_Ca_over_entire_pdb(pdb+chain)
    for res in res_nums:	# go over residues that are relevant . BUT - maybe some have consurf score but no dist..(will be taken care below)
        dists=[]
        consurf_scores_wout_problems={}
        consurf_scores_wout_problems_list=[]
        if res in dict_of_min_dist_from_res_list:
	    min_dist_from_res_list=dict_of_min_dist_from_res_list[res]
            for res2 in res_nums:
                if res2 in min_dist_from_res_list and min_dist_from_res_list[res2]<30.0:
	            dists.append(float(min_dist_from_res_list[res2]))
	            consurf_scores_wout_problems[res2]=consurf_scores[res2] # only scores for res that have a pdb coordinates
                    consurf_scores_wout_problems_list.append(consurf_scores[res2])
            cor=pearsonr(consurf_scores_wout_problems_list,dists)
            slope=linregress(consurf_scores_wout_problems_list,dists)
        if dists!=[]:
	    correlations[res]=cor[0]
	    pval_correlations[res]=cor[1]
            slopes[res]=slope[0]
    tuple_slopes,tuple_correlations=[(slopes[res],res) for res in slopes],[(correlations[res],res) for res in correlations]
    tuple_slopes.sort()
    tuple_correlations.sort()
  #  sorted_zipped_slopes,sorted_zipped_correlations=tuple_slopes.sort(),tuple_correlations.sort()
    l=float(len(slopes))
    for residue in tuple_correlations:
        res,ranked_pearson=residue[1],float(tuple_correlations.index(residue))/l
        ranked_pearsons[res]=ranked_pearson
    for residue in tuple_slopes:
        res,ranked_slope=residue[1],float(tuple_slopes.index(residue))/l
        ranked_slopes[res]=ranked_slope
    return correlations,pval_correlations,slopes,ranked_pearsons,ranked_slopes #correaltions[res]=pearson_cor   pval_correlations[res]=pval of pearson 
def spearmanr_slope_dist_cons_inPDB(pdb,chain,consurfDB_file_lines):
    consurf_scores={}
    res_nums=[]
    for line in consurfDB_file_lines:
        if line[0]!='-' and len(line.split())>0 and 'Below' not in line and 'estimated' not in line and '==' not in line and 'Amino' not in line and '3LATOM' not in line and 'normalized'  not in line and line.split()[2]!='-':
	    res_num=line.split()[2][3:].split(':')[0]
	    consurf_scores[res_num]=float(line.split()[3])
            res_nums.append(res_num)
    correlations={}
    slopes={}
    pval_correlations={}
    combined_cons_cor={}
    ranked_spearmans,ranked_slopes={},{}

    dict_of_min_dist_from_res_list=calc_dist_from_res_list_based_on_Ca_over_entire_pdb(pdb+chain)
    for res in res_nums:	# go over residues that are relevant . BUT - maybe some have consurf score but no dist..(will be taken care below)
        dists=[]
        consurf_scores_wout_problems={}
        consurf_scores_wout_problems_list=[]
        if res in dict_of_min_dist_from_res_list:
	    min_dist_from_res_list=dict_of_min_dist_from_res_list[res]
            for res2 in res_nums:
                if res2 in min_dist_from_res_list :
	            dists.append(float(min_dist_from_res_list[res2]))
	            consurf_scores_wout_problems[res2]=consurf_scores[res2] # only scores for res that have a pdb coordinates
                    consurf_scores_wout_problems_list.append(consurf_scores[res2])
            cor=spearmanr(consurf_scores_wout_problems_list,dists)
            slope=linregress(consurf_scores_wout_problems_list,dists)
        if dists!=[]:
	    correlations[res]=cor[0]
	    pval_correlations[res]=cor[1]
            slopes[res]=slope[0]
    tuple_slopes,tuple_correlations=[(slopes[res],res) for res in slopes],[(correlations[res],res) for res in correlations]
    tuple_slopes.sort()
    tuple_correlations.sort()
  #  sorted_zipped_slopes,sorted_zipped_correlations=tuple_slopes.sort(),tuple_correlations.sort()
    l=float(len(slopes))
    for residue in tuple_correlations:
        res,ranked_spearman=residue[1],float(tuple_correlations.index(residue))/l
        ranked_spearmans[res]=ranked_spearman
    for residue in tuple_slopes:
        res,ranked_slope=residue[1],float(tuple_slopes.index(residue))/l
        ranked_slopes[res]=ranked_slope
    return correlations,pval_correlations,slopes,ranked_spearmans,ranked_slopes #correaltions[res]=spearman_cor   pval_correlations[res]=pval of pearson 
def print_corrConsDist_ofres(pdb,chain,consurfDB_file_lines,list_of_res):
    consurf_scores={}
    res_nums=[]
    for line in consurfDB_file_lines[15:-4]:
	res_num=line.split()[2][3:].split(':')[0]
	consurf_scores[res_num]=float(line.split()[3])
        res_nums.append(res_num)
    correlations,slopes,pval_correlations,combined_cons_cor={},{},{},{}
    dict_of_min_dist_from_res_list=calc_dist_from_res_list_based_on_Ca_over_entire_pdb(pdb+chain)
    for res in res_nums:	# go over residues that are relevant . BUT - maybe some have consurf score but no dist..(will be taken care below)
        print pdb,chain,res
        if res in list_of_res:
            dists,consurf_scores_wout_problems,consurf_scores_wout_problems_list=[],{},[]
        if res in dict_of_min_dist_from_res_list:
	    min_dist_from_res_list=dict_of_min_dist_from_res_list[res]
            for res2 in res_nums:
                if res2 in min_dist_from_res_list :
	            dists.append(float(min_dist_from_res_list[res2]))
	            consurf_scores_wout_problems[res2]=consurf_scores[res2] # only scores for res that have a pdb coordinates
                    consurf_scores_wout_problems_list.append(consurf_scores[res2])
            zipped_points=zip(consurf_scores_wout_problems_list,dists)
            for resi in zipped_points:
                print resi[0],'  ',resi[1]
            cor=pearsonr(consurf_scores_wout_problems_list,dists)
            print cor[0]
    return
def reduce_repeat_lines(f_inp,out):
    out=open(out,'w')
    f_inp=open(f_inp)
    lines=[]
    for line in f_inp:
        if line not in lines:
            lines.append(line)
    for line in lines:
        out.write(line)
    out.close()
    return

def Bootstrap_ave(list_of_val,num_of_bootstraps):
    # turn values of dict into list: 
    l=len(list_of_val)
    # take a bootstrapped sample
    counter=0
    list_of_aves_of_bootstraps=[]
    while counter<num_of_bootstraps:
        counter=counter+1
        bootstrapped_list_of_val=list(numpy.random.choice(list_of_val,l,replace=True))
    #    print bootstrapped_list_of_val
        ave_bootstrapped_sample=0.0
        for i in bootstrapped_list_of_val:
            ave_bootstrapped_sample=ave_bootstrapped_sample+float(i)
        ave_bootstrapped_sample=ave_bootstrapped_sample/l
        list_of_aves_of_bootstraps.append(ave_bootstrapped_sample)
  #  print list_of_aves_of_bootstraps
    final_average=numpy.mean(list_of_aves_of_bootstraps)
    final_std=numpy.std(list_of_aves_of_bootstraps)
    return num_of_bootstraps,final_average,final_std

def lin_func(x,a,b):
    return a*x+b
def lin_reg_w_yerr(x_array,y_array,y_err):
    coeffs, matcov = curve_fit(lin_func, x_array, y_array,p0=None,sigma=y_err)
    return coeffs, matcov
def calc_slope_w_yerr_from_dnds_file(dnds_file,x_array):
    dnds_file=open(dnds_file)
    y_array=[]
    y_error=[]
    l=len(x_array)
    for line in dnds_file:
        if line.split()[0]=='dN/dS':
            y_array.append(float(line.split()[1]))
            y_error.append(float(line.split()[4]))
    y_array=y_array[:l]
    y_error=y_error[:l]
    coeffs, matcov=lin_reg_w_yerr(x_array,y_array,y_error)
    perr = numpy.sqrt(numpy.diag(matcov))
    fit_array=[]
    for x in x_array:
        fit_array.append(coeffs[0]*x+coeffs[1])
    return coeffs, matcov,perr,fit_array
def calc_slope_no_yerr_from_dnds_file(dnds_file,x_array):
    dnds_file=open(dnds_file)
    y_array=[]
    y_error=None
    l=len(x_array)
    for line in dnds_file:
        if line.split()[0]=='dN/dS':
            y_array.append(float(line.split()[1]))
    y_array=y_array[:l]
    coeffs, matcov=lin_reg_w_yerr(x_array,y_array,y_error)
    perr = numpy.sqrt(numpy.diag(matcov))
    fit_array=[]
    for x in x_array:
        fit_array.append(coeffs[0]*x+coeffs[1])
    return coeffs, matcov,perr,fit_array

def parse_BioLip(BioLip_file):
    BioLip_info={}
    BioLip_pdb_to_uniprot={}
    enz_list=[]
    BioLip_file=open(BioLip_file)
    for line in BioLip_file:
        pdb=line.split('\t')[0]
	chain=line.split('\t')[1]
        BS=line.split('\t')[3]
        ligand=line.split('\t')[4]
        lig_chain=line.split('\t')[5]
	if pdb+chain not in BioLip_info:
	    BioLip_info[pdb+chain]={}
	binding_res=[]
	binding_data=line.split('\t')[7].split()
	ec=line.split('\t')[11]
	if ec!='':
	    enz_list.append(pdb+chain)
        uniprot_id=line.split('\t')[17]
	for i in binding_data:
	    binding_res.append(i[1:])
        if BS+ligand+lig_chain not in BioLip_info[pdb+chain]:
	    BioLip_info[pdb+chain][BS+ligand+lig_chain]=[]
            BioLip_info[pdb+chain][BS+ligand+lig_chain]=BioLip_info[pdb+chain][BS+ligand+lig_chain]+binding_res
        else:
            BioLip_info[pdb+chain][BS+ligand+lig_chain]=BioLip_info[pdb+chain][BS+ligand+lig_chain]+binding_res
        BioLip_pdb_to_uniprot[pdb+chain]=uniprot_id
    return BioLip_info,enz_list,BioLip_pdb_to_uniprot # BioLip_info[pdb+chain]={BS+ligand:binding res, BS+ligand 2: binding res 2,...}
def parse_BioLip_pubmed_id(BioLip_file):
    BioLip_info={}
    BioLip_file=open(BioLip_file)
    for line in BioLip_file:
        pdb,chain,ligand,binding_data,pubmed_id=line.split('\t')[0],line.split('\t')[1],line.split('\t')[4],line.split('\t')[7].split(),line.split('\t')[18]
        if pdb+chain not in BioLip_info:
            BioLip_info[pdb+chain]={}
            BioLip_info[pdb+chain][ligand]=pubmed_id
        else:
            BioLip_info[pdb+chain][ligand]=pubmed_id
    return BioLip_info

def parse_yeast_mapping_uniprot_geneNm():
    f=open('/home/avitaliv/projects/ctb-yxia/avitaliv//xialab/Data_files/YEAST_559292_idmapping.dat')
    yeastNm_to_uniprot={}
    for line in f:
        uniprot_id=line.split()[0]
        if 'Gene_OrderedLocusName' in line:
            yeast_nm=line.split()[2]
            yeastNm_to_uniprot[yeast_nm]=uniprot_id
    return yeastNm_to_uniprot

def correlate_two_aligned_properties(aligned_file1,aligned_file2):
    aligned1_properties_ofy,aligned1_seq_ofy=read_aligned_res_num(aligned_file1)
    aligned2_properties_ofy,aligned2_seq_ofy=read_aligned_res_num(aligned_file2)
    reduced_prop1=[]
    reduced_prop2=[]
    for y in aligned1_properties_ofy:
        if y in aligned2_properties_ofy:
            counter=-1
            for res1 in aligned1_properties_ofy[y]:
                counter=counter+1
                res2=aligned2_properties_ofy[y][counter]
                if '--' not in res1 and '--' not in res2:
                    reduced_prop1.append(float(res1))
                    reduced_prop2.append(float(res2))

    P=pearsonr(reduced_prop1,reduced_prop2)
    S=spearmanr(reduced_prop1,reduced_prop2)
    print 'pearson',P
    print 'spearman',S
    return

def reduce_aligned_accordY(inp_f,out_f):
    inp_f=open(inp_f)
    out=open(out_f,'w')
    Y_list=[]
    for line in inp_f:
        if line.split()[0] not in Y_list:
            Y_list.append(line.split()[0])
            out.write(line)
    out.close()
    return

def read_Moad_file():
    MOAD_valid_ligands={}
    MOAD_file=open('/home/avitaliv/projects/ctb-yxia/avitaliv//xialab/Data_files/MOAD_every.csv')
    for line in MOAD_file:
        if line[0]==',' and line.split(',')[2]!='':
            pdb= line.split(',')[2].lower()

        elif line[2]==',' and line.split(',')[4]=='valid':
            ligand=line.split(',')[3].split(':')[0]
            chain=line.split(',')[3].split(':')[1]
            if pdb+chain not in MOAD_valid_ligands:
                MOAD_valid_ligands[pdb+chain]=[]
                MOAD_valid_ligands[pdb+chain].append(ligand)
            else:
                MOAD_valid_ligands[pdb+chain].append(ligand)
    return MOAD_valid_ligands # MOAD_valid_ligands[pdb+chain]=list_of_ligands

def read_pdb_seq_clusters_95():
#    f=open('/home/avitaliv/projects/ctb-yxia/avitaliv//xialab/Data_files/pdb_sequence_clusters_90id.txt')
#    f=open('/home/avitaliv/projects/ctb-yxia/avitaliv//xialab/Data_files/pdb_sequence_clusters_70id.txt')
#    f=open('/home/avitaliv/projects/ctb-yxia/avitaliv//xialab/Data_files/pdb_sequence_clusters_50id.txt')
    f=open('/home/avitaliv/projects/ctb-yxia/avitaliv//xialab/Data_files/pdb_sequence_clusters_95id.txt')
    clusters=[]
    for line in f:
        cluster=[]
        for pdb in line.split():
            pdb=pdb.split('_')[0].lower()+pdb.split('_')[1]
            cluster.append(pdb)
        clusters.append(cluster)
    return clusters
def read_pdb_to_uniprot():
    f=open('/home/avitaliv/projects/ctb-yxia/avitaliv//xialab/Data_files/pdbtosp.txt')
    pdb_to_sp={}
    for line in f:
        if 'EM' in line or 'X-ray' in line or 'NMR' in line:
            pdb_to_sp[line.split()[0].lower()]=line[41:47].strip()
    return pdb_to_sp

def read_universal_genetic_code():
    f=open('/home/avitaliv/projects/ctb-yxia/avitaliv//singen_largerYdata2/data/universal_genetic_code.dat')
    nuc_aa_conversion={}
    for line in f:
        nuc_aa_conversion[line.split()[0]]=line.split()[1]
    return nuc_aa_conversion

def align_non_synonymous(codon_alignment_file,list_species,protein_nm):
    nuc_aa_conversion=read_universal_genetic_code()
    codon_alignment_file=open(codon_alignment_file)
    nucs={}
    for line in codon_alignment_file:
        if protein_nm in line:
            for species in list_species:
                if species in line and '.nt' in line:
                    nucs[species]=line.split()[2:]
    protein_non_synonymous=[]
    index_location=-1
    for nuc in nucs[list_species[0]]:
        index_location=index_location+1
        aa_prime=nuc_aa_conversion[nuc]
        non_synonymous=0
        if '--' in nuc:
            non_synonymous='---'
        else:
            for other_species in list_species[1:]:
                if '--' not in nucs[other_species][index_location]:
                    aa_other=nuc_aa_conversion[nucs[other_species][index_location]]
                else:
                    aa_other=aa_prime
                if aa_other!=aa_prime:
                    non_synonymous=1
            protein_non_synonymous.append(non_synonymous)
    return protein_non_synonymous # list of non-synonymouse (0 or 1) for each nuc of the specified yeast protein  
def calc_id_and_coverage_of_hits(hits_file):
    id_cov={}
    f=open(hits_file)
    for line in f:
       qseqid,sseqid,qlen,slen,qstart,qend,sstart,send,nident=line.split()[0],line.split()[1],line.split()[2],line.split()[3],line.split()[5],line.split()[6],line.split()[7],line.split()[8],line.split()[9]
       scov=100*(float(send)-float(sstart))/float(slen)
       qcov=100*(float(qend)-float(qstart))/float(qlen)
       sid=100*float(nident)/float(slen)
       qid=100*float(nident)/float(qlen)
       id_cov[(qseqid,sseqid)]=[sid,qid,scov,qcov]
    return id_cov

def read_funct_res(list_file):
    list_file=open(list_file)
    funct_residues={}
    list_line_identifiers=[]
    for line in list_file:
        yeast,pdb,chain=line.split()[0],line.split()[1],line.split()[2]
        if yeast+pdb+chain not in list_line_identifiers:
            counter=1
            list_line_identifiers.append(yeast+pdb+chain)
            funct_residues[yeast+pdb+chain]={}
        else:
            list_line_identifiers.append(yeast+pdb+chain)
            counter=list_line_identifiers.count(yeast+pdb+chain)
        if (line.split()[3]).isalpha()==True:
            lig=line.split()[3]
            funct_ress=line.split()[4:]
        else:
            funct_ress=line.split()[3:]
        funct_residues[yeast+pdb+chain][counter]=funct_ress
    return funct_residues

def remove_duplicate_lines(file_inp,out):
    file_inp=open(file_inp)
    out=open(out,'w')
    list_lines=[]
    for line in file_inp:
        if line not in list_lines:
            list_lines.append(line)
            out.write(line)
    out.close()
    return

def read_pdb_seq_clusters_90():
#    f=open('/home/avitaliv/projects/ctb-yxia/avitaliv/xialab/Data_files/pdb_sequence_clusters_90id.txt')
 #   f=open('/home/avitaliv/projects/ctb-yxia/avitaliv/xialab/Data_files/pdb_sequence_clusters_70id.txt')
    f=open('/home/avitaliv/projects/ctb-yxia/avitaliv/xialab/Data_files/pdb_sequence_clusters_50id.txt')
    clusters=[]
    for line in f:
        cluster=[]
        for pdb in line.split():
            pdb=pdb.split('_')[0].lower()+'_'+pdb.split('_')[1]
            cluster.append(pdb)
        clusters.append(cluster)
    return clusters

def read_pdb_seq_clusters_100():
#    f=open('/home/avitaliv/projects/ctb-yxia/avitaliv/xialab/Data_files/pdb_sequence_clusters_90id.txt')
 #   f=open('/home/avitaliv/projects/ctb-yxia/avitaliv/xialab/Data_files/pdb_sequence_clusters_70id.txt')
    f=open('/home/avitaliv/projects/ctb-yxia/avitaliv/xialab/Data_files/pdb_sequence_clusters_100id.txt')
    clusters=[]
    for line in f:
        cluster=[]
        for pdb in line.split():
            pdb=pdb.split('_')[0].lower()+'_'+pdb.split('_')[1]
            cluster.append(pdb)
        clusters.append(cluster)
    return clusters

def read_pdb_seq_clusters_30():
    f=open('/home/avitaliv/projects/ctb-yxia/avitaliv/xialab/Data_files/pdb_sequence_clusters_30id.txt')
    clusters=[]
    for line in f:
        cluster=[]
        for pdb in line.split():
            pdb=pdb.split('_')[0].lower()+'_'+pdb.split('_')[1]
            cluster.append(pdb)
        clusters.append(cluster)
    return clusters
def read_ec_pdbs():
    f=open('/home/avitaliv/projects/ctb-yxia/avitaliv/xialab/Data_files/ec.csv')
    list_enz=[]
    for line in f:
        pdb=line.split()[0].lower()
        for chain in line.split()[1:]:
            list_enz.append(pdb+chain)
    return list_enz

def calc_WCNSC_of_pdb(pdb,chain,pdb_path):
    res_atoms_properties={}
    mass={'C':12,'O':16,'N':14,'S':32,'Se':79,'H':1}
    pdb_path=pdb_path+pdb+'.pdb.gz'
    with gzip.open(pdb_path,'r') as f:
        lines=f.readlines()
    residues_used=[]
    for line in lines:
        if len(line)>19 and line.split()[0]=='ATOM' and line[21]==chain:
        #    print line[30:37],line[38:45],line[46:53],line[13],line[13:17].strip(),line[22:26].strip(),line[17:20]
            x,y,z,m,atom_type,res_num,res_type=float(line[30:37]),float(line[38:45]),float(line[46:53]),float(mass[line[12:17].strip()[0]]),line[13:16].strip(),line[22:26].strip(),line[17:20]
            if res_num not in res_atoms_properties:
                res_atoms_properties[res_num]={}
                res_atoms_properties[res_num][atom_type]=(x,y,z,m)
            else:
                res_atoms_properties[res_num][atom_type]=(x,y,z,m)
            if res_type=='GLY' and atom_type=='CA':
                res_atoms_properties[res_num]['gly']=(x,y,z,1.0)
    # correction - remove residues without a CA
    corrected_res_atoms_properties={}
    for res in res_atoms_properties:
        if 'CA' in res_atoms_properties[res]:
            corrected_res_atoms_properties[res]=res_atoms_properties[res]
    res_atoms_properties=corrected_res_atoms_properties
    CMx_of_res,CMy_of_res,CMz_of_res={},{},{}
    WCNSC_of_res={}
    for res in res_atoms_properties:
        sc_xs,sc_ys,sc_zs,sc_ms=[],[],[],[]
        for atom_type in res_atoms_properties[res]:
            if atom_type!='CA' and atom_type!='C' and atom_type!='O' and atom_type!='N' and atom_type!='H' and atom_type!='D':
                sc_xs.append(float(res_atoms_properties[res][atom_type][0]))
                sc_ys.append(float(res_atoms_properties[res][atom_type][1]))
                sc_zs.append(float(res_atoms_properties[res][atom_type][2]))
                sc_ms.append(float(res_atoms_properties[res][atom_type][3]))
        if sc_xs!=[]: # make sure to only take into account residues with side chains (otherwise will have 'nan')
            CMx_of_res[res]=numpy.sum(numpy.array(sc_xs)*numpy.array(sc_ms))/numpy.sum(numpy.array(sc_ms)) 
            CMy_of_res[res]=numpy.sum(numpy.array(sc_ys)*numpy.array(sc_ms))/numpy.sum(numpy.array(sc_ms))
            CMz_of_res[res]=numpy.sum(numpy.array(sc_zs)*numpy.array(sc_ms))/numpy.sum(numpy.array(sc_ms))
    for res in CMx_of_res:  # only for residues that have a side chain we can calculate the wcnsc 
        WCNSC_of_res[res]=0.0
        for res_temp in res_atoms_properties:
            if res!=res_temp:
                term1=1/((CMx_of_res[res]-res_atoms_properties[res_temp]['CA'][0])**2+(CMy_of_res[res]-res_atoms_properties[res_temp]['CA'][1])**2+(CMz_of_res[res]-res_atoms_properties[res_temp]['CA'][2])**2)
                if res_temp in CMx_of_res and 'gly' not in res_atoms_properties[res_temp]: # only temp_residues that have a side chain will contribute the term2 part as well 
                    term2=1/((CMx_of_res[res]-CMx_of_res[res_temp])**2+(CMy_of_res[res]-CMy_of_res[res_temp])**2+(CMz_of_res[res]-CMz_of_res[res_temp])**2)
                else:
                    term2=0.0
                WCNSC_of_res[res]=WCNSC_of_res[res]+term1+term2
    return WCNSC_of_res
