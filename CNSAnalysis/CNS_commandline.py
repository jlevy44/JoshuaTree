#######################################################################
############################# IMPORTS #################################

import begin, ete3
from ete3 import PhyloTree,Tree,TreeStyle,NodeStyle,EvolTree
import subprocess, os
import re
from Bio import Phylo
from Bio.Phylo import Consensus as CS
import dill as pickle
import numpy as np
import sys
import pandas as pd
from collections import Counter
from bx.align import maf
from bx import interval_index_file
from itertools import combinations


############################################################################################
############################# PARSE CNS CONFIG AND CNS RUN #################################

def parseConfigFindList(stringFind,configFile):
    """parseConfigFindList inputs a particular string to find and read file after and a configuration file object
    outputs list of relevant filenames"""
    read = 0
    listOfItems = []
    for line in configFile:
        if line:
            if read == 1:
                if 'Stop' in line:
                    configFile.seek(0)
                    break # exit the function and return the list of files or list information
                listOfItems.append(line.strip('\n'))
            if stringFind in line:
                read = 1 # if find string specified, begin reading lines
    configFile.seek(0)
    return listOfItems

@begin.subcommand
def runCNSAnalysis():
    subprocess.call('python runCNSAnalysis.py')


######################################################################################################
############################# COMPARE CNS TO DIFF EXPRESSION RESULTS #################################

@begin.subcommand
def CNSvsDE(DE_genes, tree_file, shortname2name, main_species, CNS_bed):#all_genes,
    from scipy.stats import chi2_contingency
    #all_genes = np.vectorize(lambda x: x.strip('"'))(os.popen("awk -F'[,=]' '{print $1}' %s"%all_genes).read().splitlines()[1:])
    DE_genes = np.vectorize(lambda x: x.strip('"'))(os.popen("awk -F'[,=]' '{print $1}' %s"%DE_genes).read().splitlines()[1:])
    #non_DE_genes = np.setdiff1d(all_genes,DE_genes)
    if shortname2name.endswith('.txt'):
        with open(shortname2name,'r') as f:
            shortname2name = dict([tuple(line.split()) for line in f.read().splitlines() if line])
    else:
        shortname2name = {shortname:name for shortname,name in [shortname_dict.split(':') for shortname_dict in shortname2name.split(',')]}
    dist_mat = tree2matrix(tree_file).rename(index=shortname2name,columns=shortname2name)
    final_array = []
    with open(CNS_bed,'r') as f:
        for line in f:
            if line:
                ll = line.split('closestGene=')
                MASegs = ll[0].split()[-1].split('|')
                ll2 = ll[1].split(';')
                genes = ll2[0].split('|')
                distance = ll2[1].replace('distance=','')
                for MASeg in MASegs:
                    #print MASeg
                    species_cons = dict([species.split(':') for species in MASeg.split(';')[1].split(',')])
                    present_species = [species for species in species_cons if int(species_cons[species])]
                    conservation_score = np.sum(np.vectorize(lambda x: dist_mat[x][main_species] if x != main_species else 0)(present_species))
                    for gene in genes:
                        final_array.append([gene in DE_genes,distance,conservation_score])
    df = pd.DataFrame(final_array,columns=['DE_GENE','Distance','Conservation_Score'])
    print df
    df.to_csv('CNSvsDE_raw.csv')
    distanceGroups_keys = ['Present Within','Present Within 0-100','Present Within 100-1000','Present Outside 1000']
    distances = [(0,),(0,100),(100,1000),(1000,1000000)]
    final_table_keys = ['DE_Gene','No_DE_Gene']
    final_table = {key:{distance_group:0 for distance_group in distanceGroups_keys} for key in final_table_keys}

    for i,DE in [(False,'No_DE_Gene'),(True,'DE_Gene')]:
        df2 = df[df['DE_GENE']==i]
        #print df2
        #print sum(df2['Distance'].as_matrix().astype(np.int) > 0)
        for j in range(len(distances)):
            print DE,distanceGroups_keys[j]
            if len(distances[j]) == 1:
                #print 'A'
                #print sum(df2['Distance'].as_matrix().astype(np.int) == distances[j][0])
                final_table[DE]['Present Within'] = sum(df2['Distance'].as_matrix().astype(np.int) == distances[j][0])
            else:
                #print 'B'
                #print sum(((df2['Distance'].as_matrix().astype(np.int) > distances[j][0])&(df2['Distance'].as_matrix().astype(np.int) <= distances[j][1])))
                final_table[DE][distanceGroups_keys[j]] = sum(((df2['Distance'].as_matrix().astype(np.int) > distances[j][0])&(df2['Distance'].as_matrix().astype(np.int) <= distances[j][1])))
        print final_table
    df = pd.DataFrame(final_table)
    print df
    df.to_csv('CNSvsDE_final.csv',index=False)
    chi_sq = chi2_contingency(df.as_matrix())
    print chi_sq
    pd.DataFrame(np.array(chi_sq)).to_csv('CNSvsDE_final_chisq.csv',index=False)
    # FIXME DEBUG ABOVE, SO CLOSE!!!

#################################################################################################
############################## CONVERT TREE TO DISTANCE MATRIX ##################################

@begin.subcommand
def tree2matrix(tree_file,distance_matrix = 'distance_matrix.csv'):
    tree = Phylo.read(tree_file,'newick')
    allclades = list(tree.find_clades(order='level'))
    species_names = [clade.name for clade in allclades if clade.name]
    df = pd.DataFrame(np.nan, index=species_names, columns=species_names)
    for i,j in combinations(species_names,r=2):
        if i == j:
            df.set_value(i,j,0)
        if i != j:
            distance = tree.distance(i,j)
            df.set_value(i,j,distance)
            df.set_value(j,i,distance)
    df.to_csv(distance_matrix)
    return df

@begin.subcommand
def replace_species_names(csv_file, reverse_lookup_dict):
    with open(csv_file,'r') as f:
        txt = f.read()
        'FINISH!!'


#############################################################################################
############################## ETE3 AND BIOPYTHON ANALYSES ##################################

@begin.subcommand
def output_Tree(aln_output_txt,aln_file, out_fname):
    if aln_file.endswith('.phylip'):
        print 'Input must be fasta file for now'
        quit()
    elif aln_file.endswith('.fasta') or aln_file.endswith('.fa'):
        subprocess.call("awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' %s > aln.fasta"%aln_file,shell=True)
        t = PhyloTree(aln_output_txt,alignment='aln.fasta',alg_format='fasta')
    else:
        t = Tree(aln_output_txt)
    ts = TreeStyle()
    ns = NodeStyle()
    ns['size']=0
    ts.show_leaf_name = True
    ts.show_branch_length = False
    ts.show_branch_support = True
    for n in t.traverse():
        n.set_style(ns)
    #t.show(tree_style=ts)
    t.render(out_fname,tree_style = ts)

@begin.subcommand()
def find_tree_lengths(tree_file):
    trees_gen = Phylo.parse(tree_file, 'newick')
    output_tree_lengths = []
    with open(tree_file,'r') as f:
        for i in range(len(f.readlines())):
            try:
                tree = trees_gen.next()
                output_tree_lengths.append(tree.total_branch_length())
            except:
                output_tree_lengths.append(0)
    np.save('out_tree_lengths.npy',np.array(output_tree_lengths))
    return output_tree_lengths

@begin.subcommand
def evaluate_selective_pressure(maf_structure_pickle, neutral_tree_nwk):
    # FIXME codeml only works in protein coding regions https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1470900/ need to add for CNS!!! https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3725466/
    # FIXME branch specific rate changes??
    # FIXME also degree of conservation when looking at selective pressure density
    # FIXME add new metric? Conserved pressure metric, conserved ratio (# conserved/total conserved weighted) * selective pressure density?? for different elements??
    codeMLPath = next(path for path in sys.path if 'conda/' in path and '/lib/' in path).split('/lib/')[0]+'/bin/ete3_apps/bin'
    try:
        os.mkdir('./selectivePressureTests')
    except:
        pass
    maf = pickle.load(open('MAFWork.p','rb'))
    onesCondition = maf[3].keys()[np.argmax([sum(val.values()) if type(val) != type([]) else 0 for val in maf[3].values()])]
    for maseg in maf[1][onesCondition]['CS']:
        model_tree = EvolTree(open(neutral_tree_nwk,'rb').read(),binpath=codeMLPath) #FIXME is this the right tree to be using???
        aln_lines = maf[1][onesCondition]['CS'][maseg].splitlines()
        coord_info = {}
        for i in range(len(aln_lines))[::2]:
            coord_info[aln_lines[i].find('>')+1:aln_lines[i].find('.')] = aln_lines[i][aln_lines[i].find('.'):] #FIXME output these coordinates to bed file
            aln_lines[i] = aln_lines[i][:aln_lines[i].find('.')]
            alignment = '\n'.join(aln_lines)
        model_tree.link_to_alignment(alignment)
        model_tree.workdir = './selectivePressureTests'
        model_tree.run_model('M1')
        model_tree.run_model('M2')
        pval = model_tree.get_most_likely('M2','M1')
        if pval < 0.05:
            for s in range(len(model2.sites['BEB']['aa'])):
                print 'positively selected site %s at position: %s, with probability: %s' % (model2.sites['BEB']['aa'][s], s+1, model2.sites['BEB']['p2'][s])
        # FIXME first fix codeml, then print out selective pressure density after removing --, append that to masegment, bed file with coordinates, maSeg, and selective pressure density, or could just output one-base coordinate and calculate density through bedtools
        else:
            print 'fail to reject null hypothesis'
        """
        best_model = None
        best_lnl = float('-inf')
        for starting_omega in [0.2, 0.7, 1.2]: # FIXME optimize and choose parameters
            model_tree.run_model('b_free.'+str(starting_omega))
            current_model = model_treetree.get_evol_model('b_free.'+str(starting_omega))
            if current_model.lnL > best_lnl:
                best_lnl = current_model.lnL
                best_model = current_model
        """
        #print best_model
        tree.render(maseg+'.png')
        break # FIXME remove after testing
        "FIXME UNDER DEVELOPMENT"

@begin.subcommand
def reroot_Tree(tree_in,root_species,tree_out):
    t = PhyloTree(open(tree_in,'r').read())
    t.set_outgroup(root_species)
    t.write(tree_out)

@begin.subcommand
def run_phyml(phylip_in,bootstrap):
    if phylip_in.endswith('.fasta') or phylip_in.endswith('.fa'):
        subprocess.call(['perl', 'Fasta2Phylip.pl', phylip_in, phylip_in.replace('fasta','phylip')])
        phylip_in = phylip_in.replace('fasta','phylip')
    phymlLine = next(path for path in sys.path if 'conda/' in path and '/lib/' in path).split('/lib/')[0]+'/bin/ete3_apps/bin/phyml'
    subprocess.call([phymlLine, '-i', phylip_in, '-s', 'BEST', '-q', '-b', bootstrap, '-m', 'GTR'])

@begin.subcommand
def omega_analysis(aln_syn,aln_non_syn,phyml_fasttree, reference, list_species, build_tree):
    #run_phyml(aln_syn,'1')
    if int(build_tree):
        subprocess.call('ete3 build -n %s -o syn_tree/ --clearall -w none-none-none-%s_default'%(aln_syn,phyml_fasttree),shell=True)
        subprocess.call('ete3 build -n %s -o non_syn_tree/ --clearall -w none-none-none-%s_default'%(aln_non_syn,phyml_fasttree),shell=True)

    # get distances between reference and list of species
    tree2matrix('syn_tree/none-none-none-%s_default/%s.final_tree.nw'%(phyml_fasttree,aln_syn),'synonymous_dist_matrix.csv')
    tree2matrix('non_syn_tree/none-none-none-%s_default/%s.final_tree.nw'%(phyml_fasttree,aln_non_syn),'non-synonymous_dist_matrix.csv')
    list_species = list_species.split(',')
    syn_dist = pd.read_csv('synonymous_dist_matrix.csv',index_col=[0])
    non_syn_dist = pd.read_csv('non-synonymous_dist_matrix.csv',index_col=[0])
    df = pd.DataFrame(index=[aln_type+reference for aln_type in ['syn:','non-syn:']],columns=list_species)
    for species in list_species:
        df.set_value('syn:'+reference,species,syn_dist[reference][species])
        df.set_value('non-syn:'+reference,species,non_syn_dist[reference][species])
    df.to_csv('syn_non-syn_comparison_%s.csv'%reference)

    # FIXME Find way to compare two models
    neutral_tree = EvolTree('syn_tree/none-none-none-%s_default/%s.final_tree.nw'%(phyml_fasttree,aln_syn))
    with open(aln_syn,'r') as f:
        neutral_tree.link_to_alignment(f.read())
    non_syn_tree = EvolTree('non_syn_tree/none-none-none-%s_default/%s.final_tree.nw'%(phyml_fasttree,aln_non_syn))
    with open(aln_non_syn,'r') as f:
        non_syn_tree.link_to_alignment(f.read())
    neutral_tree.run_model('fb')
    non_syn_tree.run_model('fb')
    print neutral_tree.get_evol_model('fb')
    print non_syn_tree.get_evol_model('fb')

#########################################################################
############################ VCF FILTER WORK ############################

@begin.subcommand
def intersect_vcf(vcf_in, bed_regions, vcf_out):
    subprocess.call('bedtools intersect -wa -a %s -b %s > %s'%(vcf_in, bed_regions, vcf_out),shell=True)

@begin.subcommand
def annotate_snps(vcf_in,genome_name,fasta_in,gff_in,vcf_out):
    """Goal is to find synonymous and non-synonymous regions of SNPs"""
    # format gff in
    #subprocess.call('grep -v "#" %s | sort -k1,1 -k2,2n -k3,3n -t$\'\t\' | bgzip -c > %s.gz && tabix -p gff %s.gz'%(gff_in,gff_in,gff_in),shell=True)
    snp_eff_line = next(path for path in sys.path if 'conda/' in path and '/lib/' in path).split('/lib/')[0]+'/share/snpeff-4.3.1r-0/snpEff.config'
    subprocess.call('scp %s . && mkdir -p ./data/%s'%(snp_eff_line,genome_name),shell=True)
    with open('snpEff.config','r') as f:
        write = 0
        writeLines = []
        for line in f:
            if 'Databases & Genomes' in line:
                write = 1
            writeLines.append(line)
            if write and '#-----------------------------' in line:
                writeLines.append('\n# my %s genome\n%s.genome : %s\n'%(genome_name,genome_name,genome_name))
                write = 0
    with open('snpEff.config','w') as f:
        f.writelines(writeLines)
    subprocess.call('scp %s ./data/%s/sequences.fa && scp %s ./data/%s/genes.gff'%(fasta_in,genome_name,gff_in,genome_name),shell=True)
    subprocess.call('snpEff build -c snpEff.config -gff3 -v %s > snpEff.stdout 2> snpEff.stderr'%genome_name,shell=True)
    subprocess.call('snpEff -c snpEff.config %s %s > %s'%(genome_name,vcf_in,vcf_out),shell=True)

@begin.subcommand
def generate_N_S_aln(reference_species,vcf_in):
    with open(vcf_in,'r') as f:
        for line in f:
            if line.startswith('#CHROM'):
                species_list = line.split()[9:]
                break
        aln_dict = {site_type: {species : '' for species in species_list} for site_type in ['synonymous','non-synonymous']}
    #print aln_dict
    for site_type, f in zip(['synonymous','non-synonymous'],[os.popen('grep synonymous %s'%vcf_in),os.popen('grep protein_coding %s | grep -v synonymous'%vcf_in)]):
        lines = f.read().splitlines()
        for line in lines:
            ll = line.split()
            #print ll
            #aln_dict[site_type][reference_species] += ll[3]
            snp_dict = dict({'0':ll[3],'.':'-'},**{str(snp+1): snp_chr for snp,snp_chr in enumerate(ll[4].split(','))})
            #print snp_dict
            for idx, snp in enumerate(ll[9:]):
                aln_dict[site_type][species_list[idx]] += snp_dict[snp]
        with open(site_type+'_regions_aln_%s.fasta'%reference_species,'w') as f:
            f.write('\n'.join(['>%s\n%s'%(species,aln_dict[site_type][species]) for species in aln_dict[site_type]]))


def change_of_coordinates(in_file,out_file): # FIXME USE PYTABIX FIXES NEEDED, maffilter change coordinates??
    with open(in_file,'r') as f, open(out_file,'w') as f2:
        for line in f:
            if line.startswith('#') == 0:
                break
            else:
                f2.write(line)
        #    offset = f.tell()
        #f.seek(offset) #FIXME TEST
        for line in f:
            lineList = line.split()
            lineList[1] = str(int(lineList[0].split('_')[-2]) + int(lineList[1]))
            lineList[0] = lineList[0][lineList[0].find('.')+1:[m.start(0) for m in re.finditer('_',lineList[0])][-2]]
            f2.write('\t'.join(lineList)+'\n')

def check_vcf_empty(vcf_in):
    with open(vcf_in,'r') as f:
        for line in f:
            if line.startswith('#') == 0:
                break
            offset = f.tell()
        f.seek(offset)
        if f.readline():
            return False
        else:
            return True

def sort_vcf(vcf_in,vcf_out):
    subprocess.call("""bgzip -c {0} > vcfs/out.vcf.gz
                (zcat vcfs/out.vcf.gz | head -300 | grep ^#;
                zcat vcfs/out.vcf.gz | grep -v ^# | sort -k1,1d -k2,2n;) \
                | bgzip -c > {1}.gz
                bcftools index {1}.gz""".format(vcf_in,vcf_out), shell=True)

@begin.subcommand
def concat_vcf(list_vcfs,vcf_out):#,vcf_out):
    list_vcfs = list_vcfs.split(',')
    master_df = pd.DataFrame()
    header_lines = []
    for vcf_in in list_vcfs:
        with os.popen('zcat %s'%vcf_in) as f:
            for line in f:
                header_lines.append(line)
                if line.startswith('#CHROM'):
                    line_info = line.strip('/n').split()
                    break
        master_df = master_df.append(pd.DataFrame(np.hstack([np.array(os.popen("zcat %s | grep -v ^# | awk '{ print $%d }'"%(vcf_in,i+1)).read().splitlines())[:,None] for i in range(len(line_info))]),columns = line_info))
    header_lines = set(header_lines)
    master_df = master_df.sort_values(['#CHROM','POS'])
    #print master_df
    master_df.to_csv(vcf_out,sep='\t',index=False, na_rep = '.')
    with open(vcf_out.replace('.vcf','.headers.vcf'),'w') as f, open(vcf_out,'r') as f2:
        for line in [line2 for line2 in header_lines if '#CHROM' not in line2]:
            f.write(line)
        f.write(f2.read())
    subprocess.call('mv %s %s'%(vcf_out.replace('.vcf','.headers.vcf'),vcf_out),shell=True)
    #print master_df

#########################################################################
############################ MAF FILTER WORK ############################


@begin.subcommand
def maf2vcf(cns_config, reference_species, change_coordinates, out_all_species, overlaps):
    change_coordinates = int(change_coordinates)
    out_all_species = int(out_all_species)
    overlaps = int(overlaps)
    mafFiles = [file for file in os.listdir('.') if file.endswith('.maf') and file.startswith('FastaOut')]
    try:
        os.mkdir('vcfs')
    except:
        pass
    if cns_config.endswith('.txt'):
        with open(cns_config,'r') as f:
            master_species = parseConfigFindList("masterListSpecies",f)
    else:
        master_species = cns_config.split(',')
    all_species = set([species.split('_')[0] for species in master_species])
    all_species_but_one = all_species - {reference_species}
    if out_all_species:
        species = all_species
    else:
        species = all_species_but_one
    finalOutVCFFiles = []
    for i,maf in enumerate(mafFiles):
        subprocess.call("sed '/Anc/d' %s > out_maf.maf"%maf,shell=True)
        with open('maf_filter_config.bpp','w') as f:
            f.write("""
    input.file=./out_maf.maf
    input.format=Maf
    output.log=out.log
    maf.filter=\\
        Subset(\\
                strict=yes,\\
                keep=no,\\
                species=(%s),\\
                remove_duplicates=yes),\\
        Subset(\\
                strict=yes,\\
                keep=yes,\\
                species=(%s),\\
                remove_duplicates=yes),\\
        VcfOutput(\\
                file=vcfs/Out%d_all.vcf,\\
                genotypes=(%s),\\
                all=no,\\
                reference=%s)
            """%(','.join(all_species),reference_species,i,','.join(species),reference_species))
        subprocess.call('./maffilter param=maf_filter_config.bpp',shell=True)
        if check_vcf_empty('vcfs/Out%d_all.vcf'%i) == 0:
        # see if there is anything in file
        # if yes, change coordinates and sort
            finalOutVCFFiles.append('vcfs/Out%d_all.vcf.gz'%i)
            subprocess.call("sed 's/##source=Bio++/##INFO=<ID=AC>/g' %s > vcfs/temp.vcf && mv vcfs/temp.vcf %s && rm vcfs/temp.vcf"%('vcfs/Out%d_all.vcf'%i,'vcfs/Out%d_all.vcf'%i),shell=True)
            change_file = 'vcfs/Out%d_all.vcf'%i
            if change_coordinates:
                change_of_coordinates('vcfs/Out%d_all.vcf'%i,'vcfs/Out%d_new_coord_unsorted.vcf'%i)
                change_file = 'vcfs/Out%d_new_coord_unsorted.vcf'%i
            sort_vcf(change_file,'vcfs/Out%d_all.vcf'%i)
            #tabix -p vcf %s.gz
        # change the coordinate system and sort the file, also remove empty vcf files, fix merge/concat
    # FIXME problem with overlaps, --allow overlaps works with some analyses but not others, may want to just throw in contig list!!! below
    concat_vcf(','.join(finalOutVCFFiles),'vcfs/final_all_sorted.vcf')
    #subprocess.call('bcftools concat%s -O v -o vcfs/final_all.vcf %s'%(' --allow-overlaps' if overlaps else '',' '.join(finalOutVCFFiles)),shell=True)
    #subprocess.call("sed 's/##source=Bio++/##INFO=<ID=AC>/g' vcfs/final_all.vcf > vcfs/final_all_edit.vcf && mv vcfs/final_all_edit.vcf vcfs/final_all.vcf && rm vcfs/final_all_edit.vcf",shell=True)
    #sort_vcf('vcfs/final_all.vcf','vcfs/final_all_sorted.vcf')
    #subprocess.call('rm vcfs/final_all_sorted.vcf',shell=True)
    #subprocess.call('zcat vcfs/final_all_sorted.vcf.gz > vcfs/final_all_sorted.vcf',shell=True)

    #subprocess.call('bcftools sort -o vcfs/final_all_sorted.vcf vcfs/final_all.vcf',shell=True)
    #FIXME add sort function

@begin.subcommand
def index_maf(maf_file, ref_species):
    index_file = maf_file + '.idx'
    indexes = interval_index_file.Indexes()
    with open(maf_file,'r') as f:
        reader = maf.Reader(f)
        while True:
            pos = reader.file.tell()
            rec = reader.next()
            if rec is None:
                break
            for c in rec.components:
                indexes.add(c.src,c.forward_strand_start,c.forward_strand_end, pos)
    with open(index_file,'w') as f:
        indexes.write(f)

@begin.subcommand
def index_maf_2(maf_file):
    idxs = []
    with open(maf_file,'r') as f:
        offset = 0
        for line in f:
            if line.startswith('a'):
                idxs.append(offset)
            offset = f.tell()
        idxs.append(f.tell())
    return {idxs[i]:idxs[i+1] for i in range(len(idxs)-1)}

@begin.subcommand
def estimate_phylogeny(cns_config, consensus_algorithm, major_cutoff, min_block_length, concat_size, consensus_tree, bootstrapped_distance, feature_file, reference_species, feature_type):
    min_block_length = int(min_block_length)
    consensus_tree = int(consensus_tree)
    concat_size = int(concat_size)
    bootstrapped_distance = int(bootstrapped_distance)
    try:
        os.mkdir('./maf_trees')
    except:
        pass
    mafFiles = [file for file in os.listdir('.') if file.endswith('.maf') and (file.startswith('merged') + file.startswith('out_maf') == 0)]
    try:
        os.mkdir('vcfs')
    except:
        pass
    with open(cns_config,'r') as f:
        master_species = parseConfigFindList("masterListSpecies",f)
    all_species = set([species.split('_')[0] for species in master_species])
    if consensus_tree:
        fileOutTrees = []
        for i,maf in enumerate(mafFiles):
            subprocess.call("sed '/Anc/d' %s > out_maf.maf"%maf,shell=True)
            with open('maf_filter_config.bpp','w') as f:
                f.write("""
        input.file=./out_maf.maf
        input.format=Maf
        output.log=out.log
        maf.filter=\\
            Subset(\\
                    strict=yes,\\
                    keep=no,\\
                    species=(%s),\\
                    remove_duplicates=yes),\\
            MaskFilter(species=(%s)),\\
            MinBlockLength(min_length=%d), \\
            Merge(species=(%s)),\\
            Concatenate(minimum_size=%d),\\
            RemoveEmptySequences(),\\
            DistanceEstimation(\\
                    method=ml,\\
                    model=GTR,\\
                    gap_option=no_gap,\\
                    parameter_estimation=initial,\\
                    gaps_as_unresolved=no,\\
                    unresolved_as_gap=yes,\\
                    extended_names=yes),\\
            DistanceBasedPhylogeny(\\
                    method=bionj,\\
                    dist_mat=MLDistance),\\
            OutputTrees(\\
                    tree=BioNJ,\\
                    file=./maf_trees/trees_%d.nwk,\\
                    compression=none,\\
                    strip_names=yes)
                """%(','.join(all_species),','.join(all_species),min_block_length,','.join(all_species),concat_size,i))
            subprocess.call('./maffilter param=maf_filter_config.bpp',shell=True)
            with open('./maf_trees/trees_%d.nwk'%i,'r') as f:
                if f.readline():
                    fileOutTrees.append('./maf_trees/trees_%d.nwk'%i)
        with open('./maf_trees/final_trees.nwk','w') as f: #FIXME maf filter throws exception while calculating distance matrix at empty sites
            for tree in fileOutTrees:
                with open(tree,'r') as f2:
                    for line in f2:
                        if line and line.count('(') == line.count(')'):
                            f.write(line)
        #subprocess.call('for f in ./maf_trees/trees_*.nwk; do cat "$f"; echo "\\newline"; done > out && mv out ./maf_trees/final_trees.nwk',shell=True)
        trees_gen = Phylo.parse('./maf_trees/final_trees.nwk', 'newick')
        trees = []
        with open('./maf_trees/final_trees.nwk','r') as f:
            for i in range(len(f.readlines())):
                try:
                    trees.append(trees_gen.next())
                except:
                    pass
        if consensus_algorithm == 'strict':
            tree = CS.strict_consensus(trees)
        elif consensus_algorithm == 'majority':
            tree = CS.majority_consensus(trees,float(major_cutoff))
        else:
            tree = CS.adam_consensus(trees)
        Phylo.write(tree,'./maf_trees/output_tree_consensus.nh','newick')
    else:
        subprocess.call('rm merged.maf',shell=True)
        for file in mafFiles:
            subprocess.call("sed -e '/Anc/d;/#/d' %s >> merged.maf"%file,shell=True)
        subprocess.call("cat -s merged.maf > temp.maf && mv temp.maf merged.maf",shell=True)
        #with open('merged.maf','r') as f: #FIXME memory errors
        #    txt = f.read().replace('\n\n\n','\n\n')
        #with open('merged.maf','w') as f:
        #    f.write(txt)
        #del txt
        #subprocess.call("sed '/#/d' merged.maf > out_maf.maf",shell=True)
        with open('maf_filter_config.bpp','w') as f:
                f.write("""
        input.file=./merged.maf
        input.format=Maf
        output.log=out.log
        maf.filter=\\
            Subset(\\
                    strict=yes,\\
                    keep=no,\\
                    species=(%s),\\
                    remove_duplicates=yes),\\
            MaskFilter(species=(%s)),\\
            MinBlockLength(min_length=%d), \\%s
            Merge(species=(%s)),\\
            Concatenate(minimum_size=1000000000000),\\%s
            DistanceEstimation(\\
                    method=ml,\\
                    model=GTR,\\
                    gap_option=no_gap,\\
                    parameter_estimation=initial,\\
                    gaps_as_unresolved=no,\\
                    unresolved_as_gaps=yes,\\
                    extended_names=yes),\\%s
            DistanceBasedPhylogeny(\\
                    method=bionj,\\
                    dist_mat=MLDistance),\\
            OutputTrees(\\
                    tree=BioNJ,\\
                    file=./maf_trees/output_tree_consensus.nh,\\
                    compression=none,\\
                    strip_names=yes)
                """%(','.join(all_species),','.join(all_species),min_block_length,
                     ("""ExtractFeature(\\
                            ref_species=%s,\\
                            feature.file=%s,\\
                            feature.format=%s,\\
                            feature.type=%s,\\
                            feature.file.compression=none,\\
                            complete=no,\\
                            compression=none),\\"""%(reference_species,feature_file,'BedGraph' if feature_file.endswith('.bed')
                                        or feature_file.endswith('.bed3') else 'GFF', feature_type) if feature_file.endswith('.bed') or
                                        feature_file.endswith('.bed3') or feature_file.endswith('.gff') or feature_file.endswith('.gff3') else ''),
                     ','.join(all_species),
                     ("""\nWindowSplit(\\
                            preferred_size=%d,\\
                            align=adjust,\\
                            keep_small_blocks=no),\\"""%bootstrapped_distance if bootstrapped_distance else ''),
                     ("""\nOutputDistanceMatrices(\\
                            distance=MLDistance,\\
                            file=data.distmat.ph,\\
                            compression=none,\\
                            strip_names=yes),\\""" if bootstrapped_distance == 0 else ''))) #keep_small_blocks=yes
        subprocess.call('./maffilter param=maf_filter_config.bpp',shell=True)
        if bootstrapped_distance:
            trees_gen = Phylo.parse('./maf_trees/output_tree_consensus.nh', 'newick')
            trees = []
            with open('./maf_trees/output_tree_consensus.nwk','r') as f:
                for i in range(len(f.readlines())):
                    try:
                        trees.append(trees_gen.next())
                    except:
                        pass
            if consensus_algorithm == 'strict':
                tree = CS.strict_consensus(trees)
            elif consensus_algorithm == 'majority':
                tree = CS.majority_consensus(trees,float(major_cutoff))
            else:
                tree = CS.adam_consensus(trees)
            Phylo.write(tree,'./maf_trees/output_tree_consensus_bootstrap.nh','newick')

def maf_change_coordinates(segment,ref_species):
    aln_lines = segment.splitlines()
    for i,line in enumerate(aln_lines):
        if line.startswith('s'):
            lineList = line.split()
            orientation = lineList[4]
            lineList2 = lineList[1].split('.')
            lineList3 = lineList2[-1].split('_')[-2:]
            lineList2[2] = lineList2[2].replace('_'+'_'.join(lineList3),'')
            if orientation == '-':
                lineList[2] = str(int(lineList3[-1])-int(lineList[2]))#-int(lineList[3]))
            else:
                lineList[2] = str(int(lineList3[-2]) + int(lineList[2]))
            lineList[1] = '.'.join(lineList2)#FIXME  '.'.join(lineList2[::2])
            aln_lines[i] = '\t'.join(lineList)
            if lineList2[0] == ref_species:
                chrom = lineList2[2]
                position = int(lineList[2])
    return chrom,position,'\n'.join(sorted(filter(None,aln_lines)))+'\n\n'

@begin.subcommand
def selective_pressure_statistics(cns_config,reference_species, min_block_length, dist_max, window_size, root_species): # FIXME add ingroup outgroup options, maybe add ability to output unrooted tree !!!!!!  12/18/17
    min_block_length = int(min_block_length)
    window_size = int(window_size)
    dist_max = int(dist_max)
    try:
        os.mkdir('./maf_trees')
    except:
        pass
    if cns_config.endswith('.txt'):
        with open(cns_config,'r') as f:
            master_species = parseConfigFindList("masterListSpecies",f)
    else:
        master_species = cns_config.split(',')
    master_species = set([species.split('_')[0] for species in master_species])
    print master_species
    mafFiles = [file for file in os.listdir('.') if file.startswith('FastaOut') and file.endswith('.maf')]
    #FIXME tests
    subprocess.call('rm merged.maf',shell=True)
    for file in mafFiles:
        subprocess.call("sed -e '/Anc/d;/#/d' %s >> merged.maf"%file,shell=True)
    subprocess.call("cat -s merged.maf > temp.maf && mv temp.maf merged.maf",shell=True)#sed 's/\\n\\n\\n/\\n\\n/g'
    #with open('merged.maf','r') as f: #FIXME memory error
    #    txt = f.read().replace('\n\n\n','\n\n')
    #with open('merged.maf','w') as f:
    #    f.write(txt)
    #del txt
    with open('maf_filter_config.bpp','w') as f:
                f.write("""
        input.file=./merged.maf
        input.format=Maf
        output.log=out.log
        maf.filter=\\
            Subset(\\
                    strict=yes,\\
                    keep=no,\\
                    species=(%s),\\
                    remove_duplicates=yes),\\
            MaskFilter(species=(%s)),\\
            MinBlockLength(min_length=%d),\\
            Output(\\
                    file=merged.filtered.maf,\\
                    compression=none,\\
                    mask=yes)
                """%(','.join(master_species),','.join(master_species),min_block_length))
    subprocess.call('./maffilter param=maf_filter_config.bpp',shell=True)
    # FIXME sort them by coordinates then output, and keep order of species
    maf_idx = index_maf_2('merged.filtered.maf')
    count = 0
    maf_sort_structure = []
    with open('merged.filtered.maf','r') as f, open('merged.filtered.new_coords.maf','w') as f2:
        for idx in maf_idx:
            f.seek(idx)
            chrom, position, segment = maf_change_coordinates(f.read(maf_idx[idx] - idx),reference_species)
            maf_sort_structure.append((chrom, position, count))
            count += 1
            f2.write(segment)

    maf_sort_structure = pd.DataFrame(maf_sort_structure).sort_values([0,1])
    maf_idx = index_maf_2('merged.filtered.new_coords.maf')
    maf_idx_sorted = np.array(maf_idx.keys())[maf_sort_structure.as_matrix([2])]
    with open('merged.filtered.new_coords.maf','r') as f, open('merged.filtered.new_coords.sorted.maf','w') as f2:
        for idx in maf_idx_sorted:
            f.seek(idx)
            f2.write(f.read(maf_idx[idx] - idx))
    with open('maf_filter_config.bpp','w') as f:
        f.write("""
        input.file=./merged.filtered.new_coords.sorted.maf
        input.format=Maf
        output.log=out.log
        maf.filter=\\
            Merge(\\
                    species=(%s),\\
                    dist_max=%d),\\
            WindowSplit(\\
                    preferred_size=%d,\\
                    align=adjust,\\
                    keep_small_blocks=no),\\
            RemoveEmptySequences(),\\
            Output(\\
                    file=merged.filtered.new_coords.syntenic.maf,\\
                    compression=none,\\
                    mask=no)
                """%(reference_species,dist_max,window_size))
    subprocess.call('./maffilter param=maf_filter_config.bpp',shell=True)
    maf_idx = index_maf_2('merged.filtered.new_coords.syntenic.maf')
    with open('merged.filtered.new_coords.syntenic.maf','r') as f, open('merged.filtered.new_coords.syntenic.fixed.maf','w') as f2:
        for idx in maf_idx:
            f.seek(idx)
            segment = f.read(maf_idx[idx] - idx)
            if segment.strip('a\n'):
                f2.write(segment+'\n\n')
    maf_configs = []
    maf_configs.append("""
    input.file=./merged.filtered.new_coords.syntenic.fixed.maf
    input.format=Maf
    output.log=out.log
    maf.filter=\\
        XFullGap(species=(%s)),\\
        RemoveEmptySequences(),\\
        Subset(\\
                strict=yes,\\
                keep=no,\\
                species=(%s),\\
                remove_duplicates=yes),\\
        MinBlockLength(min_length=%d), \\
        OutputCoordinates(\\
                file=coordinates.txt,\\
                compression=none,\\
                species=(%s),\\
                output_src_size=yes),\\
        DistanceEstimation(\\
                method=ml,\\
                model=GTR,\\
                gap_option=no_gap,\\
                parameter_estimation=initial,\\
                gaps_as_unresolved=no,\\
                unresolved_as_gaps=yes,\\
                extended_names=yes),\\
        DistanceBasedPhylogeny(\\
                    method=bionj,\\
                    dist_mat=MLDistance),\\
        NewOutgroup(\\
                    tree_input=BioNJ,\\
                    tree_output=BioNJ,\\
                    outgroup=%s),\\
        SequenceStatistics(\\
            statistics=(\\
                BlockCounts,\\
                AlnScore,\\
                BlockLength,\\
                CountClusters(\\
                    tree=BioNJ,\\
                    threshold=0.001),\\
                SiteFrequencySpectrum(\\
                    bounds=(-0.5,0.5,1.5,2.5,3.5,4.5),\\
                    ingroup=(%s),\\
                    outgroup=%s),\\
                SiteStatistics(\\
                    species=(%s)),\\
                DiversityStatistics(ingroup=(%s)),\\
                    ),\\
            ref_species=%s,\\
            file=data.statistics.csv),\\
        OutputTrees(\\
                    tree=BioNJ,\\
                    file=./maf_trees/output_trees.nwk,\\
                    compression=none,\\
                    strip_names=yes)"""%(','.join(master_species),','.join(master_species),int(window_size/2),reference_species,root_species,','.join(master_species),root_species,','.join(master_species),','.join(master_species),reference_species))
    # FIXME master species on 4th entry
    for maf_config in maf_configs:
        with open('maf_filter_config.bpp','w') as f:
            f.write(maf_config)
        subprocess.call('./maffilter param=maf_filter_config.bpp',shell=True)
    # output will be trees and statistics, use pandas and grab tree lengths to pandas, encode array, then add polymorphism density from vcf possibly, and then pca the results and cluster the regions
    # neural network model???
    df = pd.read_csv('data.statistics.csv')
    df['tree_lengths'] = find_tree_lengths('./maf_trees/output_trees.nwk')

#############################################
################### START ###################

@begin.start
def main():
    pass

####################################################
################### DEPRECATED #####################

"""
        if c.src.startswith(ref_species):

    rec_info = []
    with open(maf_file) as f:
        reader = maf.Reader(f)
        while True:
            pos = reader.file.tell()
            rec = reader.next()
            if rec is None:
                break
            rec_info.append((rec,pos))
    rec_info.sort()
    """


"""
sed '/Anc/d' FastaOut1.maf > out_maf.maf

input.file=./out_maf.maf
input.format=Maf
output.log=out.log
maf.filter=\
        Subset(\
                strict=yes,\
                keep=no,\
                species=(Sbicolor,Osativa,Bdistachyon,Bstacei,ABRD,BdtwD,BhybsixD,BhybzeroD,ABRS,BdtwS,BhybsixS,BhybzeroS,Phallii,PvirgatumK,PvirgatumN,ZmaysAGPv3),\
                remove_duplicates=yes),\
        Subset(\
                strict=yes,\
                keep=yes,\
                species=(Bdistachyon),\
                remove_duplicates=yes),\
        VcfOutput(\
                file=vcfs/all.vcf,\

                genotypes=(Sbicolor,Osativa,Bstacei,ABRD,BdtwD,BhybsixD,BhybzeroD,ABRS,BdtwS,BhybsixS,BhybzeroS,Phallii,PvirgatumK,PvirgatumN,ZmaysAGPv3),\
                all=no,\
                reference=Bdistachyon)

        MaskFilter(\
                species=(%s),\
                window.size=10,\
                window.step=1,\
                max.masked=3),\
        MinBlockLength(min_length=%d), \
            XFullGap(species=(%s)),\
            MinBlockLength(min_length=25),\

                CAN ADD TO ABOVE with .gz extension included
                compression=gzip,\
nohup ./maffilter param=filterTest.bpp &
"""
"""
maf_sort_structure = []
with open('merged.filtered.maf','r') as f: # FIXME consider indexing this file if too large and sorting that way
    for segment in f.read().split('\n\n'): # FIXME can turn this to generator in future, especially with large maf size
        if segment:
            maf_sort_structure.append(maf_change_coordinates(segment,reference_species))
maf_sort_structure = pd.DataFrame(maf_sort_structure).sort_values([0,1])
with open('merged.filtered.new_coords.maf','w') as f2:
    for seg in maf_sort_structure.itertuples():
        f2.write(seg[3])
del maf_sort_structure #FIXME may be bad when maf file is hundreds of gb large

segments = f.read().split('\n\n')
    with open('merged.filtered.new_coords.syntenic.maf','w') as f:
        for segment in segments: # FIXME can turn this to generator in future, especially with large maf size
            if segment.strip('a\n'):
                f.write(segment+'\n\n')
"""