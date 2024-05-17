# Import needed modules
import os
import sys
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo

#create file input and output variables of paths 
input = "/shared/forsythe/BB485/Week06/Brass_CDS_seqs/"
output = "/scratch/rollers/week06/phylo_out/" 

#create a list of all input files 
inputfiles = glob.glob(input+"*fasta")

### Step 1: align sequences using mafft 
for file in inputfiles: 
    new_file_path = file.replace(input, output)
    #print(new_file_path)
    aln_cmd = 'mafft --auto --quiet '+file+' > '+new_file_path
    #print(aln_cmd)
    os.system(aln_cmd)

### Step 2: use tree IQ using system callfor phylogenomic analysis 
#create a list of the aligned file names from above  
alnfiles = glob.glob(output+"*fasta")
#print(alnfiles)

for aln in alnfiles: 
    tree_command = f"iqtree -s {aln} -m TEST -nt 2"
    #print(tree_command)
    os.system(tree_command)

###Step 3: read the trees in and test the topology
#read in the trees by creating a list 
treefiles = glob.glob(output+"*treefile")
#print(treefiles)

#create an empty list 
topology_list=[]

#create a loop to test the topology 
for tree in treefiles: 
    #Read in the tree and store as phylo object
    temp_tree = Phylo.read(tree, "newick")

    #Loop through the tips in the tree to find which one contains Es (the outgroup)
    for tip in temp_tree.get_terminals():
	    if "Es_" in tip.name:
		    es_tip = tip
		    #Stope the loop once we found the correct tip
		    break
    
    #Root the tree by the outgroup taxon
    temp_tree.root_with_outgroup(es_tip)
    
    #Get a list of all terminal (aka tips) branches
    all_terminal_branches = temp_tree.get_terminals()
    
    #Loop through the branches and store the names of the tips of each
    for t in all_terminal_branches:
        if "Bs_" in t.name:
            Bs_temp=t 
        elif "Cr_" in t.name:
            Cr_temp=t
        elif "At_" in t.name:
            At_temp=t
        else:
            out_temp=t
        
    #Make lists of pairs of branches, so that we can ask which is monophyletic
    P1_and_P2=[Bs_temp, Cr_temp]
    P1_and_P3=[Bs_temp, At_temp]
    P2_and_P3=[Cr_temp, At_temp]
    
    #Use series of if/else statemetns to ask which pair in monophyletic
    if bool(temp_tree.is_monophyletic(P1_and_P2)):
        topo_str = "12top"
    elif bool(temp_tree.is_monophyletic(P1_and_P3)):
        topo_str = "13top"
    elif bool(temp_tree.is_monophyletic(P2_and_P3)):
        topo_str = "23top"
    else:
        topo_str = "Unknown"

    topology_list.append(topo_str) 

#save the number of each type of topology as a new variable
total12top = topology_list.count("12top")
total13top = topology_list.count("13top")
total23top = topology_list.count("23top")

#create print statements showing the number of each topology for the different sisters 
with open("output.txt", "a") as f:
    print(f"The number of trees with Bs and Cr sister is {total12top}", file = f)
    print(f"The number of trees with Bs and At sister is {total13top}",file = f)
    print(f"The number of trees with Cr and At sister is {total23top}",file = f)
