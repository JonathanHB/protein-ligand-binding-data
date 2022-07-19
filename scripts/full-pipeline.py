import os
#from operator import itemgetter
import numpy as np
import csv
import pickle

upperdir = "/project/bowmanlab/borowsky.jonathan"
directory = f"{upperdir}/docking-project"

#directory structure
#binding-data-extraction/
#   scripts/
#   moad-fasta-files/
        #every_part_a/
        #every_part_b/
#   usearch-files/
#   output-files/

#-------------------------------------------------------------------------------
#       Part 1: assemble FASTA files from MOAD and cluster by sequence
#-------------------------------------------------------------------------------

#download MOAD FASTA
#get the PDB IDs of all the MOAD structures
#structures = np.hstack([np.unique([i[0:4].upper() for i in os.listdir(f"{directory}/every_part_a/BindingMOAD_2020/")]), np.unique([i[0:4].upper() for i in os.listdir(f"{directory}/every_part_b/BindingMOAD_2020/")])])
#print(f"{len(structures)} structures")

#existing_fasta = [i[0:4] for i in os.listdir(f"{directory}/rcsb_fasta_{part}/")]

#for s in structures:
#    if s not in existing_ids:
#        os.system(f"wget https://www.rcsb.org/fasta/entry/{s}/display -O {directory}/rcsb_fasta_{part}/{s}.fasta")


#note that this requires FASTA files downloaded by step 8. (blastp-moad2pdb.py) of the cryptic pocket filtering pipeline
#https://github.com/JonathanHB/cryptic-pocket-filtering-3
#if you are not using existing bowmore directories, change fasta_directory to
#/path/to/cryptic-pocket-filtering-3/iofiles-blast/rcsb-fasta-ab/
fasta_directory = f"{upperdir}/FAST-cs/protein-sets/new_pockets/iofiles/rcsb_fasta_all"

#path to usearch executable
#install usearch if necessary
usearch_dir = f"{upperdir}/software/"

usearch_io_dir = f"{directory}/usearch-files"
#get fasta files
combined_fasta_name = "moad-fasta-all"

run_usearch = True

if run_usearch:

    #bracketed brackets are an escape char to obtain a single pair of curly brackets
    os.system(f"find {fasta_directory} -type f -exec cat {{}} \; > {usearch_io_dir}/{combined_fasta_name}")
    #^^uncomment^^

    #run usearch clustering
    os.system(f"export PATH=$PATH:{usearch_dir}")

    seqid_threshold = 1
    centroid_name = f"moad-centroids"
    msa_name = f"moad-msa"

    #trying to feed os.system a path for the -centroids argument did not work
    #note that os.system("cd [directory]") does not work
    os.chdir(usearch_io_dir)
    os.system(f"usearch11.0.667_i86linux32 -cluster_fast {combined_fasta_name} -id {seqid_threshold} -centroids {centroid_name} -msaout {msa_name}")
    #^^uncomment^^

#-------------------------------------------------------------------------------
#             Part 1.1: parse USEARCH clusters into lists of PDB IDs

#get a list of the msa files for each cluster
msa_files = [f for f in os.listdir(usearch_io_dir) if f[0:len(msa_name)] == msa_name]
#print(msa_files)

#get a list of the pdb ids of the centroids of each cluster
centroid_ids = []
centroid_pdb_ids = []

for line in open(centroid_name):
    if line[0] == ">":
        segments = line.split("|")
        pdbid = segments[0][1:5]
        chainids = [s.strip() for s in segments[1][7:].split(",")]
        centroid_ids.append([pdbid, chainids])
        centroid_pdb_ids.append(pdbid)
#print(centroid_ids)

clusters = []

#get lists of the proteins in each cluster
print("Processing msa files to collect cluster contents.")

for x, file in enumerate(msa_files):
    if x%500 == 0:
        print(x/len(msa_files))

    #extract the cluster size and proteins in each cluster

    proteins = []
    pdbids = [] #usearch msa files sometimes contain duplicate entries, which I don't want

    for line in open(file):
        if line[0] == ">":
            segments = line.split("|")
            pdbid = segments[0][1:5]
            chainids = [s.strip() for s in segments[1][7:].split(",")]

            if pdbid not in pdbids: #avoid duplicates
                proteins.append([pdbid, chainids])

            pdbids.append(pdbid)

    clusters.append([proteins, int(file[len(msa_name):]), len(proteins)])


#sort clusters in order of descending size and save them
indexer = lambda x : x[2]
clusters_sorted = sorted(clusters, key = indexer, reverse = True)

print(clusters_sorted[0])

np.save(f"{directory}/output-files/moad-clusters-by-size", clusters_sorted)

#-------------------------------------------------------------------------------
#            Part 2: Download binding data for the proteins in each
#                    cluster from MOAD
#-------------------------------------------------------------------------------

moadxmldir = f"{directory}/moad-xml"

#get a list of already-downloaded ligand lists to avoid downloading extra copies
existing_xids_moad = [i[0:4] for i in os.listdir(moadxmldir)]

#-------------------------------------------------------------------------------
#           Part 2.1: define method for extracting binding data from the
#                     csv files available from MOAD

#get biological ligands and their binding constants for a given protein from MOAD
def getligs_binding(struct, moadxmldir, existing_xids_moad): #, chain="", bindingdict):

    #WARNING: The use of --no-check-certificate compromises security
    #ligand information (contains validity information not found in pdb structure)
    if struct.lower() not in existing_xids_moad:
        try:
            os.system(f"wget https://www.bindingmoad.org/files/csv/{struct.lower()}.csv -O {moadxmldir}/{struct.lower()}.csv --no-check-certificate")
            existing_xids_moad.append(struct)
            #keeps the list of existing xml files up to date for when the code encounters apo candidates which are in moad and were previously loaded as holo candidates
        except:
            print(f"{struct} not in moad")
            return []

    bindingdict = {"Ki":[],"Kd":[],"":[],"other":[]}

    #chain_elems = chain.split("_")

    with open(f'{moadxmldir}/{struct.lower()}.csv') as csvfile:
        reader = csv.reader(csvfile)
        firstrow = True
        for row in reader:

            if firstrow: #skip first row, which lacks ligand information
                firstrow = False
                continue

            contents = row[3].split(":")
            if row[4]=="valid":

            #the latter case (chain == "") is for use of this method in get_struct(),
            #where it is used to obtain a list of all chains in a pdb structure containing biologically relevant ligands
                datapoint = [struct] + contents + [row[5]] + row[7:9]

                if row[5] in bindingdict.keys():
                    bindingdict[row[5]].append(datapoint)
                else:
                    bindingdict["other"].append(datapoint)

    return bindingdict

#-------------------------------------------------------------------------------
#            Part 2.2: compile binding data for the largest n clusters
#                      and sort by the number of unique ligands

#the maximum number of clusters to process
n_clusters = 100
#the minimum cluster size to process
min_cluster_size = 50
#Note that there is not a 1:1 mapping between cluster size and number
#of unique ligands with k_d data, but cluster size is the best available
#heuristic.

binding_all_allclusters = []

for x, cluster in enumerate(clusters_sorted[0:n_clusters]):

    if cluster[2] >= min_cluster_size:
        print("processing cluster")

        binding_all = {"Ki":[],"Kd":[],"":[],"other":[]}

        for prot_chain in cluster[0]:
            prot = prot_chain[0]
            #print(prot)
            bindingdict = getligs_binding(prot, moadxmldir, existing_xids_moad)
            #print(bindingdict)

            for k in binding_all.keys():
                 binding_all[k]+=bindingdict[k]

        uniqueligs = np.unique([i[1] for i in binding_all['Kd']])

        binding_all_allclusters.append([binding_all, len(uniqueligs), x, cluster[0][0][0]])

    else:
        print(f"Cluster {x} has only {cluster[2]} proteins; skipping it.")

#sort clusters by the number of unique ligand molecules with Kd data
kd_indexer = lambda x : x[1]

binding_all_allclusters_sorted = sorted(binding_all_allclusters, key = kd_indexer, reverse = True)

for x, cluster_bdata in enumerate(binding_all_allclusters_sorted):

    output_file = open(f"{directory}/output-files/{x}-{cluster_bdata[3]}", "wb")
    pickle.dump(cluster_bdata[0], output_file)
    output_file.close()
