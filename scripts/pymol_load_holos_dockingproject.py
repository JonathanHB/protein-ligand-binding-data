import sys
import string
import pickle

upperpath = "/Users/jonathanborowsky/mount"
directory = f"{upperpath}/bowmanlab/borowsky.jonathan/docking-project"

#set these to match a filename in output-directory
rank = 6
pdbid = "2bsm"

#0 carbonic anhydrase
#1 brd4 bromodomain
#2 something with <10 nonpolymeric ligands
#3 carbonic anhydrase
#4 carbonic anhydrase
#5 2 beta sheets, interesting
#6 hsp90, interesting

input_file = open(f"{directory}/output-files/{rank}-{pdbid.upper()}", "rb")
protlist = pickle.load(input_file)["Kd"]
input_file.close()


cmd.delete("all")

#duplicates sometimes emerge if there are multiple valid ligands in a
#single chain and cause the script to crash
pnames = []

#used for alignment and centering
prot0 = protlist[0][0]+protlist[0][2]

#fetch and align proteins
for x, prot in enumerate(protlist):

    print(prot)

    pname = prot[0]+prot[2]
    if pname in pnames:
        print("Skipping duplicate protein chain, it has another ligand which may not currently be displayed.")
        continue
    pnames.append(pname)

    if len(prot[1].split(" ")) > 1:
        print("skipping polymeric ligand")
        continue

    cmd.fetch(pname)

    if x != 0:
        cmd.align(pname, prot0)

    cmd.hide("sticks", f"{pname} and not resn {prot[1]}")


cmd.center(prot0)

#graphics settings and hiding highly soluble ligands
print(f"hide sticks, resn GOL+EDO+FMT+DMS+ACT+NO3 or elem H; hide spheres, resn NA+CL+K+PO4+SO4+NO3; \
hide nonbonded; hide cartoon, not {prot0}; util.cbag; center {prot0}")

#less soluble but still MOAD-invalid substances
print("hide sticks, resn PEG+PGE+BU3")

#more important ions
print("hide spheres, resn CA+NI+IOD")
