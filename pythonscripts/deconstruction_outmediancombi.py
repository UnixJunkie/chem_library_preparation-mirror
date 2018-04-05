#!/bin/python

import sys
import os
import rdkit.Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import rdMolDescriptors
from networkx import *
import colorsys
from operator import itemgetter

import numpy
## Few functions to have less text after
def median(lst):
    return numpy.percentile(a=numpy.array(lst), q=50, interpolation="lower")
def std(lst):
    return numpy.std(numpy.array(lst))
def mean(lst):
    return numpy.mean(numpy.array(lst))
def minimum(lst):
    return numpy.amin(numpy.array(lst))

if len(sys.argv)!=5:
    print "Usage = python script.py experimentalvalues.csv correctivevalues.csv file.mol2 out.mol2"

if not os.path.exists("output_mol2"):
    os.makedirs("output_mol2")

## Read the file with experimental values
smartsfile = open(sys.argv[1], 'r')
smartscontent = smartsfile.readlines()
smartsfile.close()

smarts_energy_value = {}
fragment_objects = []
fragment_name_list = []
frag_counter_dict = {}
counter = 47
outmol=""

for a in range(0,len(smartscontent)):
     smarts = rdkit.Chem.MolFromSmarts(smartscontent[a].split('/')[1])
     counter += 1
     fragment_objects.append(smarts)
     fragment_name_list.append(counter)
     frag_counter_dict[counter] = smartscontent[a].split('/')[0]
     smarts_energy_value[counter] = float(smartscontent[a].split('/')[2])

### REMARK ###
### Sources for the experimental values are ###
### cabani 1981 journal of solution chemistry ###
### wolfenden 1981 biochemistry ###
### marten 1996 j phys chem ###
### FreeSolv "Mobley" database Nov 2015 ###


## Read the file with corrective values
a = open(sys.argv[2], 'r')
smartscorrections = a.readlines()
a.close()

corr_energy_value = {}
corr_objects = []
corr_name_list = []
corr_counter_dict = {}
counter = 0

for a in range(0,len(smartscorrections)):
    smarts = rdkit.Chem.MolFromSmarts(smartscorrections[a].split('/')[1])
    counter += 1
    corr_objects.append(smarts)
    corr_name_list.append(counter)
    corr_counter_dict[counter] = smartscorrections[a].split('/')[0]
    corr_energy_value[counter] = float(smartscorrections[a].split('/')[2])

#print corr_counter_dict
#print corr_name_list
#print corr_energy_value

## Read the mol2 file
a = open(sys.argv[3], 'r')
mol2file = a.readlines()
mol2file.append("@<TRIPOS>MOLECULE\n")
a.close()
mol2string = ''
termination = False
substructure_corr_dict=None

for mol2line in mol2file:
    if len(mol2line.split())==0:
        continue
    elif mol2line[0]=="#":
        continue
    if '@<TRIPOS>MOLECULE' in mol2line:
        """
        Read the mol2 block, renumber the elements to have H at the end
        """
        if len(mol2string) != 0:
            #try:
            mol2 = rdkit.Chem.MolFromMol2Block(mol2string, sanitize=True, removeHs=False)
            #rdkit.Chem.SanitizeMol(mol2,rdkit.Chem.SanitizeFlags.SANITIZE_PROPERTIES|rdkit.Chem.SanitizeFlags.SANITIZE_SYMMRINGS)
            #rdkit.Chem.Kekulize(mol2, clearAromaticFlags=True)
            heavy_atom_list = []
            hydrogen_atom_list = []
            if mol2 != None:
                for atom in mol2.GetAtoms():
                    atomnum = atom.GetAtomicNum()
                    if atomnum > 1 : 
                        heavy_atom_list.append(atom.GetIdx())
                    else:
                        hydrogen_atom_list.append(atom.GetIdx())
                    if atomnum == 15:
                        print "P atom detected, that's not gonna work.", mol2.GetProp("_Name")
                        termination = False
                nbofatoms=len(hydrogen_atom_list)+len(heavy_atom_list)
                #if termination == True:
                #        mol2string = ''
                #        mol2string += mol2line
                #        mol2 = None
                #        continue
            #except:
            #    mol2string = ''
            #    mol2string += mol2line
            #    mol2 = None
            #    print "bug"
            #    continue
            #if rdkit.Chem.GetFormalCharge(mol2) != 0:
            #    print "molecule is charged, skipped.", mol2.GetProp("_Name")
            if mol2 != None:# and rdkit.Chem.GetFormalCharge(mol2) == 0:
                    """try:"""
                    #both_lists = heavy_atom_list + hydrogen_atom_list
                    #mol2_clean_index = None
                    #mol2_clean_index = rdkit.Chem.RenumberAtoms(mol2, both_lists)
                    molname = ""
                    molname = mol2.GetProp("_Name")
                    full = 0
                    """
                    Find fusion carbons
                    """
                    fusionC = []
                    itis=0
                    ringbond=rdkit.Chem.MolFromSmarts('[$(*)]~&@[$(*)]')
                    ringb = mol2.GetSubstructMatches(ringbond)
                    ringbonds=[]
                    for i in ringb:
                        ringbonds.append(list(i))
                    for atom in mol2.GetAtoms():
                        if atom.IsInRing():
                            if atom.GetAtomicNum()==6:
                                if atom.GetDegree() >= 3:
                                    check_bondtype = 0
                                    for neighbor in atom.GetNeighbors():
                                        if sorted([atom.GetIdx(), neighbor.GetIdx()]) in ringbonds:
                                            check_bondtype+=1
                                            if check_bondtype >= 3:
                                                fusionC.append(atom.GetIdx())
                    """
                    Substructure search in the experimental values
                    """
                    substructure_results_dict = {}
                    fragment = None
                    substructure = None
                    fragment_smarts_combined = ""
                    for f in range(0,len(fragment_objects)):
                        fragment = fragment_objects[f]
                        substructure = mol2.GetSubstructMatches(fragment)
                        for k in range(0, len(substructure)):
                            fragment_smarts_combined = str(fragment_name_list[f]) + '_' + str(k)
                            substructure_results_dict[fragment_smarts_combined] = set(substructure[k])

                    """
                    substructure_results_dict contains:
                        Keys = fragment indice number (as in the csv file) + an indice corresponding to a set. A molecule can contain several times
                    a substructure. This indice relates to the "number" of the substructure in the molecule (if it's occurence 0, 1, 2...)
                    The corresponding set is in the value of this key.
                        Value = Set of atom indices of the molecule containing the substructure.
                    
                    Checking which atoms of the molecules are included in some substructures of the experimental data set.
                    Finding missing ones. Used at the end.
                    """
                    present=[]
                    for jr in range(0, len(substructure_results_dict.keys())):
                        present.extend(list(substructure_results_dict.values()[jr]))

                    dedupl = ()
                    dedupl=set(present)
                    """
                    find substructures with fusion carbons
                    """
                    list_withfu=[]
                    for i in fusionC:
                        for n in range (0, len(substructure_results_dict.keys())):
                            if i in substructure_results_dict.values()[n]:
                                list_withfu.append(substructure_results_dict.keys()[n])
                    set_withfu=set(list_withfu)
                    list_uniqfu=list(set_withfu)
                    dict_of_fuC={}
                    """
                    Creation of a dictionary containing as a key a fusion C and as values the substructures it is in.
                    """
                    for i in fusionC:
                        for j in range(0, len(list_uniqfu)):
                            if i in substructure_results_dict[list_uniqfu[j]]:
                                dict_of_fuC.setdefault(i, []).append(list_uniqfu[j])

                    """
                    Find a way to have substructures that have in common only those fusion C
                    """
                    intersec=None
                    fusubs={}
                    """
                    Create the graph with those substructures.
                    
                    First, we find for each substructure which other substructure it has no intersection with.
                    ie, no atom in common.
                    graph_dict contains:
                        keys = a key from substructure_results_dict. Each key will be a node of the graph
                        values = a list of keys from substructure_results_dict with which the key (of graph_dict)
                    has no intersection with, as defined just above in the description.
                    The values will be "non edges" of the graph.
                    """
                    graph_dict = {}
                    for n in range (0, len(substructure_results_dict.keys())):
                        graph_dict[substructure_results_dict.keys()[n]] = []
                        for m in range(0, len(substructure_results_dict.keys())):
                            if len(substructure_results_dict[substructure_results_dict.keys()[n]].intersection(substructure_results_dict[substructure_results_dict.keys()[m]])) == 0:
                                graph_dict[substructure_results_dict.keys()[n]].append(substructure_results_dict.keys()[m])
                    """
                    Add putative redondant structures of fusion C
                    """
                    verify=0
                    chk=0
                    for n in range (0, len(substructure_results_dict.keys())):
                        for m in range(0, len(substructure_results_dict.keys())):
                            if len(substructure_results_dict[substructure_results_dict.keys()[n]].intersection(substructure_results_dict[substructure_results_dict.keys()[m]])) != 0:
                                intersec=substructure_results_dict[substructure_results_dict.keys()[n]].intersection(substructure_results_dict[substructure_results_dict.keys()[m]])
                                verify=len(intersec)
                                chk=0
                                for truc in intersec:
                                    if truc in fusionC:
                                        chk+=1
                                    if verify==chk:
                                        graph_dict[substructure_results_dict.keys()[n]].append(substructure_results_dict.keys()[m])

                    """
                    Now we build the graph and find the cliques.
                    Each node = a substructure
                    Edges are between all non intersected substructures.
                    """
                    graph = None
                    cliques = None
                    graph = networkx.from_dict_of_lists(graph_dict)
                    cliques =  list(find_cliques(graph))
                    """
                    Cliques contains a list of lists of connected nodes (substructures).
                    Object clique will be set to an indice of cliques. 
                    
                    min_clique_length_param is set to the nb of HA, but before being used it's set to 
                    the number of nodes in the cliques (len of clique)
                    
                    then we print some data about how many atoms are covered by the smallest number of maximum cliques.
                    There is a small checker: Only cliques that cover the maximum number of atoms are retained. It improves drastically speed.
                    Then, among those cliques, select those with the fewer nodes (biggest substructures)
                    """
                    min_clique_length = 0
                    min_clique_length = len(heavy_atom_list)
                    all_energies=[]
                    all_cliquelength=[]
                    clique = None
                    coverage_tuple=[]
                    for c in range(0, len(cliques)):
                        coverage = []
                        clique = cliques[c]
                        for t in range(0,len(clique)):
                            coverage = coverage + list(substructure_results_dict[clique[t]])
                        coverage_tuple.append([len(set(coverage))-len(list_uniqfu), len(clique), c])
                    maxcoverage = 0
                    """
                    Introduction of a new flag. If there's 0 experimental value for the molecule, it crashes.
                    This flag is here to overcome this.
                    Several modifications later to solve this
                    """
                    noexpdata=False
                    if len(coverage_tuple) == 0:
                        noexpdata = True
                    if noexpdata == False:
                        maxcoverage = max(coverage_tuple,key=itemgetter(0))[0]
                    #print "The coverage of atoms with experimental values is", maxcoverage, "out of", len(heavy_atom_list)
                    #print "That means we need", len(heavy_atom_list)-maxcoverage, "corrective values."
                    mincliques = 0
                    nbofcliques = []
                    for z in range(0, len(coverage_tuple)):
                        if coverage_tuple[z][0] == maxcoverage:
                            nbofcliques.append(coverage_tuple[z][1])
                    if noexpdata == False:
                        mincliques = min(nbofcliques)
                    #print "The smallest number of maximum cliques is", mincliques
                    interestingcliques = []
                    if noexpdata == False:
                        for u in range(0, len(coverage_tuple)):
                            if coverage_tuple[u][0] == maxcoverage and coverage_tuple[u][1] == mincliques:
                                interestingcliques.append(coverage_tuple[u][2])
                    """
                    checklist is set to 0 for each new clique of cliques.
                    It appends all atoms in each substructure of clique.
                    
                    N is the number of substructures in clique of cliques.
                    """
                    winner=None
                    corr_cliques=[]
                    if noexpdata == False:
                        for e in interestingcliques:
                            clique = cliques[e]
                            #print clique
                            checklist = []
                            N = len(clique)
                            HSV_tuples = [(x*1.0/N, 1.0, 1.0) for x in range(N)]
                            RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
                            atom_to_color_dict = {}
                            test_tri_tuple = []
                            for t in range(0,len(clique)):
                                clique_member = clique[t]
                                test_tri_tuple.append([frag_counter_dict[int(clique[t].split('_')[0])],smarts_energy_value[int(clique[t].split('_')[0])],RGB_tuples[t]])
                                for testindex in list(substructure_results_dict[clique_member]):
                                    atom_to_color_dict[testindex] = RGB_tuples[t]
                                    checklist.append(testindex)
                                if set(checklist) == set(heavy_atom_list):
                                    mol2filename = molname +'_'+str(e)  +'.png' 
                                    full = 1
                                    energy = 0.00
                                    for value in clique:
                                        energy += smarts_energy_value[int(value.split('_')[0])]
                                    min_clique_length_param = len(clique)
                                    all_energies.append(energy)
                                    all_cliquelength.append(N)
                                    """
                                    If all atoms are covered with experimental values, we print the energy.
                                    Otherwise, we try to find corrective values for remaining atoms and directly add their values to the energy total.
                                    Once a corrective value has been found for a missing atom/group, it deletes this atom/group from the missing indices list.
                                    Watch out, I've sorted the corrective smarts csv file in order to have groups first and single atomic contributions at the end.
                                    That could make a difference.
                                    """
                            if full == 0:
                                energy = 0.00
                                clique_withcorr=clique
                                #print clique_withcorr
                                for value in clique:
                                    energy += smarts_energy_value[int(value.split('_')[0])]
                                stillmissing = list(set(heavy_atom_list) - set(checklist))
                                substructure_corr_dict={}
                                for da in range(0, len(corr_objects)):
                                    fragment = corr_objects[da]
                                    substructure = mol2.GetSubstructMatches(fragment)
                                    addendum=0
                                    for k in range(0, len(substructure)):
                                        if len(set(substructure[k]).intersection(stillmissing)) == len(substructure[k]):
                                            energy += corr_energy_value[da+1] 
                                            #print da+1
                                            for tmp in clique_withcorr:
                                                if str(da+1)+"_"+str(addendum) in tmp:
                                                    addendum+=1
                                            clique_withcorr.append(str(da+1)+"_"+str(addendum))
                                            substructure_corr_dict[str(da+1)+"_"+str(addendum)]=set(substructure[k])
                                            #print rdkit.Chem.MolToSmarts(fragment)
                                            #print clique_withcorr
                                            #print substructure_corr_dict
                                            checklist = checklist + list(set(substructure[k]).intersection(stillmissing))
                                            stillmissing = list(set(heavy_atom_list) - set(checklist))
                                if set(checklist) == set(heavy_atom_list):
                                    mol2filename = molname +'_'+str(e)  +'.png' 
                                    #tmp=rdkit.Chem.Draw.MolToImage(mol2, size=(500, 500), fitImage=True,highlightMap=atom_to_color_dict, wedgeBonds=False, legend=str(energy),legend_colorscheme=test_tri_tuple)
                                    #tmp.save(mol2filename)
                                    all_energies.append(energy)
                                    all_cliquelength.append(len(clique_withcorr))
                                    full = 1
                                else:
                                    None
                                    #print "Something is strange. It should be because it didn't find a corrective SMARTS."
                                corr_cliques.append(clique_withcorr)
                                #print corr_cliques
                    else:
                        energy = 0
                        checklist = []
                        stillmissing = heavy_atom_list
                        substructure_corr_dict={}
                        for da in range(0, len(corr_objects)):
                            fragment = corr_objects[da]
                            substructure = mol2.GetSubstructMatches(fragment)
                            addendum=0
                            clique_withcorr=[]
                            for k in range(0, len(substructure)):
                                if len(set(substructure[k]).intersection(stillmissing)) == len(substructure[k]):
                                    energy += corr_energy_value[da+1] 
                                    for tmp in clique_withcorr:
                                        if str(da+1)+"_"+str(addendum) in tmp:
                                            addendum+=1
                                    clique_withcorr.append(str(da+1)+"_"+str(addendum))
                                    substructure_corr_dict[str(da+1)+"_"+str(addendum)]=set(substructure[k])
                                    checklist = checklist + list(set(substructure[k]).intersection(stillmissing))
                                    stillmissing = list(set(heavy_atom_list) - set(checklist))
                        if set(checklist) == set(heavy_atom_list):
                            #mol2filename = molname +'_'+str(e)  +'.png' 
                            #tmp=rdkit.Chem.Draw.MolToImage(mol2, size=(500, 500), fitImage=True,highlightMap=atom_to_color_dict,
                            #tmp.save(mol2filename)
                            all_energies = energy
                            full = 1
                        else:
                            print "Something is strange. It should be because it didn't find a corrective SMARTS."
                    if full == 1 and noexpdata == False:
                        dico={}
                        for ba in range(0, len(all_cliquelength)):
                            dico[all_energies[ba]] = [all_cliquelength[ba], ba]
                        mini=minimum(all_cliquelength)
                        interest_ener=[]
                        for ener, clisize in dico.items():
                            if clisize[0] == mini:
                                interest_ener.append(ener)
                        #print molname + ' is parameterizable!'
                        #print "The median value is", median(interest_ener)
                        for i in dico.items():
                            if i[0]==median(interest_ener):
                                winner=i[1][1]
                        #print winner, interestingcliques[winner], interestingcliques, cliques[220]
                        if len(corr_cliques) != 0:
                            winningclique=corr_cliques[winner]
                        else:
                            winningclique=cliques[interestingcliques[winner]]
                        winningdico={}
                        if substructure_corr_dict != None :
                            substruct_combined=dict(substructure_results_dict, **substructure_corr_dict)
                        else:
                            substruct_combined=substructure_results_dict
                        #print substruct_combined
                        if len(corr_cliques) != 0:
                            for cliqitem in range(0, len(winningclique)):
                                winningdico[corr_cliques[winner][cliqitem]]=substruct_combined[corr_cliques[winner][cliqitem]]
                        else:
                            for cliqitem in range(0, len(winningclique)):
                                winningdico[winningclique[cliqitem]]=substruct_combined[winningclique[cliqitem]]
                        #print winningdico
                        yetanotherdico={}
                        for i in winningdico.values():
                            for j in range(0, nbofatoms):
                                    if j in i:
                                        for keys in winningdico.keys():
                                            if winningdico[keys]==i:
                                                yetanotherdico[j]=keys.split("_")[0]
                        for i in range(0, nbofatoms):
                            if i in yetanotherdico:
                                continue
                            else:
                                yetanotherdico[i]=0
                        #print winningdico
                        tmpdico={}
                        compte=0
                        for keys in winningdico.keys():
                            for item in winningdico[keys]:
                                tmpdico[item]=compte
                            compte+=1
                        for i in range(0, nbofatoms):
                            if i in tmpdico:
                                continue
                            else:
                                tmpdico[i]=compte
                                compte+=1 
                        #dicoofsubsnb={}
                        #print tmpdico
                        #print winningdico
                        #for i in range(0, nbofatoms):
                        #    for keys in winningdico.keys():
                        #        if i in winningdico[keys]
                        print molname.replace(" ",""), median(all_energies)#, realvalue
                        #print "ici tocard", molname, median(interest_ener)
                        writeflag=False
                        #outmol=""
                        atomcounter=0
                        for yetanotherline in mol2string.split("\n"):
                            if len(yetanotherline.split())==0:
                                continue
                            elif yetanotherline[0]=="#":
                                continue
                            if "BOND" in yetanotherline:
                                writeflag=False
                            if writeflag == True:
                                outmol+=('%7s %-8s%10s%10s%10s %-7s%5s %-11s%8s\n') %(yetanotherline.split()[0], yetanotherline.split()[1], yetanotherline.split()[2],
                                    yetanotherline.split()[3], yetanotherline.split()[4], yetanotherline.split()[5], yetanotherdico[atomcounter],
                                    tmpdico[atomcounter], yetanotherline.split()[8])
                                atomcounter+=1
                            if "ATOM" in yetanotherline:
                                writeflag=True
                                outmol+=yetanotherline+"\n"
                            if writeflag == False:
                                outmol+=yetanotherline+"\n"
                    elif full == 1 and noexpdata == True:
                        print molname+" is parameterizable, but only with corrective values:", all_energies
                        yetanotherdico={}
                        for i in substructure_corr_dict.values():
                            for j in range(0, nbofatoms):
                                    if j in i:
                                        for keys in substructure_corr_dict.keys():
                                            if substructure_corr_dict[keys]==i:
                                                yetanotherdico[j]=keys.split("_")[0]
                        for i in range(0, nbofatoms):
                            if i in yetanotherdico:
                                continue
                            else:
                                yetanotherdico[i]=0
                        tmpdico={}
                        compte=0
                        for keys in substructure_corr_dict.keys():
                            for item in substructure_corr_dict[keys]:
                                tmpdico[item]=compte
                            compte+=1
                        for i in range(0, nbofatoms):
                            if i in tmpdico:
                                continue
                            else:
                                tmpdico[i]=compte
                                compte+=1  

                        writeflag=False
                        #outmol=""
                        atomcounter=0
                        for yetanotherline in mol2string.split("\n"):
                            if len(yetanotherline.split())==0:
                                continue
                            elif yetanotherline[0]=="#":
                                continue
                            if "BOND" in yetanotherline:
                                writeflag=False
                            if writeflag == True:
                                outmol+=('%7s %-8s%10s%10s%10s %-7s%5s %-11s%8s\n') %(yetanotherline.split()[0], yetanotherline.split()[1], yetanotherline.split()[2],
                                    yetanotherline.split()[3], yetanotherline.split()[4], yetanotherline.split()[5], yetanotherdico[atomcounter],
                                    tmpdico[atomcounter], yetanotherline.split()[8])
                                atomcounter+=1
                            if "ATOM" in yetanotherline:
                                writeflag=True
                                outmol+=yetanotherline+"\n"
                            if writeflag == False:
                                outmol+=yetanotherline+"\n"

                    else:
                        print molname+" bugged because of reasons. Work harder"
                        print stillmissing, "atom indices are still missing."
                        a=0
                        for atom in mol2.GetAtoms():
                            if a in stillmissing:
                                print "Atomic number", atom.GetAtomicNum()
                            a+=1                            
                    ## Cleaning 
                    mol2string = ''
                    mol2string += mol2line
                    """except:
                    print "it bugged.", molname"""
            mol2string = ''
            mol2string += mol2line
        else:
            mol2string = ''
            mol2string += mol2line
    else:
        mol2string += mol2line

a=open("output_mol2/"+sys.argv[4], "w")
a.write(outmol)
a.close()
