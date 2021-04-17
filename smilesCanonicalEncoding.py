from rdkit import Chem
import copy
from collections import defaultdict

def initializeRingData(m):
    branchEdgesByStart = defaultdict(list)
    branchEdgesByEnd = defaultdict(list)

    # store the mappings from index to atomic number
    mapIndxToAtomicNum = defaultdict(list)

    ringEdges = set()

    for bond in m.GetBonds():
        # This encoding does not support aromatic rings so request a new molecule.
        if str(bond.GetBondType()) == "AROMATIC":
            x = input("This molecule is aromatic. Please try again.\n")
            return encodeSMILES(x)

        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        startAtomNum, endAtomNum = bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()
        
        # add this item to the map:
        mapIndxToAtomicNum[start] = startAtomNum
        mapIndxToAtomicNum[end] = endAtomNum
        if bond.IsInRing():
            ringEdges.add((start,end))           
        else:
            branchEdgesByStart[start].append(end)
            branchEdgesByEnd[end].append(start)

    return branchEdgesByStart, branchEdgesByEnd, ringEdges, mapIndxToAtomicNum


# CREATING INITIAL RING REPRESENTATIONS:
# O: ordering of rings and the traversal + starting point for each ring + map from fragment id to list of rings it belongs to
# for each ring: <------- gives the best representation for each ring 
#   for each start:
#    for each direction:
#       ":" 
#        for each point on the ring traversal, add the point, for each fragment rooted at that point (dont forget ot remove th point itself), create all the fragments and attach them in lexicographic order, each in separate parenths
#        (map the fragment to the list of rings it is attached to -- branch to list of rings mapping)
#        (keep track of the traversal order)
#        ":"
# SORT RING REPRESENTATIONS (LEXICOGRAPHIC) -> returns ordering of rings (branches duplicated)
# PERFORM BRANCH TRIMMING AND MERGING INTO ONE MOLECULE:
# for each ring (in order):
#     do the established traversal, skipping the used branches (keep track of the fragments that have been used already)
#     + "-" when you finish adding the ring


def recursiveJoin():
    pass

# creates individual string molecules from dictionary for each branch submolecule
# creates a map to map each molecule back to the indices of the atoms in it
# TODO: not incldue the atom if it is on the ring (check if an atom is attached to something else in the entire molecule)
# TODO: MolFragmentToSmiles
# 
def createMoleculesFromBranches(allTrees, mapIndxToAtomicNum):
    allBranchMols = []
    for branchDict in allTrees.values():
        # create empty editable mol object
        mol = Chem.RWMol()
        node_to_idx = {}
        for anAtomIdx in branchDict.keys():
            a = Chem.Atom(mapIndxToAtomicNum[anAtomIdx])
            molIdx = mol.AddAtom(a)
            node_to_idx[anAtomIdx] = molIdx

        # add bonds between adjacent atoms
        for anAtom in branchDict.keys():
            for neibAtom in branchDict[anAtom]:
                # TODO: doing only single bonds -- hopefully not a problem
                bond_type = Chem.rdchem.BondType.SINGLE
                if (anAtom < neibAtom):
                    mol.AddBond(node_to_idx[anAtom], node_to_idx[neibAtom], bond_type)
                else:
                    mol.AddBond(node_to_idx[neibAtom], node_to_idx[anAtom], bond_type)
        # Convert RWMol to Mol object
        mol = mol.GetMol()
        allBranchMols.append((Chem.MolToSmiles(mol), list(branchDict.keys())))
    print("all mols: ", allBranchMols)
    return allBranchMols

# FIX: Need to backtrack for branches that diverge
# branchEdgesByStart/End values are lists

# TODO: change the branches, such that the string is unique and also keep track of all the tree's endpoints
# make sure it gives me smiles canonical strings
# TODO: figure out which branch the bond belongs to
# list of atoms in this branch; when you cant add anything else, then you know the branch is done
# keep the list of atoms and the bonds and you pass the list of bonds to smiles and it does it for you
# look at winsotn's implementation and the visitor thing to figure out whether I can use it here
def getJoinedBranches(branchEdgesByStart, branchEdgesByEnd, mapIndxToAtomicNum):
    print("branchEdgesByStart:", branchEdgesByStart)
    # while there are still non-ring edges we have not used
    allTrees = dict()
    ctr = 0
    while len(branchEdgesByStart) > 0:
        thisTreeNodeDict = dict() # a set to keep track of the nodes in the current tree
        
        # add the starting item for this tree:
        (s,neighbors) = branchEdgesByStart.popitem()
        thisTreeNodeDict[s] = neighbors
        queue = set(copy.copy(neighbors)) # keep the elements in the queue that we wish to explore
        # keep adding all the edges that belong to this tree
        while (len(queue) > 0):
            elem = queue.pop()
            neibs = branchEdgesByStart[elem]
            thisTreeNodeDict[elem] = neibs
            # remove these edges from branchEdgesByStart
            del branchEdgesByStart[elem]
            for i in neibs:
                if i not in thisTreeNodeDict:
                    queue.add(i)
        allTrees[ctr] = thisTreeNodeDict
        ctr += 1
    # get all the molecular representations of the branches; allBranchMols is a list of tuples (u,v)
    # where u is the string representation of the branch and v is the list of atom indecies involved in that branch
    allBranchMols = createMoleculesFromBranches(allTrees, mapIndxToAtomicNum)
    return allBranchMols

def getRingData(m, ringVectors):
    branchEdgesByStart, branchEdgesByEnd, ringEdges, mapIndxToAtomicNum = initializeRingData(m)
    joinedBranches = getJoinedBranches(branchEdgesByStart, branchEdgesByEnd, mapIndxToAtomicNum)

    ringEdgesByStart = defaultdict(dict)
    ringEdgesByEnd = defaultdict(dict)
    while len(ringEdges) != 0:
        (s,e) = ringEdges.pop()
        for ring in ringVectors:
            if s in ring:
                ringEdgesByStart[ring][s] = e
                ringEdgesByEnd[ring][e] = s
                    
    return joinedBranches, ringEdgesByStart, ringEdgesByEnd

def getRingTraversals(joinedBranches, ringEdgesByStart, ringEdgesByEnd):
    ringVectors = ringEdgesByStart.keys()
    ringBranches = defaultdict(list) # Values are all branches ordered in the same way as the values for ringLongestBranch below
    #FIX: NEED TO CHANGE ringLongestBranch FOR TIES 
    ringLongestBranch = defaultdict(list) # Values are longest branch ordered so the left most element is NOT in 
                                            # contact with the ring and the right most element IS in contact with ring
    # TODO: implement the queue structure to ensure that you try all the branches if needed 
    for branch in joinedBranches: #[1,2,3] "CCO", (1,3)
        for ring in ringVectors:
            # check if the branch is attached to the ring
            if (branch[0] in ring):
                branchRev = branch[::-1]
                ringBranches[ring].append(branchRev) # standard order
                if len(branchRev) > len(ringLongestBranch[ring]):
                    ringLongestBranch[ring] = branchRev
            elif (branch[-1] in ring):
                ringBranches[ring].append(branch)
                if len(branch) > len(ringLongestBranch[ring]):
                    ringLongestBranch[ring] = branch
        
    minTraversals = defaultdict(list)
    # try different directions of traversing the ring
    for ring in ringVectors:
        traversal1 = copy.copy(ringLongestBranch[ring])
        traversal2 = copy.copy(ringLongestBranch[ring])
        # TODO: incorporate all branches stored ringBranches, make sure that you do not add the longest branch more than once
        while (len(traversal1) < (len(ring)+1)):
            traversal1.append(ringEdgesByStart[ring][traversal1[-1]])
            traversal2.append(ringEdgesByEnd[ring][traversal2[-1]])
        
        # TODO: do the encoding and then do the min traversal
        minTraversals[ring] = min(traversal1,traversal2) # FIX: After taking min, need to embed remaining branches from ringBranches when converting to string
    return minTraversals

# Converts SMILES to canonical SMILES
def encodeSMILES(s):
    m = Chem.MolFromSmiles(s) # gets the molecule
    # Asks for input again if molecule not valid SMILES
    if m == None:
        x = input("Please try again.\n")
        return encodeSMILES(x)              
    
    ri = m.GetRingInfo()
    n = ri.NumRings()

    # If there are any rings in structure, we need to perform encoding so 
    # molecule can be parsed using context-free grammar.
    if n > 0:
        ringVectors = ri.AtomRings() # atoms on the ring
        joinedBranches, ringEdgesByStart, ringEdgesByEnd = getRingData(m, ringVectors) # joinedBranches -- figured out 
        # which atoms were not in the ring and keep extending until you hit a ring atom; used up all the edges
        # ringEdgesByStart/End -- gives you the edges and their directions as dictionaries
        minRingTraversals = getRingTraversals(joinedBranches, ringEdgesByStart, ringEdgesByEnd)
    else:
        pass
        # Output SMILES

# x = input("Please enter a valid SMILES encoding of a molecule.")
# encodeSMILES(x)
cubane = 'C12C3C4C1C5C4C3C25'
bicyclohexyl = 'C1CCCCC1C1CCCCC1'
cyclohexane = 'C1CCCCC1'
_3_propyl_4_isopropyl_1_heptene = 'C=CC(CCC)C(C(C)C)CCC'
_3_ringed_multiple_branches = 'C1(C2=CC(F)=C(C=C2N(C2CC2)C=C1C(=O)O)N1CCNCC1)=O'
_1_methyl_3_bromo_cyclohexene_1 = 'CC1=CC(CCC1)Br'
encodeSMILES(_1_methyl_3_bromo_cyclohexene_1)
