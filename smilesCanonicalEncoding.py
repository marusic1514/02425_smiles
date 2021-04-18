from rdkit import Chem
import copy
from collections import defaultdict

def initializeRingData(m):
    branchEdgesByStart = defaultdict(list)
    branchEdgesByEnd = defaultdict(list)

    ringEdges = set()

    for bond in m.GetBonds():
        # This encoding does not support aromatic rings so request a new molecule.
        if str(bond.GetBondType()) == "AROMATIC":
            x = input("This molecule is aromatic. Please try again.\n")
            return encodeSMILES(x)

        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        startAtomNum, endAtomNum = bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()
        
        # add this item to the map:
        if bond.IsInRing():
            ringEdges.add((start,end))           
        else:
            branchEdgesByStart[start].append(end)
            branchEdgesByEnd[end].append(start)

    return branchEdgesByStart, branchEdgesByEnd, ringEdges

# creates individual string molecules from dictionary for each branch submolecule
# creates a map to map each molecule back to the indices of the atoms in it
# TODO: not incldue the atom if it is on the ring (check if an atom is attached to something else in the entire molecule)
# TODO: MolFragmentToSmiles
# 
def moleculesFromBranches(allTrees,molecule,ringAtoms):
    print(ringAtoms, "HERE")
    allTreeMolecules = dict()
    for tree in allTrees:
        atoms, bonds = tree
        for root in atoms:
            a = []
            b = []
            for atom in atoms:
                if atom == root or atom not in ringAtoms:
                    a.append(atom)
            for (x,y) in bonds:
                if x == root or y == root:
                    b.append(bond)
                elif x not in ringAtoms or y not in ringAtoms:
                    b.append(bond)

            m = Chem.MolFragmentToSmiles(molecule, atomsToUse=a, bondsToUse=b, rootedAtAtom=root)
            allTreeMolecules[root] = m[1:]
    #print(allTreeMolecules)
    return allTreeMolecules


# FIX: Need to backtrack for branches that diverge
# branchEdgesByStart/End values are lists

# TODO: change the branches, such that the string is unique and also keep track of all the tree's endpoints
# make sure it gives me smiles canonical strings
# TODO: figure out which branch the bond belongs to
# list of atoms in this branch; when you cant add anything else, then you know the branch is done
# keep the list of atoms and the bonds and you pass the list of bonds to smiles and it does it for you
# look at winsotn's implementation and the visitor thing to figure out whether I can use it here
def getAllTreeBranches(branchEdgesByStart,m):
    # while there are still non-ring edges we have not used
    allTrees = []
    allTreesPrintView = []
    while len(branchEdgesByStart) > 0:
        thisTreeNodeAtoms = set() # a set to keep track of the nodes in the current tree
        thisTreeNodeBonds = set()
        # add the starting item for this tree:
        (s,neighbors) = branchEdgesByStart.popitem()
        thisTreeNodeAtoms.add(s)
        thisTreeNodeAtoms.union(set(neighbors))
        thisTreeNodeBonds.union(set([(s,x) for x in neighbors]))
        queue = set(copy.copy(neighbors)) # keep the elements in the queue that we wish to explore
        # keep adding all the edges that belong to this tree
        while (len(queue) > 0):
            elem = queue.pop()
            neibs = branchEdgesByStart[elem]
            thisTreeNodeAtoms.add(elem)
            thisTreeNodeAtoms.union(set(neibs))
            thisTreeNodeBonds.union(set([(s,x) for x in neibs]))
            # remove these edges from branchEdgesByStart
            del branchEdgesByStart[elem]
            for i in neibs:
                if i not in thisTreeNodeAtoms:
                    queue.add(i)
        allTrees.append((thisTreeNodeAtoms,thisTreeNodeBonds))
    # get all the molecular representations of the branches; allBranchMols is a list of tuples (u,v)
    # where u is the string representation of the branch and v is the list of atom indecies involved in that branch
    return allTrees

def getRingData(m, ringVectors, ringEdges):
    ringEdgesByStart = defaultdict(dict)
    ringEdgesByEnd = defaultdict(dict)
    while len(ringEdges) != 0:
        (s,e) = ringEdges.pop()
        for ring in ringVectors:
            if s in ring:
                ringEdgesByStart[ring][s] = e
                ringEdgesByEnd[ring][e] = s             
    return ringEdgesByStart, ringEdgesByEnd


# finds the ring representation in one direction (direction is set by the ringEdges)
def getTraversalForStartingPoint(m, startingIndex, ring, ringEdges, branchesByRoot):
    order = [startingIndex]
    currIndex = startingIndex
    traversal = ":"
    traversal = traversal + str(m.GetAtomWithIdx(currIndex).GetSymbol())
    if currIndex in branchesByRoot.keys():
        traversal = traversal + "(" + branchesByRoot[currIndex] + ")"
    currIndex = ringEdges[currIndex]
    #print("all trees: ", allTrees)
    # iterate traverse the ring until you circle back around, finding all the branches incident
    # on that point in the ring adding each one in parenths, and then add the atom itself
    #print(traversal, currIndex)
    while (currIndex != startingIndex):
        order.append(currIndex)
        traversal = traversal + str(m.GetAtomWithIdx(currIndex).GetSymbol())
        if currIndex in branchesByRoot.keys():
            traversal = traversal + "(" + branchesByRoot[currIndex] + ")"
        currIndex = ringEdges[currIndex]
    traversal = traversal + (":")
    return traversal, order
    #print("getTraversalForStartingPoint all molecules", allMolecules)
    

# creates the ordering of the rings by attaching all the relevant branches to each ring (allows for branch reusal)
# keeps track of the best traversal (the one that is lexicographically shortest)
# To do so, it tries each possible starting point on the ring and each one of the two rotations (clockwise and coutnerclockwise)
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

# INPUT: ringEdgesByStart -- maps a tuple representing atoms in a ring to the dicitonary of the forward edges
#        ringEdgesByEnd -- maps a tuple representing atoms in a ring to a dictionary of the reversed edges
#        mapIndxToAtomicNum -- maps the index of an atom in a molecule to the atomic number of that atom
#        allTreeBranches -- maps branch index to the edges in it
def getRingTraversals(m, branchesByRoot, ringEdgesByStart, ringEdgesByEnd):
    completeMolecule = []
    stringToRingOrder = dict()
    ringVectors = ringEdgesByStart.keys()
    for ring in ringVectors:
        allTraversals = []
        traversalToOrder = dict()
        for startingIndex in ring:
            # try forward direction:
            ringTraversalForward,atomOrderingF = getTraversalForStartingPoint(m, startingIndex, ring, ringEdgesByStart[ring], branchesByRoot)
            # try reverse direction: 
            ringTraversalReversed,atomOrderingR = getTraversalForStartingPoint(m, startingIndex, ring, ringEdgesByEnd[ring], branchesByRoot)
            
            allTraversals.append(ringTraversalForward)
            allTraversals.append(ringTraversalReversed)
            traversalToOrder[ringTraversalForward] = atomOrderingF
            traversalToOrder[ringTraversalReversed] = atomOrderingR

        allTraversals.sort()
        minTraversal = allTraversals[0]
        stringToRingOrder[minTraversal] = traversalToOrder[minTraversal]
        completeMolecule.append(allTraversals[0])

    completeMolecule.sort()
    moleculeOrderings = [stringToRingOrder[mol] for mol in completeMolecule]
    result = ""
    print(moleculeOrderings)
    print(branchesByRoot)
    for j in range(len(moleculeOrderings)):
        traversal = moleculeOrderings[j]
        result = result + ":"
        for i in range(len(traversal)):
            atomIdx = traversal[i]
            atom = m.GetAtomWithIdx(atomIdx).GetSymbol()
            result = result + atom
            if i == len(traversal)-1:
                result = result + ":"
            if atomIdx in branchesByRoot.keys():
                result = result + "(" + branchesByRoot[atomIdx] + ")"
                del branchesByRoot[atomIdx]
        if j != len(moleculeOrderings)-1:
            result = result + "-"        
    print(result)
    return result

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
        print(ringVectors)
        ringAtoms = set()
        for ring in ringVectors:
            for r in ring:
                ringAtoms.add(r)
        branchEdgesByStart, branchEdgesByEnd, ringEdges = initializeRingData(m)
        branchEdgesByStartCopy = copy.copy(branchEdgesByStart) # make copy because getAllTreeBranches is destructive
        allTreeBranches = getAllTreeBranches(branchEdgesByStartCopy,m)
        ringEdgesByStart, ringEdgesByEnd = getRingData(m, ringVectors, ringEdges)
        branchesByRoot = moleculesFromBranches(allTreeBranches, m, ringAtoms)
        # which atoms were not in the ring and keep extending until you hit a ring atom; used up all the edges
        # ringEdgesByStart/End -- gives you the edges and their directions as dictionaries
        minRingTraversals = getRingTraversals(m, branchesByRoot, ringEdgesByStart, ringEdgesByEnd)
    else:
        return s
        # Output SMILES

if __name__ == "__main__":
    # x = input("Please enter a valid SMILES encoding of a molecule.")
    # encodeSMILES(x)
    cubane = 'C12C3C4C1C5C4C3C25'
    bicyclohexyl = 'C1CCCCC1C1CCCCC1'
    cyclohexane = 'C1CCCCC1'
    _3_propyl_4_isopropyl_1_heptene = 'C=CC(CCC)C(C(C)C)CCC'
    _3_ringed_multiple_branches = 'C1(C2=CC(F)=C(C=C2N(C2CC2)C=C1C(=O)O)N1CCNCC1)=O'
    _1_methyl_3_bromo_cyclohexene_1 = 'CC1=CC(CCC1)Br'
    encodeSMILES(_1_methyl_3_bromo_cyclohexene_1)
