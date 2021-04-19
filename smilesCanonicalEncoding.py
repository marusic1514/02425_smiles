from rdkit import Chem
import copy
from collections import defaultdict

def initializeRingData(m, ringAtoms):
    branchEdgesAll = defaultdict(list)

    ringEdges = set()

    for bond in m.GetBonds():
        # This encoding does not support aromatic rings so request a new molecule.
        # NOTE: I think this rejects aromatic molecules encoded with
        # double/single bonds, so maybe should not do this check
        if str(bond.GetBondType()) == "AROMATIC":
            x = input("This molecule is aromatic. Please try again.\n")
            return encodeSMILES(x)

        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        startAtomNum, endAtomNum = bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()

        # add this item to the map:
        if bond.IsInRing():
            ringEdges.add((start,end))
        else:
            # Exclude mapping branch edges from ring atoms
            # If included, may combine branches with same root into one branch
            if end not in ringAtoms or\
                    (end in ringAtoms and start in ringAtoms):
                branchEdgesAll[end].append(start)
            if start not in ringAtoms or\
                    (end in ringAtoms and start in ringAtoms):
                branchEdgesAll[start].append(end)

    return ringEdges, branchEdgesAll

# creates individual string molecules from dictionary for each branch submolecule
# creates a map to map each molecule back to the indices of the atoms in it

def moleculesFromBranches(allTrees,molecule,ringAtoms):
    allTreeMolecules = defaultdict(list)
    branchAtoms = []
    idxSyms = []
    for i in range(molecule.GetNumAtoms()):
        idxSyms.append(molecule.GetAtomWithIdx(i).GetSymbol() + str(i))
    #print(idxSyms)
    #print(allTrees)
    for tree in allTrees:
        atoms, bonds = tree
        for root in atoms:
            a = list(atoms)
            b = list(bonds)
            branchAtoms.extend(a)
            m = Chem.MolFragmentToSmiles(molecule, atomsToUse=a, bondsToUse=b,
                    rootedAtAtom=root, atomSymbols=idxSyms)
            #print(m)
            mstr = Chem.MolFragmentToSmiles(molecule, atomsToUse=a, bondsToUse=b,
                    rootedAtAtom=root)
            #print(mstr)
            i = 1
            while m[i].isdigit():
                i += 1
            allTreeMolecules[root].append((m[i:],atoms, mstr[1:]))
    #print(allTreeMolecules)
    return allTreeMolecules, set(branchAtoms)


# list of atoms in this branch; when you cant add anything else, then you know the branch is done
# keep the list of atoms and the bonds and you pass the list of bonds to smiles and it does it for you
# look at winsotn's implementation and the visitor thing to figure out whether I can use it here
def getAllTreeBranches(branchEdgesAll, ringAtoms, m):
    # while there are still non-ring edges we have not used
    allTrees = []
    while len(branchEdgesAll) > 0:
        thisTreeNodeAtoms = []
        thisTreeNodeBonds = []
        # add the starting item for this tree:
        (s,neighbors) = branchEdgesAll.popitem()
        thisTreeNodeAtoms.append(s)
        thisTreeNodeAtoms.extend(neighbors)
        thisTreeNodeBonds.extend([m.GetBondBetweenAtoms(s,x).GetIdx() for x in neighbors])
        queue = set(copy.copy(neighbors)) # keep the elements in the queue that we wish to explore
        # keep adding all the edges that belong to this tree
        while (len(queue) > 0):
            elem = queue.pop()
            neibs = branchEdgesAll[elem]
            thisTreeNodeAtoms.append(elem)
            thisTreeNodeBonds.extend([m.GetBondBetweenAtoms(elem,x).GetIdx() for x in neibs])
            # remove these edges from branchEdgesByStart
            del branchEdgesAll[elem]
            for i in neibs:
                if i not in thisTreeNodeAtoms and i not in ringAtoms:
                    queue.add(i)
            # Moved after check to add to queue
            thisTreeNodeAtoms.extend(neibs)
        allTrees.append((set(thisTreeNodeAtoms),set(thisTreeNodeBonds)))
    # get all the molecular representations of the branches; allBranchMols is a list of tuples (u,v)
    # where u is the string representation of the branch and v is the list of atom indecies involved in that branch
    return allTrees

def getRingData(m, ringVectors, ringEdges):
    ringEdgesByStart = defaultdict(dict)
    ringEdgesByEnd = defaultdict(dict)

    # All ring edges implied by rings in ringVectors
    for ring in ringVectors:
        for i in range(len(ring)):
            s = ring[i]
            e = ring[(i+1) % len(ring)]
            ringEdgesByStart[ring][s] = e
            ringEdgesByEnd[ring][e] = s
            ringEdges.difference_update([(s,e),(e,s)])

    # All edges should be found by processing ringVectors
    assert(len(ringEdges) == 0)

    # OLD IMPLEMENTATION
    # while len(ringEdges) != 0:
    #     # NOTE: direction of edges may not be consistent around ring
    #     (s,e) = ringEdges.pop()
    #     for ring in ringVectors:
    #         if s in ring and e in ring:
    #             ringEdgesByStart[ring][s] = e
    #             ringEdgesByEnd[ring][e] = s
    #     if s in [0,1,6,7] and e in [0,1,6,7]:
    #         print(s,e)
    #         print("END: ",ringEdgesByEnd)
    #         print("START: ",ringEdgesByStart)
    return ringEdgesByStart, ringEdgesByEnd


# finds the ring representation in one direction (direction is set by the ringEdges)
def getTraversalForStartingPoint(m, startingIndex, ring, ringEdges, branchesByRoot):
    order = [startingIndex]
    currIndex = startingIndex
    traversal = ":"
    traversal = traversal + str(m.GetAtomWithIdx(currIndex).GetSymbol())
    if currIndex in branchesByRoot.keys():
        branchStrs = []
        for _,_,b in branchesByRoot[currIndex]:
            branchStrs.append(b)
        branchStrs.sort()
        for b in branchStrs:
            traversal = traversal + "(" + b + ")"
    currIndex = ringEdges[currIndex]
    #print("all trees: ", allTrees)
    # iterate traverse the ring until you circle back around, finding all the branches incident
    # on that point in the ring adding each one in parenths, and then add the atom itself
    #print(traversal, currIndex)
    prevIndex = currIndex
    while (currIndex != startingIndex):
        bond = ""
        if prevIndex != currIndex:
            bondType = m.GetBondBetweenAtoms(prevIndex, currIndex).GetBondType()
            if bondType == Chem.BondType.DOUBLE:
                bond = "="
            elif bondType == Chem.BondType.TRIPLE:
                bond = "#"
        order.append(currIndex)
        traversal = traversal + bond + str(m.GetAtomWithIdx(currIndex).GetSymbol())
        if currIndex in branchesByRoot.keys():
            branchStrs = []
            for _,_,b in branchesByRoot[currIndex]:
                branchStrs.append(b)
            branchStrs.sort()
            for b in branchStrs:
                traversal = traversal + "(" + b + ")"
        prevIndex = currIndex
        currIndex = ringEdges[currIndex]
    traversal = traversal + (":")
    bond = ""
    bondType = m.GetBondBetweenAtoms(prevIndex, startingIndex).GetBondType()
    if bondType == Chem.BondType.DOUBLE:
        bond = "="
    elif bondType == Chem.BondType.TRIPLE:
        bond = "#"
    traversal = traversal + bond

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

            allTraversals.append((ringTraversalForward,ring))
            allTraversals.append((ringTraversalReversed, ring))
            traversalToOrder[(ringTraversalForward, ring)] = atomOrderingF
            traversalToOrder[(ringTraversalReversed, ring)] = atomOrderingR

        allTraversals.sort(key=lambda t: t[0])
        minTraversal = allTraversals[0]
        stringToRingOrder[minTraversal] = traversalToOrder[minTraversal]
        completeMolecule.append(allTraversals[0])

    completeMolecule.sort(key=lambda t: t[0])
    #print(completeMolecule)

    moleculeOrderings = [stringToRingOrder[mol] for mol in completeMolecule]

    # since branches only appear once, there is redundance in overlapping rings and exposed endpoints of branches
    # sharedAtoms = []
    # seenAtoms = set()
    # for ring in ringVectors:
    #     for r in ring:
    #         if r in seenAtoms:
    #             sharedAtoms.append(r)
    #         else:
    #             seenAtoms.add(r)

    # for atom in branchAtoms:
    #     if atom in seenAtoms:
    #         sharedAtoms.append(atom)

    # once you use one atom in a branch, need to remove all other shared endpoints to avoid duplicate tagging

    tagDict = dict() # atom index to tag
    freqDict = defaultdict(lambda: 0)
    idxOrder = []

    # Initial string, with indices as tags for everything
    result = ""
    for j in range(len(moleculeOrderings)):
        traversal = moleculeOrderings[j]
        result = result + ":"
        prevIdx = traversal[0]
        for i in range(len(traversal)):
            atomIdx = traversal[i]
            idxOrder.append(atomIdx)
            freqDict[atomIdx] += 1
            bond = ""
            if prevIdx != atomIdx:
                bondType = m.GetBondBetweenAtoms(prevIdx, atomIdx).GetBondType()
                if bondType == Chem.BondType.DOUBLE:
                    bond = "="
                elif bondType == Chem.BondType.TRIPLE:
                    bond = "#"

            atom = m.GetAtomWithIdx(atomIdx).GetSymbol()
            result = result + bond + atom + str(atomIdx)
            #result = tagIfShared(result, atomIdx, sharedAtoms, tagDict, branch, treesByAtoms)

            if atomIdx in branchesByRoot.keys():
                branchStrs = []
                branchAlists = []
                for b,alist,_ in branchesByRoot[atomIdx]:
                    branchStrs.append(b)
                    branchAlists.append(alist)
                branchStrs.sort()
                for b in branchStrs:
                    result = result + "(" + b + ")"
                    i = 0
                    idxStr = ""
                    while i < len(b):
                        if b[i].isalpha() and len(idxStr) > 0:
                            idx = int(idxStr)
                            idxOrder.append(idx)
                            freqDict[idx] += 1
                            idxStr = ""
                        elif b[i].isdigit():
                            idxStr += b[i]

                        i += 1
                    if len(idxStr) > 0:
                        idx = int(idxStr)
                        idxOrder.append(idx)
                        freqDict[idx] += 1
                        idxStr = ""



                for alist in branchAlists:
                    for a in alist:
                        if a in branchesByRoot:
                            del branchesByRoot[a]
            prevIdx = atomIdx
        bond = ""
        bondType = m.GetBondBetweenAtoms(prevIdx, traversal[0]).GetBondType()
        if bondType == Chem.BondType.DOUBLE:
            bond = "="
        elif bondType == Chem.BondType.TRIPLE:
            bond = "#"

        result = result + ":" + bond
        if j != len(moleculeOrderings)-1:
            result = result + "-"

    print(result)
    curTag = 0
    for idx in idxOrder:
        if freqDict[idx] > 1 and idx not in tagDict:
            tagDict[idx] = curTag
            curTag += 1
    print(tagDict)
    minTagResult = ""
    idxStr = ""
    i = 0
    while i < len(result):
        if result[i].isdigit():
            idxStr += result[i]
        else:
            if len(idxStr) > 0:
                idx = int(idxStr)
                if idx in tagDict:
                    minTagResult += str(tagDict[idx])
                idxStr = ""
            minTagResult += result[i]

        i += 1
    if len(idxStr) > 0:
        idx = int(idxStr)
        if idx in tagDict:
            minTagResult += str(tagDict[idx])
        idxStr = ""


    return minTagResult

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
        ringAtoms = set()
        for ring in ringVectors:
            for r in ring:
                ringAtoms.add(r)

        ringEdges, branchEdgesAll = initializeRingData(m, ringAtoms)
        allTreeBranches = getAllTreeBranches(branchEdgesAll, ringAtoms, m)
        ringEdgesByStart, ringEdgesByEnd = getRingData(m, ringVectors, ringEdges)
        branchesByRoot, branchAtoms = moleculesFromBranches(allTreeBranches, m, ringAtoms)
        # ringEdgesByStart/End -- gives you the edges and their directions as dictionaries
        result = getRingTraversals(m, branchesByRoot, ringEdgesByStart, ringEdgesByEnd)
        return result
    else:
        return s
        # Output SMILES

if __name__ == "__main__":
    # x = input("Please enter a valid SMILES encoding of a molecule.")
    # encodeSMILES(x)
    cubane = 'C12C3C4C1C5C4C3C25'
    bicyclohexyl = 'C1CCCCC1C1CCCCC1'
    cyclohexane = 'C1CCCCC1'
    #Changed into cycle to test branch code
    _3_propyl_4_isopropyl_1_heptene = 'C1=CC(CCC)C(C(=C)C)CCC1'
    testMultBranchSameRoot1 = 'C1=CC(CCC)(C(=C)C)CCC1'
    testMultBonds1 = 'C1=C=C=C=C=C=C=1'
    _3_ringed_multiple_branches = 'C1(C2=CC(F)=C(C=C2N(C2CC2)C=C1C(=O)O)N1CCNCC1)=O'
    _1_methyl_3_bromo_cyclohexene_1 = 'CC1=CC(CCC1)Br'
    testSymHack1 = 'C1OSN1CCCCCC1PSN1'
    testSymHack2 = 'C1CCC1C1CCC1'
    testCurr = cubane
    #print(testCurr)
    print(encodeSMILES("CC1SC2NCNC(NCC3CCCO3)C2C1C"))
