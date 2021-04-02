from rdkit import Chem
import copy
from collections import defaultdict

def initializeRingData(m):
    branchEdgesByStart = dict()
    branchEdgesByEnd = dict()
    ringEdges = set()

    for bond in m.GetBonds():
        # This encoding does not support aromatic rings so request a new molecule.
        if str(bond.GetBondType()) == "AROMATIC":
            x = input("This molecule is aromatic. Please try again.\n")
            return encodeSMILES(x)

        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bond.IsInRing():     
            ringEdges.add((start,end))           
        else:
            branchEdgesByStart[start] = end
            branchEdgesByEnd[end] = start

    return branchEdgesByStart, branchEdgesByEnd, ringEdges
    
def getJoinedBranches(branchEdgesByStart,branchEdgesByEnd):
    (s,e) = branchEdgesByStart.popitem()
    del branchEdgesByEnd[e]
    joinedBranches = [[s,e]]
    while len(branchEdgesByStart) != 0:   
        start = joinedBranches[-1][0]
        end = joinedBranches[-1][-1]

        if end in branchEdgesByStart:
            joinedBranches[-1].append(branchEdgesByStart[end])
            del branchEdgesByEnd[branchEdgesByStart[end]]
            del branchEdgesByStart[end]
        elif start in branchEdgesByEnd:
            joinedBranches[-1] = [branchEdgesByEnd[start]] + joinedBranches[-1]
            del branchEdgesByStart[branchEdgesByEnd[start]]
            del branchEdgesByEnd[start]
        else:
            (s,e) = branchEdgesByStart.popitem()
            del branchEdgesByEnd[e]
            joinedBranches.append([s,e])

    return joinedBranches

def getRingData(m, ringVectors):
    branchEdgesByStart, branchEdgesByEnd, ringEdges = initializeRingData(m)
    joinedBranches = getJoinedBranches(branchEdgesByStart, branchEdgesByEnd)

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
    for branch in joinedBranches:
        for ring in ringVectors:
            if (branch[0] in ring):
                branchRev = branch[::-1]
                ringBranches[ring].append(branchRev)
                if len(branchRev) > len(ringLongestBranch[ring]):
                    ringLongestBranch[ring] = branchRev
            elif (branch[-1] in ring):
                ringBranches[ring].append(branch)
                if len(branch) > len(ringLongestBranch[ring]):
                    ringLongestBranch[ring] = branch
        
    minTraversals = defaultdict(list)
    for ring in ringVectors:
        traversal1 = copy.copy(ringLongestBranch[ring])
        traversal2 = copy.copy(ringLongestBranch[ring])

        while (len(traversal1) < (len(ring)+1)):
            traversal1.append(ringEdgesByStart[ring][traversal1[-1]])
            traversal2.append(ringEdgesByEnd[ring][traversal2[-1]])
        
        minTraversals[ring] = min(traversal1,traversal2) # FIX: After taking min, need to embed remaining branches from ringBranches when converting to string
    return minTraversals

# Converts SMILES to canonical SMILES
def encodeSMILES(s):
    m = Chem.MolFromSmiles(s)

    # Asks for input again if molecule not valid SMILES
    if m == None:
        x = input("Please try again.\n")
        return encodeSMILES(x)              
    
    ri = m.GetRingInfo()                    
    n = ri.NumRings()

    # If there are any rings in structure, we need to perform encoding so 
    # molecule can be parsed using context-free grammar.
    if n > 0:
        ringVectors = ri.AtomRings()
        joinedBranches, ringEdgesByStart, ringEdgesByEnd = getRingData(m, ringVectors)
        minRingTraversals = getRingTraversals(joinedBranches, ringEdgesByStart, ringEdgesByEnd)



# x = input("Please enter a valid SMILES encoding of a molecule.")
# encodeSMILES(xs)
cubane = 'C12C3C4C1C5C4C3C25'
bicyclohexyl = 'C1CCCCC1C1CCCCC1'
encodeSMILES(bicyclohexyl)
