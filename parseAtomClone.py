from arpeggio import Optional, ZeroOrMore, OneOrMore, EOF, StrMatch
from arpeggio import RegExMatch as _
from arpeggio import ParserPython
import sys

# From Arpeggio tutorial, prevents string from appearing in parse tree
class SuppressStr(StrMatch):
    suppress = True

# Parser Functions
def atomUnit():
    # Second option lists accepted atoms, can add more if necesary
    return [(SuppressStr('('), atomList, SuppressStr(')')),
            (_(r'[BCHNOPS]'), Optional(_(r'\d*'))) ]

def atomList():
    return ZeroOrMore(atomUnit)

def ring():
    return SuppressStr(':'), atomUnit, atomList, SuppressStr(':')

def fragment():
    return atomList, ring

def fragmentList():
    return fragment, ZeroOrMore(SuppressStr('-'),fragment), EOF

def initParser():
    """
    Set parser global variable to arpeggio parser for our format
    """
    global parser
    parser = ParserPython(fragmentList)

def strToAdjacencyMatrix(s):
    """
    Parses string and returns adjacency matrix and label array

    arguments:
    s - string to be parsed

    returns:
    adjacency matrix as list of lists, label matrix as list of strings
    """
    tree = parser.parse(s)
    print(tree.tree_str())
    return tree

def printMatrix(mat):
    """
    Prints matrix (list of lists) in nice format
    (assumes all entries take same numeber of characters)

    arguments:
    mat - list of lists representing matrix
    """
    for row in mat:
        print(' '.join([str(x) for x in row]))

assert(len(sys.argv) == 2)

initParser()


with open(sys.argv[1]) as f:
    for line in f:
        cleanLine = line.strip()
        print("String: ", cleanLine)
        matrix = strToAdjacencyMatrix(cleanLine)
        print(matrix)
