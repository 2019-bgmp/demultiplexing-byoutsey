#!/usr/bin/env python3
def isReverseComplement(templateSeq, querySeq):
    '''Takes two DNA Sequences (No N's). If they are reverse complements, function will return True. Otherwise, False.'''
    queryIndex = 0
    for base in templateSeq[len(templateSeq)::-1]:
        if base == "G":
            compBase = "C"
        elif base == "C":
            compBase = "G"
        elif base == "T":
            compBase = "A"
        elif base == "A":
            compBase = "T"
        else:
            print("isReverseComplement: not a normal DNA seq")
            return False

        if querySeq[queryIndex] != compBase:
            #at least one base is off. not reverse complement
            return False
        queryIndex +=1
    #all bases match. is reversed complement
    return True

testTemp = "AAAAAGGGGGCCCCCTTTTT"
testQuery = "AGGGGGGGAAACCCCCC"

#Should print False
print(isReverseComplement(testTemp, testQuery))

testTemp = "AAAAAGGGGGCCCCCTTTTT"
testQuery ="TTTTTCCCCCGGGGGAAAAA"
testQuery = testQuery[len(testQuery)::-1]

#should print True
print(isReverseComplement(testTemp, testQuery))
