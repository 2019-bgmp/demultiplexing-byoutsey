#!/usr/bin/env python3

def getReadQuality(readSeq):
    '''Calculates average read base quality in Phred33'''
    qScoreSum = 0
    baseCount = 0
    for base in readSeq:
        qScoreSum += ord(base)-33
        baseCount += 1

    readQuality = qScoreSum/baseCount
    return readQuality

#output should be 32
print(getReadQuality("AAAAAAAAAA"))
