#!/usr/bin/env python
import argparse
import numpy as np
import gzip
def get_args():
    parser = argparse.ArgumentParser(description="provides by base pos statistics for fastq input")
    parser.add_argument("-f", "--filename", help="fastq file", required=True)
    return parser.parse_args()

args = get_args()
inFile = args.filename

def convert_phred(character):
    '''Input quality character. Ouput phred33 score.'''
    return (ord(character)-33)

#array to store qscores for each location
qualArray = np.zeros((1,101))
with gzip.open(inFile, "rt") as file:
    lineCount = 1
    #counts number of reads every four lines
    readCount = 0
    for line in file:

        if lineCount == 4:
            #If quality line
            readCount += 1
            lineCount = 1
            for charIndex in range(len(line.strip())):
                #Add quality score for each location
                qualArray[0,charIndex] += convert_phred(line[charIndex])
        else:
            lineCount += 1
    #calculate average quality for each position
    averageScores = qualArray/readCount

    import matplotlib.pyplot as plt
    plt.bar(range(101), averageScores[0,:])
    plt.title(inFile)
    plt.xlabel("Base Positon")
    plt.ylabel("Mean Quality Score (Phred33)")
    plt.savefig(inFile + ".png")
