#!/usr/bin/env python3
import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="Demultiplexes Paired-End Reads into separate fastq.gz files. Outputs fastq files (R1 & R2) for: each barcode [barcode.fq], reads with different indexes [mistmatched.fq], & reads of low quality (below input threshold)[unknown.fq]")
    parser.add_argument("-f", "--forwardRead", help="forward read fastq", required=True)
    parser.add_argument("-r", "--reverseRead", help="reverse read fastq", required=True)
    parser.add_argument("-b", "--forwardBarcodes", help="Foward index fastq", required=True)
    parser.add_argument("-B", "--reverseBarcodes", help="Reverse index fastq", required=True)
    parser.add_argument("-q", "--minimumQuality", help="minimum average quality score for indexes (>q goes into unknown.fq)", required=True)
    parser.add_argument("-o", "--output", help="Output directory path/ (MUST BE EMPTY & ALREADY EXIST)", required=True)
    parser.add_argument("-m", "--mappingFile", help="Tsv file of indexes. Two columns: Symbol    Barcode", required=True)


    return parser.parse_args()

args = get_args()
forInsertFile = args.forwardRead
revInsertFile = args.reverseRead
forIndexFile = args.forwardBarcodes
revIndexFile = args.reverseBarcodes
minQual = args.minimumQuality
outDirectory = args.output
mapFile = args.mappingFile

##### READ IN BARCODES
barcodeDict = dict()

with open(mapFile, "r") as mapInput:
    for line in mapInput:
        cols = line.split(" ")
        barcodeDict.setdefault(cols[1].strip(), cols[0])

#####  DEFINE FUNCTIONS
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
        elif base == "N":
            compBase = "N"
        else:
            print("isReverseComplement: not a normal DNA seq")
            return False

        if querySeq[queryIndex] != compBase:
            #at least one base is off. not reverse complement
            return False
        queryIndex +=1
    #all bases match. is reversed complement
    return True

def getReadQuality(readSeq):
    '''Calculates average read base quality in Phred33'''
    qScoreSum = 0
    baseCount = 0
    for base in readSeq:
        qScoreSum += ord(base)-33
        baseCount += 1

    readQuality = qScoreSum/baseCount
    return readQuality

##### CHECK IF FILES ARE ZIPPED
if "gz" in forInsertFile:
    zipped = True
else:
    zipped = False
    print("gzip your input fastq's")

##### OPEN FILES WITH OPEN OR GZIP
if zipped:
    with gzip.open(forInsertFile, "rt") as forInsert, gzip.open(revInsertFile, "rt") as revInsert, gzip.open(forIndexFile, "rt") as forIndex, gzip.open(revIndexFile, "rt") as revIndex:
        #Counters for reads based on index matching/quality
        uknownCount = 0
        mismatchedCount = 0
        totalReadCount = 0

        #Counter for line within each record (4 lines per record)
        lineCount = 1

        for forInsertLine in forInsert:
            if lineCount ==1:
                #AT THE INSTANCE OF A NEW RECORD:

                #save headers as their own variable
                forHeader = forInsertLine.strip()
                revHeader = revInsert.readline().strip()

                #start assuming nothing is "wrong" with the read
                indexUnknown = False
                mismatch = False

                #don't need index headers. Skip these lines
                forIndex.readline()
                revIndex.readline()

            elif lineCount == 2:
                #SEQUENCE LINE
                #check is indexes are mismatched
                forIndexSeq = forIndex.readline().strip()
                revIndexSeq = revIndex.readline().strip()

                if isReverseComplement(forIndexSeq, revIndexSeq):
                    pass
                else:
                    mismatchedCount += 1
                    mismatch = True

                #create index header titles
                indexCombo = forIndexSeq + "-" + revIndexSeq
                forwardRead = forHeader + " " + indexCombo + "\n"
                reverseRead = revHeader + " " + indexCombo + "\n"

                #save insert sequences
                forwardRead += forInsertLine
                reverseRead += revInsert.readline()


            elif lineCount == 3:
                #saving plus line in inserts and moving line for barcodes
                forwardRead += forInsertLine
                reverseRead += revInsert.readline()
                forIndex.readline()
                revIndex.readline()

            elif lineCount == 4:
                #QUALITY LINE AND END OF RECORD
                #check read quality and compare it to miniumum
                forIndexQ = forIndex.readline().strip()
                revIndexQ = revIndex.readline().strip()

                if getReadQuality(forIndexQ) < int(minQual) or getReadQuality(revIndexQ) < int(minQual):
                    indexUnknown = True
                    uknownCount += 1

                #Determining whether or not quality barcode pair is in mapping file
                if not indexUnknown and not mismatch:
                    if forIndexSeq in barcodeDict:
                        #symbol in mapping file will be used as output
                        barcodeSymbol = barcodeDict[forIndexSeq]
                    else:
                        indexUnknown = True
                        uknownCount += 1



                #save insert quality lines. FULL READS ARE NOW STORED
                forwardRead += forInsertLine
                reverseRead += revInsert.readline()

                ##### write forward and reverse reads to correct output files
                #name output file
                if mismatch:
                    outputFileName = "mismatch"
                elif indexUnknown:
                    outputFileName = "unknown"
                else:
                    outputFileName = barcodeSymbol

                sendForTo = outDirectory + outputFileName + "_R1.fq"
                sendRevTo = outDirectory + outputFileName + "_R2.fq"

                with open(sendForTo, "a") as forOutput, open(sendRevTo, "a") as revOutput:
                    forOutput.write(forwardRead)
                    revOutput.write(reverseRead)

                ##### END OF RECORD
                lineCount = 0
                totalReadCount += 1

            ##### END OF LINE
            lineCount += 1

print("DEMULTIPLEXING COMPLETE")
print("\t Total Reads Demultiplexed: " + str(totalReadCount))
print("\t # of Reads unkown reads (low quality OR not in mapping File): " + str(uknownCount))
print("\t # of Reads with mistmatched Indexes: " + str(mismatchedCount))
print("\t # barcoded outputFiles: " + str(len(barcodeDict)))
