#!/usr/bin/env python2.7

##############################################################################################################
# RATTwithGFF.py                                                                                             #
# Author: Preston Basting                                                                                    #
# Last Changed 06/21/18                                                                                      #
# Contact: pjb68507@uga.edu                                                                                  #
# Function: This program is designed to run the rapid annotation transfer tool (RATT) using gff as input.    #
#   RATT is a tool for transferrin annotations from one genome to another. Normally, RATT requires embl as   #
#   the annotation format. This program uses other scripts to handle the file conversions and corrects for   #
#   the errors that often occur. Produces  gff files as output.                                              #
##############################################################################################################

import subprocess
import sys
import timeit

def main():
    startTime = timeit.default_timer()
    if (len(sys.argv) != 6):
        sys.stderr.write('<RATTwithGFF.py>USAGE ERROR: SetupEmbls.py <reference gff> <reference fasta> <query fasta> <runID> <Transfer type>\n')
        return 0

    gffFileName = sys.argv[1]
    fastaFileName = sys.argv[2]
    queryFastaFile = sys.argv[3]
    sampleID = sys.argv[4]
    rattType = sys.argv[5]
    contigs = []

    if (validArgs(gffFileName, fastaFileName, queryFastaFile, rattType)):
        
        #ensures end of line is LF not CRLF
        #creates temp files with correct line endings
        gffFileName = fixLineEndings(gffFileName)
        fastaFileName = fixLineEndings(fastaFileName)
        queryFastaFile = fixLineEndings(queryFastaFile)
        try:
            splitGenomicFiles(gffFileName, fastaFileName, contigs)
        except OSError:
            sys.stderr.write("***********************************************************************************\n")
            sys.stderr.write("<RATTwithGFF.py> FATAL-ERROR\n")
            sys.stderr.write('OSError: Could not parse reference files, check that samtools is installed properly\n')
            sys.stderr.write("***********************************************************************************\n")
            return 0
        
        try:
            gffsToEmbls(contigs)
        except OSError as e:
            sys.stderr.write("***********************************************************************************\n")
            sys.stderr.write("<RATTwithGFF.py> FATAL-ERROR\n")
            sys.stderr.write('OSError: Could not convert gff to embl, check that EMBLmyGFF is installed properly\n')
            sys.stderr.write("***********************************************************************************\n")
            print str(e)
            return 0
        
        global ratt_dir 
        
        ratt_dir = sampleID+"_RATT" #this is the directory that the ratt results will be stored
        
        ratt_worked = runRatt(fastaFileName,queryFastaFile, sampleID, rattType)
        if (ratt_worked):
            processResults_worked = processRattResults(gffFileName)
            if (processResults_worked):
                printTransferStats()
                subprocess.call(["rm",gffFileName, fastaFileName, queryFastaFile]) #removes temporary files with fixed line endings
                elapsedTime = (timeit.default_timer() - startTime)/60
                print("\n<RATTwithGFF.py> COMPLETE in "+format(elapsedTime,'.2f')+" minutes")
    else:
        return 0

#checks the input arguments to make sure they are vaid
#called by main
def validArgs(gff, fasta1, fasta2, transType):
    print "\n<RATTwithGFF.py> checking for valid input files..."
    with open(gff, "r") as gffFile:
        line = gffFile.readline()
        if (line.find("#gff-version 3") == -1):
            sys.stderr.write('ERROR: first argument is not a gff file\n')
            return False

    with open(fasta1, "r") as fastaFile:
        line = fastaFile.readline()
        if (line.find('>') == -1):
            sys.stderr.write('ERROR: second argument is not a fasta file\n')
            return False
    
    with open(fasta2, "r") as fastaFile:
        line = fastaFile.readline()
        if (line.find('>') == -1):
            sys.stderr.write('ERROR: second argument is not a fasta file\n')
            return False
    
    #these are the options for RATT transfer type
    validTypes = ["Assembly", "Assembly.Repetitive", "Strain", "Strain.Repetitive", "Species", "Species.Repetitive", "Multiple"]
    if transType in validTypes:
        print "\n<RATTwithGFF.py> input files are valid"
        return True
    else:
        sys.stderr.write("ERROR: '"+transType+"' is not a valid RATT transfer type\n")
        sys.stderr.write('VALID TRANSFER TYPES: Assembly, Assembly.Repetitive, Strain, Strain.Repetitive, Species, Species.Repetitive, Multiple\n')
        return False

#checks line ending format
#creates a temp file with LF line endings
def fixLineEndings(fileName):
    newFileName = "temp_"+fileName
    with open(newFileName,"wb") as outFile:
        with open(fileName) as file:
            for line in file:
                if line.find("\r") != -1:
                    newLine = line.replace("\r","")
                else:
                    newLine = line    
                outFile.write(newLine)
    return newFileName



#converts a genomic gff into seperate gffs for each contig
#converts a genoic fasta into a fasta for each contig
#called by main
def splitGenomicFiles(gff, fasta, contigNames):
    print "\n<RATTwithGFF.py> Indexing fasta.."
    subprocess.call(["samtools","faidx",fasta])
    index = fasta+".fai"
    
    with open(index,"r") as fai:
        for line in fai:
            contigNames.append(line.split("\t")[0]) #makes a list of all contigs from fasta index
    
    subprocess.call(["mkdir","contig_gff"])
    subprocess.call(["mkdir","contig_fasta"])
    subprocess.call(["mkdir","contig_embl"])

    print "\n<RATTwithGFF.py> Generating contig fastas and gffs..."
    for contig in contigNames:
        output = subprocess.check_output(["samtools","faidx", fasta, contig]) #gets sequence specific for a contig
        
        with open("contig_fasta/"+contig+".fa", "w") as contigFasta:
            contigFasta.write(output)
    
        with open(gff,"r") as gffFile:
            with open ("contig_gff/"+contig+"_tmp1.gff","w") as contigGff:
                contigGff.write("##gff-version 3\n")
                for line in gffFile:
                    if (line.find(contig) != -1):
                        contigGff.write(line)
        
        addInfoToGff(contig, "contig_gff/"+contig+"_tmp1.gff")
        subprocess.call(["rm","contig_gff/"+contig+"_tmp1.gff"])

#adds ncRNA_class attribute to ncRNA features
#adds notes that include the original parents
#called by main        
def addInfoToGff(contig, gff):
    print "\n<RATTwithGFF.py> adding ncRNA_class=other to "+contig+".gff"
    output = 'contig_gff/'+contig+".gff"
    with open(gff) as gffFile:
        with open(output, "w") as outputFile:
            for line in gffFile:
                #adds ncRNA_class attribute to ncRNA
                if( line.find("\tncRNA\t") != -1 and line.find("ncRNA_class") == -1 ):
                    line = line[:len(line)-1]
                    line = line + ";ncRNA_class=other\n"
                #adds the parent attribute as a note
                #this makes sure the parent attribute isn't removed during file conversions
                if( line.find("Parent=") >=0 ):
                    line = line[:len(line)-1]
                    parent = line[line.find("Parent=")+7:]
                    if ( parent.find(";") != -1):
                        parent = parent[:parent.find(";")]
                    else:
                        parent = parent[:len(parent)]
                    line = line + ";note=Parent:"+parent+"\n"
                outputFile.write(line)

#calls EMBLmyGFF3.py to convert the contig gff to an embl file
#called by main
def gffsToEmbls(contigNames):
    for contig in contigNames:
        print "\n<RATTwithGFF.py> converting "+contig+".gff to "+contig+".embl...."
        subprocess.call(["EMBLmyGFF3",\
                            "contig_gff/"+contig+".gff",\
                            "contig_fasta/"+contig+".fa",\
                            "-o", "contig_embl/"+contig+"_tmp1.embl",\
                            "-i", "tag",\
                            "-p", contig,\
                            "-s", "unknown",\
                            "-t", "linear",\
                            "-d", "STD",\
                            "-m", "genomic DNA",\
                            "-x", "UNC",\
                            "--rg","none",\
                            "-r","1",\
                            "-a",contig,\
                            "--keep_duplicates",\
                            "-q",\
                            "--shame"])
        cleanEmbl(contig)

#checks if str2 is found in str1
#called by isBrokenLine
def containsStr(str1, str2):
    if ( str1.find(str2) == -1 ):
        return False
    else:
        return True

#checks if the line is wrapped info from the previous line
#called by fixBrokenLines
def isBrokenLine(line):
    if (containsStr(line,"FT                   ") and not containsStr(line, "/")):
        return True
    else:
        return False

#conbines two embl lines
#called by fixBrokenLines
def combineLines(line1,line2):
    first = line1[:line1.find("\n")]
    second = line2[21:]
    return (first+second)

#checks if a value is in a list
#called by fixBrokenLines
def isInList(val,alist):
    for i in alist:
        if( i == val ):
            return True
    return False

#checks an embl file for lines that are wrapped/broken
#takes the wrapped line and adds it to the line above
#this is necessary as RATT does not recognize the wrapped line as being a part of the line above
#called by cleanEmbl
def fixBrokenLines(fileName):
    toRemove = []
    fixedLines = []
    with open(fileName) as file:
        lines = file.readlines()
        for x in range(len(lines)-1,-1,-1): #starts from last line and moves up ##necessary to fix lines broken into 3+ lines
            if isBrokenLine(lines[x]):
                lines[x-1] = combineLines(lines[x-1], lines[x])
                toRemove.append(x) #keeps track of lines to be removed
    for i in range(0,len(lines)):
        if( isInList(i, toRemove) == False):
            fixedLines.append(lines[i])
    return fixedLines

#fixes broken/wrapped lines in the embl files for each contig and writes to new file
#called by gffsToEmbls
def cleanEmbl(contig):
    print "\nfixing line-breaks for: "+contig+".embl...."
    embl = "contig_embl/"+contig+"_tmp1.embl"
    cleanEmblLines = fixBrokenLines(embl)
    with open("contig_embl/"+contig+".embl","w") as outputFile:
        for line in cleanEmblLines:
            outputFile.write(line)

    subprocess.call(["rm",embl])

#calls the RATT script to transfer the annotations
#also creates directories to organize the RATT output files
#called by main
def runRatt(subFa,queryFa, sampleID, parameter):
    print "\n<RATTwithGFF.py>***RUNNING QUAST on both sequences....."
    subprocess.call(["mkdir",ratt_dir]) #makes directory where the ratt results are stored
    try:
        subprocess.call(["quast.py","-o",ratt_dir+"/ref_quast","--fast","-s","--silent",subFa])
        subprocess.call(["quast.py","-o",ratt_dir+"/query_quast","--fast","-s","--silent",queryFa])
    except OSError:
        sys.stderr.write("***********************************************************************************\n")
        sys.stderr.write("<RATTwithGFF.py> ERROR\n")
        sys.stderr.write("OSError: Could not call QUAST, Check that it is installed correctly\n")
        sys.stderr.write("***********************************************************************************\n")
#         return False

    print "\n<RATTwithGFF.py>***RUNNING RATT....."
    try:
        subprocess.call(["start.ratt.sh","../contig_embl","../"+queryFa, sampleID, parameter],cwd=ratt_dir)
    except OSError:
        sys.stderr.write("***********************************************************************************\n")
        sys.stderr.write("<RATTwithGFF.py> FATAL-ERROR\n")
        sys.stderr.write("OSError: could not call RATT, check that RATT is installed correctly\n")
        sys.stderr.write("***********************************************************************************\n")
        return False
    
    subprocess.call(["mkdir","final_embl","Report_gff","Report_txt","NOTTransfered_embl","nucmer","tmp2_embl","uncorrected_embl","final_gff"],cwd=ratt_dir)
    subprocess.call(["mv *.final.embl final_embl"], shell=True, cwd=ratt_dir)
    subprocess.call(["mv *.Report.gff Report_gff"], shell=True, cwd=ratt_dir)
    subprocess.call(["mv *.Report.txt Report_txt"], shell=True, cwd=ratt_dir)
    subprocess.call(["mv *.NOTTransfered.embl NOTTransfered_embl"], shell=True, cwd=ratt_dir)
    subprocess.call(["mv nucmer.* nucmer"], shell=True, cwd=ratt_dir)
    subprocess.call(["mv *tmp2.embl tmp2_embl"], shell=True, cwd=ratt_dir)
    subprocess.call(["mv *.embl uncorrected_embl"], shell=True, cwd=ratt_dir)
    return True


#converts the RATT results from EMBL to gff using EMBOSS seqret
#fixes the errors that emberge due to the conversion
#generates a new genomic gff by combining the annotations from all contigs
#calculates the transfer stats for each feature
def processRattResults(origGff):
    print "\n<RATTwithGFF.py> processing RATT results..."
    genomicGff = ratt_dir+"/final_gff/genomic.final.gff"
    with open(genomicGff,"w") as genomic:
        genomic.write("##gff-version 3\n")

    embls = subprocess.check_output(["ls", ratt_dir+"/final_embl"]).splitlines()
    for embl in embls:
        outFile = ratt_dir+"/final_gff/"+embl[:embl.find(".embl")]+".gff"
        temp = "temp01.gff"
        try:
            emblToGff(ratt_dir+"/final_embl/"+embl)
        except OSError:
            sys.stderr.write("***********************************************************************************\n")
            sys.stderr.write("<RATTwithGFF.py> FATAL-ERROR\n")
            sys.stderr.write("OSError: Could not call EMBOSS seqret, make sure EMBOSS is installed correctly\n")
            sys.stderr.write("***********************************************************************************\n")
            return False

        print "\n<RATTwithGFF.py> Fixing embl-gff conversion errors in: "+outFile
        gffLines = parseGff(temp)
        gffLines = fixBiologicalRegions(gffLines)
        gffLines = cleanChromAndSource(gffLines)
        gffLines = renumberIDs(gffLines)
        gffLines = addAllParents(gffLines)
        gffLines = fixCdsPhase(gffLines)
        gffLines = cleanAttributes(gffLines)
        gffLines = fixCdsPos(gffLines)
        writeToFile(gffLines, outFile)
        addToGenomicGff(outFile,genomicGff)
    
    makeTransferStats(origGff,genomicGff)
    return True

#calls seqret function to convert embl to gff
def emblToGff(fileName):
    print "\n<RATTwithGFF.py> converting: "+fileName+" to gff..."
    subprocess.call(["seqret", "-sequence", fileName, "-feature",\
    "-fformat", "embl", "-fopenfile", fileName, "-osformat", "gff",\
    "-osname", "temp01", "-auto"])

# [0]:chrom [1]:source [2]:feature [3]:start [4]:end [5]:score [6]:strand [7]:phase [8]:info
#seperates gff into a list[lines] of lists[line values]
def parseGff(fileName):
    allParsedLines = []
    parsedLine = []
    with open(fileName) as file:
        for line in file:
            parsedLine = line.split("\t")
            allParsedLines.append(parsedLine) #adds each parsed line to a list
    subprocess.call(["rm", fileName])
    return allParsedLines

#finds the ID attribute in the given line
def getID(line):
    idString = line[line.find("ID=")+3:]
    idString = idString[:idString.find(";")]
    return idString

#finds and returns the Parent attribute in the given line
def getParent(line):
    parent = line[line.find("Parent=")+7:]
    if (parent.find(";") != -1):
        parent = parent[:parent.find(";")]
    else:
        parent = parent[:len(parent)-1]
    return parent

#takes a line and adds the given idString to the line as a parent attribute
def addParent(idString, line):
    front = line[:line.find(";")+1]
    back = line[line.find(";")+1:]
    newString = front+"Parent="+idString+";"+back
    return newString

#reformats 'biological_region' features to be consistent with the standard gff format
#these 'biological regions' feature are a consequence of using EMBOSS seqret embl->gff conversion 
def fixBiologicalRegions(lines):
    x = 0
    newLines = []
    while (x < len(lines)):
        #I have to check length of lines to distinguish between header and actual feature lines
        if((len(lines[x]) == 9) and (lines[x][2] == "biological_region") and 
        (lines[x][8].find("featflags=type:CDS") != -1)):
            bio_reg = x
            x+=1 #skips biological_region
            while (len(lines[x]) == 9 and lines[x][2] == "CDS"):
                lines[x][8] = lines[bio_reg][8] #assigns biological_region attributes to CDS
                newLines.append(lines[x])
                x+=1

        elif ((len(lines[x]) == 9) and (lines[x][2] == "biological_region") and 
            (lines[x][8].find("featflags=type:mRNA") != -1)):
            bio_reg = x
            lines[x][2] = "mRNA"
            newLines.append(lines[x])
            x+=1
            while (len(lines[x]) == 9 and lines[x][2] == "mRNA" and lines[x][8].find("ID=") == -1):
                x+=1 #doesnt add the split up mRNA features unless it is a seperate mRNA from biological feature

        elif (len(lines[x]) == 9) and (lines[x][2] == "databank_entry"):
            x+=1 #doesnt add 'databank_entry' feature

        else:
            newLines.append(lines[x])
            x+=1
    return newLines

#generates new IDs to replace the old ones
def renumberIDs(lines):
    count = 1
    for x in range(0,len(lines)):
        if (len(lines[x]) == 9):
            oldID = lines[x][8][:lines[x][8].find(";")] #assumes ID is first attribute
            lines[x][8] = lines[x][8][lines[x][8].find(";"):] #removes ID attribute
            newID = lines[x][0]+"."+str(count) #generates new ID based on feature order
            lines[x][8] = "ID="+newID+lines[x][8]

            if (x+1 < len(lines)-1 and len(lines[x+1]) == 9):
                if(oldID != lines[x+1][8][:lines[x+1][8].find(";")]): #if the next feature ID is different
                    count+=1
    return lines

#calculates the phase value for a given feature
#must provide the phase of the previous feature
def getPhase(lines,index,prevPhase):
    if (lines[index][2] == "CDS"):
        if (lines[index][6] == "+"):
            size = int(lines[index-1][4]) - (int(lines[index-1][3])-1) - int(prevPhase)
            remainder = size%3

        else:
            size = int(lines[index+1][4]) - (int(lines[index+1][3])-1) - int(prevPhase)
            remainder = size%3

        if (remainder == 0):
            return str(0)
        
        else:
            return str(3-remainder)
    
    else:
        return "."

#using information stored from the original gff
#finds if the parent feature was transferred
#adds the approprate parent attribute to link the features
def addAllParents(lines):
    for x in lines:
        if (len(x) == 9):
            if (x[8].find("note=Parent:") != -1):
                parent = x[8][x[8].find("note=Parent:")+12:]
                if (parent.find(";") != -1):
                    parent = parent[:parent.find(";")]
                
                else:#parent is last feature
                    parent = parent[:len(parent)-1]
                
                #finds parent feature
                for y in lines:
                    if (len(y) == 9):
                        if(y[8].find("note=ID:") != -1 ):
                            noteID = y[8][y[8].find("note=ID:")+8:]
                            noteID = noteID[:noteID.find(";")]
                            if (noteID == parent):
                                parent = getID(y[8])
                                x[8] = addParent(parent,x[8]) #adds corrected parent ID
                                break
    return lines

                

#For all CDS features, the phase is calculated and set
def fixCdsPhase(lines):
    x=0
    while(x < len(lines)):
        if ( len(lines[x]) == 9 and lines[x][2] == "CDS"):
            if lines[x][6] == "+":
                phase = str(0)
                lines[x][7] = phase
                x+=1
                while (len(lines[x]) == 9 and lines[x][2] == "CDS"):
                    phase = getPhase(lines,x,phase)
                    lines[x][7] = phase
                    x+=1
            else: #CDS phase on opposite strand must be calculated bottom->top
                negCDS = []
                while (len(lines[x]) == 9 and lines[x][2] == "CDS"):
                    negCDS.append(x)
                    x+=1
                i=negCDS[len(negCDS)-1]
                phase = str(0)
                lines[i][7] = phase
                i-=1
                while (i >= negCDS[0]):
                    phase = getPhase(lines,i,phase)
                    lines[i][7] = phase
                    i-=1
        else:
            x+=1
    return lines            

#removes a given attribute from a given line
def removeAttribute(line, attribute):
    if (line.find(attribute) != -1):
        startPos = line.find(attribute)
        substr = line[startPos:]
        if (substr.find(";") != -1):
            endPos = substr.find(";")+startPos+1

        else:
            endPos = len(substr)+startPos

        newLine = line[:startPos] + line[endPos:] #makes new line excluding the given attribute
        if (newLine.find("\n") == -1):
            if(newLine[len(newLine)-1] == ";"):
                newLine = newLine[:len(newLine)-1] #removes tailing ";"
            newLine = newLine+"\n" #adds missing endline

        if newLine.find(attribute) != -1:
            newLine = removeAttribute(newLine, attribute) #runs recursively until all of the given attribute are removed
            return newLine
        else:
            return newLine
    
    else:
        return line


#removes all attributes that are not needed
def cleanAttributes(lines):
    for line in lines:
        if (len(line) == 9):
            line[8] = removeAttribute(line[8],"locus_tag")
            # line[8] = removeAttribute(line[8],"note=")
            line[8] = removeAttribute(line[8],"transl_table")
            line[8] = removeAttribute(line[8],"codon_start")
            line[8] = removeAttribute(line[8],"featflags")
            line[8] = line[8].replace("standard_name=","Name=")
            if (line[2] == "ncRNA"):
                line[8] = line[8].replace("ncrna_class=","ncRNA_class=")        
    return lines

#makes sure all CDS start and end positions are within mRNA start and end positions
def fixCdsPos(lines):
    for x in lines:
        if (len(x) == 9 and x[2] == "CDS"):
            parent = getParent(x[8])
            for y in lines:
                if (len(y) == 9 and y[2] == "mRNA" and getID(y[8]) == parent): #finds the mRNA parent
                    if (x[3] < y[3]): #makes sure the CDS start and end are within the mRNA start and end
                        x[3] = y[3]
                    if (x[4] > y[4]):
                        x[4] = y[4]
                    break
    return lines




#removes the source column info
#cleans the chromosome column to only include the contig id
def cleanChromAndSource(lines):
    for line in lines:
        if (len(line) == 9):
            line[0] = line[0][:line[0].find(".final")]
            line[0] = line[0][line[0].find(".")+1:]
            line[1] = "."
    return lines

#outputs reformatted lines to file
def writeToFile(lines,outFileName):
    with open(outFileName, "w") as outFile:
        for line in range(0,len(lines)):
            newLine = "\t".join(lines[line])
            if (newLine.find("\tgap\t") ==-1): #removes gap features, it doesn't make sense to transfer gap annotations between assemblies
                outFile.write(newLine)

#adds the annotations from a contig gff to the genomic gff
def addToGenomicGff(inputGff,genomicGff):
    print "\nadding: "+inputGff+" annotations to: "+genomicGff
    with open(inputGff, 'r') as inFile:
        with open(genomicGff, 'a') as outFile:
            for line in inFile:
                if line.find("##FASTA") != -1:
                    break
                else:
                    if (line.find("#") == -1):
                        outFile.write(line)

#counts the amount of unique features in each category
#generates a table that compares the original feature counts to the counts in the final output
def makeTransferStats(originalGFF, newGFF):
    nameArray = ["CDS","exon","gene","mRNA","tRNA","ncRNA","rRNA","total features"]
    origCounts = countUniqueOccurences(originalGFF,nameArray)
    finalCounts = countUniqueOccurences(newGFF,nameArray)
    writeStatsToFile(origCounts,finalCounts, nameArray)

#generates a list of sets that contain unique IDs categorized by feature
def countUniqueOccurences(gff, names):
    sortedUniqueIDs = [set() for x in range(0,len(names))] #generates a list of sets
    with open(gff) as file:                                 #list index lines up with name array
        for line in file:
            values = line.split("\t")
            for i in range(0,len(names)-1):
                if (len(values) == 9 and values[2] == names[i]):
                    sortedUniqueIDs[i].add(getID(values[8])) #adds ID to set that corresponds to the feature
    total = 0
    for IDs in sortedUniqueIDs:
        total+=len(IDs) #counts up the total unique features
    for x in range(0,total):
        sortedUniqueIDs[len(sortedUniqueIDs)-1].add(x) #generates a set with unqiue values that equal the total
    return sortedUniqueIDs

#converts tab-delimited quast report into a 2d list
def parseQuast(results):
    row = []
    table = []
    with open(results) as file:
        for line in file:
            row = line.split("\t")
            if (len(row) == 3):
                row[2] = row[2][:len(row[2])-1]
            else:
                row[1] = row[1][:len(row[1])-1]
            table.append(row)
    return table
    
#writes out the unique feature count data as a comma-delimited file
def writeStatsToFile(oCounts,fCounts, names):
    # refQuastResults = parseQuast(ratt_dir+"/ref_quast/report.tsv")
    # queryQuastResults = parseQuast(ratt_dir+"/query_quast/report.tsv")
    with open(ratt_dir+"/transferStats.csv","w") as file:
        file.write("Feat.,Orig.,Final\n")
        for x in range(0,len(names)):
            file.write(names[x])
            file.write(",")
            file.write(str(len(oCounts[x])))
            file.write(",")
            file.write(str(len(fCounts[x])))
            file.write("\n")
        
        # file.write("Assembly Length,"+refQuastResults[15][1]+","+queryQuastResults[15][1]+"\n")

        # #if the given genomes have contigs within scaffolds
        # if len(refQuastResults[0]) == 3 and len(queryQuastResults[0]) == 3:    
        #     file.write("Scaffolds,"+refQuastResults[13][1]+","+queryQuastResults[13][1]+"\n")
        #     file.write("Contigs,"+refQuastResults[13][2]+","+queryQuastResults[13][2]+"\n")
        #     file.write("Ns per 100kbp,"+refQuastResults[20][2]+","+queryQuastResults[20][2]+"\n")
        # else:
        #     file.write("Contigs,"+refQuastResults[13][1]+","+queryQuastResults[13][1]+"\n")
        #     file.write("Ns per 100kbp,"+refQuastResults[20][1]+","+queryQuastResults[20][1]+"\n")

        # file.write("N50,"+refQuastResults[16][1]+","+queryQuastResults[16][1]+"\n")
        # file.write("L50,"+refQuastResults[18][1]+","+queryQuastResults[18][1]+"\n")

#outputs the transfer stats to stdout            
def printTransferStats():
    print "\n<RATTwithGFF.py> RATT TRANSFER STATISTICS"
    statsFile = ratt_dir+"/transferStats.csv"
    with open(statsFile) as stats:
        for line in stats:
            vals = line.split(",")
            spaces1 = (20-len(vals[0]))*' '
            spaces2 = (20-len(vals[1]))*' '
            output = vals[0]+spaces1+vals[1]+spaces2+vals[2]
            sys.stdout.write(output)

main()

