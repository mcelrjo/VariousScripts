'''This script will change the fasta names of trinity scripts prior to processing 
with fasta_insert_query_creator.py.  Query creator will then pull names and seqeunces from the 
fasta file to create a text file that can then be used in the SQL database.
'''
#Start the test block to confirm regex is functional
import re, sys
import pandas

fasta = sys.argv[1]
trinotate = sys.argv[2]
output = sys.argv[3]

def fastaIdentifierRenamer(fastaFile, trinotateFile, outputFile):
    '''To rename the fasta identifiers, open the fasta file and find an identifier line,
    next loop through the trinotateFile until a match is found, then rename the identifier.
    Trinotate file was used to create a list oftuples containing identifierName and geneName
    '''
    
    trinotate = pandas.read_table(trinotateFile)
    
    geneDescriptionRegex = '\S+ Full=([\S\s]*?);'

    geneNameList = []
    tranIDList = [] 
    
    for anot in trinotate['sprot_Top_BLASTX_hit']:
        if anot == '.':
            geneNameList.append(anot)
        else:
            geneName = re.search(geneDescriptionRegex, anot)
            gene = geneName.groups()
            geneNameList.append(gene[0])
        
    for anot in trinotate['transcript_id']:
        tranIDList.append(anot)
    
    combinedNames = []    
    for i in range(len(geneNameList)):
        combinedNames.append((tranIDList[i], geneNameList[i]))
        
    fileToWrite = open(outputFile, 'w')
    
    trinityRegex2 = '>(\S*) (.*)'

    fasta = open(fastaFile)
        
    for line in fasta:
        if line.startswith('>'):
            trinitysearch2 = re.search(trinityRegex2, line)
            trinityNamesDivided = trinitysearch2.groups()
            
            for name in combinedNames:
                if name[0] == trinityNamesDivided[0]:
                    if name[1] == '.':
                        fileToWrite.write(">"+trinityNamesDivided[1]+" | "+trinityNamesDivided[0]+ " | unidentified transcript \n")
                        break
                    else:
                        fileToWrite.write(">"+trinityNamesDivided[1]+" | "+trinityNamesDivided[0]+ " | "+name[1]+'\n')
                        break
                    #if it is just a period then need to keep original or just skip the whole thing
                    
                    
        else:
            fileToWrite.write(line)
            #If the line is a sequence it write the sequence to the file
    
    fileToWrite.close()
    fasta.close()
    
        
fastaIdentifierRenamer(fasta, trinotate, output)

def addNewGeneDescriptionWithPipe(startFile, outputFile, wrongName = '---NA---', correctName = 'unidentified transcript'):
    '''Takes a fasta file that has an ID and gene description separated by a pipe 
    and looks for a wrong gene description (wrongName) and replaces it with the correct name 
    (correctName).  THE IDENTIFIER MUST CONTAIN A PIPE IN THIS INSTANCE.
    '''
    
    regex = r'>(\S+)\|([\S\s]+)'
    
    fileToWrite = open(outputFile, 'w')
    fasta = open(startFile)
        
    for line in fasta:
        if line.startswith('>'):
            identifierSearch = re.search(regex, line.strip())
            identifierGroups = identifierSearch.groups()
            
            if identifierGroups[1] == wrongName:
                fileToWrite.write(">" + identifierGroups[0] +" | "+ correctName + '\n')
                
            else:
                fileToWrite.write(">"+ identifierGroups[0] +" | "+ identifierGroups[1] + '\n')
                    
                    
        else:
            fileToWrite.write(line)
    
    fileToWrite.close()
    fasta.close()
    

#addNewGeneDescriptionWithPipe(r'C:\Users\jsm0010\Dropbox\WeedGenomics\databases from old site\Eleusine_indica_leaf_evigene_main.fasta', r'C:\Users\jsm0010\Dropbox\WeedGenomics\databases from old site\EleusineIndicaFinal.fasta')

def addNewGeneDescriptionWithoutPipe(startFile, outputFile, wrongName = '---NA---', correctName = 'unidentified transcript'):
    '''Takes a fasta file that has an ID and gene description separated by a pipe 
    and looks for a wrong gene description (wrongName) and replaces it with the correct name 
    (correctName).  THE IDENTIFIER MUST CONTAIN A PIPE IN THIS INSTANCE.
    '''
    
    regex = r'>(\S+)\s([\S\s]+)'
    
    fileToWrite = open(outputFile, 'w')
    fasta = open(startFile)
        
    for line in fasta:
        if line.startswith('>'):
            identifierSearch = re.search(regex, line.strip())
            identifierGroups = identifierSearch.groups()
            
            if identifierGroups[1] == wrongName:
                fileToWrite.write(">" + identifierGroups[0] +" | "+ correctName + '\n')
                
            else:
                fileToWrite.write(">"+ identifierGroups[0] +" | "+ identifierGroups[1] + '\n')
                    
                    
        else:
            fileToWrite.write(line)
    
    fileToWrite.close()
    fasta.close()
