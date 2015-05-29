#This python script will take the output from ebot obtain summaries and extract the title line and the gi line
from sys import argv
import re



ebotFile = open(argv[1], 'r')

currentLine = ""
gi = ""
title = ""
journal = ""
titleTrue = False
titleCounter = 0
JournalCounter = 0
currentLine = ""
firstEntry = True
sequence = False
fasta = ""


for line in ebotFile:
    if line.startswith("VERSION") and titleTrue == False:
        if firstEntry == False and journal != "Unpublished" and title != "Direct Submission ":
		#print "a." + fasta + ".b"
		print gi + "\t" + title + "\t" + journal + "\t" + fasta
		#print "\t" + "fail"
		#print "\n" + journal
		title = ""
		currentLine = ""
		fasta = ""
		#exit()
	else:
	    title = ""
	    currentLine = ""
	    fasta = ""
		#exit()
        line = line.split()
        gi = line[2].replace("GI:","")
        titleCounter = 0
	JournalCounter = 0
	firstEntry = False

    elif line.startswith("  TITLE") and titleCounter == 0:
        title = line
        title = title.replace("  TITLE     ","")
        title = title.replace("\n"," ")
	title = title.strip("\n")
        titleTrue = True
        titleCounter +=1
	firstEntry = False

    

    elif line.find("  JOURNAL") != -1 and titleTrue == True:
        titleTrue = False
	journal = line
	journal = journal.replace("  JOURNAL   ", "")
	journal = journal.strip("\n")
	#print "a." + journal + ".b"


    elif titleTrue == True:
        currentLine = line.replace("\n","")
        title += currentLine.replace("            ","")
        title += " "
	title = title.strip("\n")
	
    elif line.startswith("ORIGIN"):
	sequence = True
    elif sequence == True and line.find("//") == -1:
	splitLine = line.split(" ")
        for element in splitLine:
            if re.match("\D", element):
		#print element
                fasta += element.strip("\t\n")
    elif line.find("//") != -1:
	sequence = False



            
    
#print gi + "\t" + title