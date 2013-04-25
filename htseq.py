# -*- coding: utf-8 -*-

# Author: Pedro Furió Tarí
# Data: 07/03/2013
#
# This script requires the following tools: Samtools & htseq-count

import getopt, sys, os.path
import threading
import Queue
import re

# Global variables
temporaryFiles = []
queueLock = threading.Lock()
workQueue = Queue.Queue(0)

class myThread (threading.Thread):
    
    def __init__(self, threadID, name, q):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.q = q
        
    def run(self):
        while not self.q.empty():
            data = self.q.get() # Get a job
            queueLock.acquire() # Acquire the lock
            print "Starting " + self.name
            queueLock.release() # Release the lock
            runHTSeq(data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],True)
            workQueue.task_done() # Let the queue know the job is finished

        queueLock.acquire() # Acquire the lock
        print "Exiting " + self.name
        queueLock.release() # Release the lock

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hm:s:a:t:i:nd:cb:g:p:", ["help", "mode=", "stranded=", "minaqual=", "type=", "idattr=", "sort", "destiny=", "clean_up", "bam=", "gtf=", "threads="])
    except getopt.GetoptError as err:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    mode = "union"
    stranded = "yes"
    minaqual = 0
    typed = "exon"
    idattr = "gene_id"
    sort = False
    destiny = "./"
    clean = True
    bam = None
    gtf = None
    threads = 1
    
    for o, a in opts:
        if o in ("-h","--help"):
            usage()
            sys.exit()
        elif o in ("-m", "--mode"):
            mode = a
        elif o in ("-s", "--stranded"):
            stranded = a
        elif o in ("-m", "--minaqual"):
            minaqual = int(a)
        elif o in ("-t", "--type"):
            typed = a
        elif o in ("-i", "--idattr"):
            idattr = a
        elif o in ("-n", "--sort"):
            sort = True
        elif o in ("-d", "--destiny"):
            if os.path.isdir(a):
                destiny = a
        elif o in ("-c", "--clean_up"):
            clean = False
        elif o in ("-b", "--bam"):
            if a.split(".")[-1] in ("bam","sam") and os.path.isfile(a):
                bam = a
        elif o in ("-g", "--gtf"):
            if os.path.isfile(a):
                gtf = a
        elif o in ("-p", "--threads"):
            threads = int(a)
        else:
            assert False, "Unhandled option"
    
    if destiny[-1] != "/":
        destiny = destiny + "/"
    
    if gtf != None and bam != None:
        run (mode, stranded, minaqual, typed, idattr, sort, destiny, clean, bam, gtf, threads)
    else:
        usage()

def usage():
    print "\nUsage: python htseq.py [options] <mandatory>"
    print "Options:"
    print "\t-h, --help:\n\t\t show this help message and exit"
    print "\t-m, --mode:\n\t\t mode to handle reads overlapping more than one feature(choices: union, intersection-strict, intersection-nonempty; default: union)"
    print "\t-s, --stranded:\n\t\t whether the data is from a strand-specific assay. Specify 'yes', 'no', or 'reverse' (default: yes). 'reverse' means 'yes' with reversed strand interpretation"
    print "\t-a, --minaqual:\n\t\t skip all reads with alignment quality lower than the given minimum value (default: 0)"
    print "\t-t, --type:\n\t\t feature type (3rd column in GFF file) to be used, all features of other type are ignored (default, suitable for Ensembl GTF files: exon)"
    print "\t-i, --idattr:\n\t\t GFF attribute to be used as feature ID (default, suitable for Ensembl GTF files: gene_id)"
    print "\t-o, --samout:\n\t\t write out all SAM alignment records into an output SAM file called SAMOUT, annotating each line with its feature assignment (as an optional field with tag 'XF')"
    print "\t-n, --sort:\n\t\t sort the bam file by name (necessary for paired-end reads"
    print "\t-d, --destiny:\n\t\t output directory (default, the directory of execution)"
    print "\t-c, --clean_up:\n\t\t Do not remove the intermediate files generated (default: remove intermediate files)"
    print "\t-p, --threads:\n\t\t Number of threads to run (default: 1)"
    print "Mandatory:"
    print "\t-b, --bam:\n\t\t bam/sam file to read"
    print "\t-g, --gtf:\n\t\t gtf file"
    print "\n07/03/2013. Pedro Furió Tarí.\n"

def bam2sam(bam, destiny):
    global temporaryFiles
    
    filename = bam.split("/")[-1].split(".")[-2]
    os.system("samtools view -h " + bam + " > " + destiny + filename + ".sam")

    temporaryFiles.append(destiny + filename + ".sam")
    
    return destiny + filename + ".sam"

def sortSam(sam, destiny):
    global temporaryFiles
    
    print "Sorting the BAM file..."
    filename = sam.split("/")[-1].split(".")[-2]
    os.system("samtools sort -n " + sam + " " + destiny + filename)
    
    temporaryFiles.append(destiny + filename + ".bam")
    
    print "BAM file sorted."
    return destiny + filename + ".bam"

def detectMaxMultiHits(sam, destiny):
    
    print "Retrieving the maximum number of multihits in the sam file..."
    
    # Create a file containing a list with the different multihits (uniq)
    os.system("awk -F'NH:i:' '{print $2}' " + sam + " | cut -f1 | sort -n | uniq | tail -1 > " + destiny + "tmp.tmp.maxNH")
    
    f = open(destiny + "tmp.tmp.maxNH", 'r')
    maxim = int(f.next())
    f.close()
    # Remove the temporal file generated
    os.system("rm " + destiny + "tmp.tmp.maxNH")
    
    print "The maximum number of multihits detected is " + str(maxim)
    
    return maxim  

def generateNHSam (sam, destiny, n_hit):
# This function will create a new sam file containing only those reads with n multihits
# It also will take only those which have pairs correctly paired

    global temporaryFiles
    
    print "Generating NH:i:" + str(n_hit) + " sam file..."   
    
    filename = sam.split("/")[-1]
    
#    filein = open(sam,'r')
#    output = open(destiny + filename + "." + str(n_hit), 'w')
#     
#    for linea in filein:
#        if linea[0] == "@":    # Write the headers
#            output.write(linea)
#        else:                  # Filter the mapped reads
#            l = linea.split()
#            cigar = int(l[1])
#            multi = l[12]
#            # We check whether they are paired end properly paired or single-end
#            if multi == "NH:i:" + str(n_hit) and (bin(cigar)[-2::] == "11" or bin(cigar)[-1::] == '0'):
#                l[12] = "NH:i:1"
#                for w in l:
#                    output.write(w + "\t")
#                output.write("\n")
#    output.close()
#    filein.close()
         
    # Put the header in the new file
    os.system("samtools view -SH " + sam + " > " + destiny + filename + "." + str(n_hit))
    
    # Now we generate the new file and change for NH:i:1
    #os.system("grep -v '^@' " + sam + " | grep -E $'NH:i:" + str(n_hit) + "\t' >> " + destiny + filename + "." + str(n_hit) + ".aux")
    #os.system("sed -i 's/NH:i:" + str(n_hit) + "/NH:i:1/g' " + destiny + filename + "." + str(n_hit) + ".aux")
    out = open(destiny + filename + "." + str(n_hit), 'a')
    infile = open(sam, 'r')

    for i in infile:
        if re.search("NH:i:" + str(n_hit) + "[\t\n]", i):
            nuevalinea = i.replace("NH:i:" + str(n_hit), "NH:i:1")

            l = nuevalinea.split()
            cigar = int(l[1])

            if bin(cigar)[-2::] == "11" or bin(cigar)[-1::] == '0':
                for w in l:
                    out.write(w + "\t")
                out.write("\n")

    out.close()
    infile.close()

    #filein  = open(destiny + filename + "." + str(n_hit) + ".aux", "r")
    #fileout = open(destiny + filename + "." + str(n_hit), "w")
    # If they are paired end -> Check that they are properly paired.
    # If single-end -> Perfect
    #for line in filein:
    #    if line[0] == "@":
    #        fileout.write(line)
    #    else:
    #        l = line.split()
    #        cigar = int(l[1])
            
    #        if bin(cigar)[-2::] == "11" or bin(cigar)[-1::] == '0':
    #            for w in l:
    #                fileout.write(w + "\t")
    #            fileout.write("\n")
    #fileout.close()
    #filein.close()
    
    # Remove the first file generated
    #os.system("rm " + destiny + filename + "." + str(n_hit) + ".aux")
    
    # Write the file that has been created in a file
    #temporalFiles(destiny + filename + "." + n_hit, destiny)
    temporaryFiles.append(destiny + filename + "." + str(n_hit))
    
    print "NH:i:" + str(n_hit) + " sam file generated"
    
    return destiny + filename + "." + str(n_hit)

def runHTSeq (mode, stranded, minaqual, typed, idattr, destiny, sam, gtf, n_hit, multihits):
    global temporaryFiles
    
    # Generate the Auxiliar sam file
    if multihits:
        newSam = generateNHSam(sam, destiny, n_hit)
        filename = newSam.split("/")[-1]
    else:
        newSam = sam
        filename = sam.split("/")[-1]
    
    print "Running HTSeq..."
    
    # Run HTSeq for that new sam file
    os.system("htseq-count -m " + mode + " -s " + stranded + " -a " + str(minaqual) + " -t " + typed + " -i " + idattr + " " + newSam + " " + gtf + " > " + destiny + filename + ".counts")

    if multihits:
        temporaryFiles.append(destiny + filename + ".counts")

        # Correct the file generated
        f = open(destiny + filename + ".counts", 'r')
        o = open(destiny + filename + ".counts2", 'w')

        for line in f:
            l = line.split()
            o.write(l[0] + "\t" + str(float(l[1])/n_hit) + "\n")

        f.close()
        o.close()

        # Move the last file generated to overwrite the one given by htseq
        os.system("mv " + destiny + filename + ".counts2" + " " + destiny + filename + ".counts")
    else:
        os.system("mv " + destiny + filename + ".counts" + " " + destiny + "htseq_multihit.counts")
    
    print "HTSeq has finished running."

def run(mode, stranded, minaqual, typed, idattr, sort, destiny, clean, bam, gtf, n_threads):
    
    # 1. Check whether is necessary to sort the sam file
    try:
        if sort and bam.split(".")[-1] == "bam":
            bam = sortSam(bam, destiny)
    except:
        print "ERROR: Could not sort the sam file"
        sys.exit(1)
    
    # 2. Check whether is a bam or a sam file
    try:
        if bam.split(".")[-1] == "bam":
            bam = bam2sam(bam, destiny)
    except:
        print "ERROR: The conversion bam to sam could not be properly done"
        sys.exit(1)

    # 3. Get the maximum number of multihits
    try:
        # If there's more than 30 multihits, we won't take them into account
        maxim = min(detectMaxMultiHits(bam, destiny), 30)
    except:
        os.system("rm " + destiny + "tmp.tmp.maxNH")
        print "WARNING: Could not detect the maximum number of hits in the mapping file " + bam
        print "HTSeq will be run as unihits"
        maxim = 1
        #exit()

    if maxim > 1:
        # 4. Run the different HTSeq functions with each of the multihits
        # Create the stack of jobs
        for i in range(1,maxim+1):
            argsVect = [mode, stranded, minaqual, typed, idattr, destiny, bam, gtf, i]
            workQueue.put(argsVect)

        # Create and lunch the threads
        for n_thread in range(n_threads):
            thread = myThread(n_thread + 1, "Thread-" + str(n_thread + 1), workQueue)
            thread.start()

        queueLock.acquire()  # Acquire the lock so we can print
        print "Waiting for threads to finish"
        queueLock.release()  # Release the lock

        workQueue.join()     # Wait till all the threads finish

        print "Starting to merge the counts into a single file..."

        # 5. Join the different files with counts
        conteos = {}
        # Initialize the dictionary of genes with the real counts of uniquely mapped reads
        file_counts1 = open(destiny + bam.split("/")[-1] + ".1.counts", 'r')
        print "Reading the file " + destiny + bam.split("/")[-1] + ".1.counts"
        for i in file_counts1:
            linea = i.split()
            conteos[linea[0]] = float(linea[1])
        file_counts1.close()

        # And add the multimapped reads
        for fileC in range(2,maxim+1):
            fileI = open(destiny + bam.split("/")[-1] + "." + str(fileC) + ".counts", 'r')
            print "Reading the file " + destiny + bam.split("/")[-1] + "." + str(fileC) + ".counts"
            for i in fileI:
                linea = i.split()
                conteos[linea[0]] = conteos[linea[0]] + float(linea[1])
            fileI.close()

        # Now, we write the results to a file
        file_counts = open(destiny + "htseq_multihit.counts", 'w')
        for i in conteos:
            file_counts.write(i + "\t" + str(conteos[i]) + "\n")
        file_counts.close()

    else:
        # No multihits
        # mode, stranded, minaqual, typed, idattr, sort, destiny, clean, bam, gtf, n_threads
        runHTSeq (mode, stranded, minaqual, typed, idattr, destiny, bam, gtf, 1,  False)

    # 6. If selected the option to remove all the intermediate files ...
    if clean:
        for i in temporaryFiles:
            os.system("rm " + i)
        print "Intermediate files removed."
    else:
        print "The following are the temporaryFiles that have been created: "
        print temporaryFiles

if __name__ == "__main__":
    main()
