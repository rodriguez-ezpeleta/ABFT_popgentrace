
# NAME

#         getSNPreads.py

# SYNOPSIS

#         getSNPreads.py -s selectedSNPs -c catalogName

# DESCRIPTION

#         From a stacks () catalog, and a list of snps, retrieves flanking region from genome

# PARAMETERS

#         -s
#                 selectedSNPs.
#         -c
#                 catalogName.
#         -g
#                 genomeName.

# AUTHOR

#         Naiara Rodriguez-Ezpeleta, nrodriguez@azti.es
#         December 2015



def  main():

        import string
        import sys
        import glob, os
        import re
        import gzip
        from Bio.Data.IUPACData import ambiguous_dna_values
        from Bio.Seq import Seq

# Test if the number of arguments is rigth

        if len(sys.argv) < 3:
                print __doc__
                sys.exit(1)

# Process arguments

        snpFile = ""
        catName = ""
        argv = sys.argv[1:]

        while len(argv) > 0 :
                arg = argv[0]
                if arg == "-s" :
                        if len(argv) > 1 :
                                snpFile = argv[1]
                        if os.path.exists(snpFile) == 0 :
                                print("\n\tERROR: SNP file \"" + snpFile + "\" could not be found!\n")
                                sys.exit(0)
                        argv = argv[2:]
                        continue

                if arg == "-c" :
                        if len(argv) > 1 :
                                catName = argv[1]
                        if os.path.exists("../" + catName + ".catalog.snps.tsv.gz") == 0 :
                                print("\n\tERROR: No SNP catalog file starting by \"" + catName + "\" could  be found!\n")
                                sys.exit(0)
                        if os.path.exists("../" + catName + ".catalog.tags.tsv.gz") == 0 :
                                print("\n\tERROR: No TAG catalog file starting by \"" + catName + "\" could  be found!\n")
                                sys.exit(0)
                        argv = argv[2:]
                        continue

                if arg == "-g" :
                        if len(argv) > 1 :
                                genomeFile = argv[1]
                        if os.path.exists(genomeFile) == 0 :
                                print("\n\tERROR: Genome file \"" + snpFile + "\" could not be found!\n")
                                sys.exit(0)
                        argv = argv[2:]
                        continue


        rev = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}
#       ambiguous_dna_values_inv = {v: k for k, v in ambiguous_dna_values.items()}
        ambiguous_dna_values_inv = {'A': 'A', 'C': 'C', 'GT': 'K', 'ACG': 'V', 'AG': 'R', 'G': 'G', 'CG': 'S', 'AC': 'M', 'CT': 'Y', 'T': 'T', 'ACT': 'H', 'AT': 'W', 'A
GT': 'D', 'CGT': 'B', 'ACGT': 'X'}

# read selected snp file and create a dictionary of tags/snps (seltags)

        snpLines = open(snpFile, "r").read().splitlines()
        seltags = {}
        for line in snpLines:
                if line.split("_")[0] not in seltags:
                        seltags[line.split("_")[0]] = []
                seltags[line.split("_")[0]].append(line.split("_")[1])

# read the catalog.snps.tsv file

        catSNP = "../" + catName + ".catalog.snps.tsv.gz"
        with gzip.open(catSNP, 'rb') as f:
                catsnpLines = f.read().splitlines()[1:]

# read the catalog.tags.tsv file
        catTAG = "../" + catName + ".catalog.tags.tsv.gz"
        with gzip.open(catTAG, 'rb') as f:
                cattagLines = f.read().splitlines()[1:]

# create a file to store the sequence of each tag
        tagsFasta = "selectedTags.fasta"
        tagsF = open(tagsFasta, "w")

# search the sequence for each tag in the catalog.tags.file; store in dictionary (tagSeq) and write in fasta file

        tagSeq = {}
        for line in cattagLines:
                col=line.split("\t")
                if col[2] in seltags:
                        tagSeq[col[2]] = col[9]
                        tagsF.write(">" + col[2] + "\n" + col[9] + "\n")
        tagsF.close()

# execute blast of tags against genome

        print "Executing blast.....\n"

        cmd = "blastn  -query selectedTags.fasta -db" + genomeFile + "-out selectedTags_blast.txt  -outfmt 6 -max_target_seqs 1"
        os.system(cmd)

        print "Blast finished.\n"

# parse blast result to get the contigs that matched to a tag; store in array blastContigs

        blastLines=open("selectedTags_blast.txt", "r").read().splitlines()

        blastContigs = {}
        contigsWtags = []
        for line in blastLines:
                contigsWtags.append(line.split("\t")[1])
                blastContigs[line.split("\t")[0]] = []
                blastContigs[line.split("\t")[0]].append(line.split("\t")[1])
                blastContigs[line.split("\t")[0]].append(line.split("\t")[8])
                blastContigs[line.split("\t")[0]].append(line.split("\t")[9])
                blastContigs[line.split("\t")[0]].append(line.split("\t")[6])
                blastContigs[line.split("\t")[0]].append(line.split("\t")[7])

# get sequence of contigs in blastsContigs; store in dictionary contigs

        genomeFasta=genomeFile + ".fa"
        genome=open(genomeFasta).read().splitlines()

        contigs = {}

        i=0

        while i < len(genome):
                if genome[i][0] == ">" and genome[i].split(" ")[0][1:] in contigsWtags:
                        contigs[genome[i].split(" ")[0][1:]] = genome[i+1]
                i = i + 2


        done = []
        for line in catsnpLines:
                col=line.split("\t")
                if col[2] in seltags:
                        bases=''.join(sorted(col[6:10])).replace("-", "")
                        pos = int(col[3])
                        if col[2] not in done:
                                x = list(tagSeq[col[2]])
                        if col[3] == seltags[col[2]][0]:
                                x[pos] = str("[" + bases[0] + "/" + bases[1] + "]")

                        else:
                                x[pos] = ambiguous_dna_values_inv[bases]

                        tagSeq[col[2]] = "".join(x)
                        done.append(col[2])

        snpList = "SNP_fluidigm.txt"
        snpL = open(snpList, "w")

        snpL.write("Tag\tPostion\tContig\tSense\tStart\tEnd\tSequence\n")
        for tag in blastContigs:
                t = list(tagSeq[tag])
                g = list(contigs[blastContigs[tag][0]])

                if int(blastContigs[tag][1]) <  int(blastContigs[tag][2]):
                        blastContigs[tag].append("F")

                        if (int(blastContigs[tag][1]) - int(blastContigs[tag][3])) - 1 < 250:
                                a = 0
                        else:
                                a = int(blastContigs[tag][1]) - int(blastContigs[tag][3]) - 250

                        if len(g) - int(blastContigs[tag][2]) + 90 - int(blastContigs[tag][4]) <  250:
                                d = len(g)

                        else:
                                d = int(blastContigs[tag][2]) + 90 - int(blastContigs[tag][4]) + 250

                        b = int(blastContigs[tag][1]) - int(blastContigs[tag][3])
                        c = int(blastContigs[tag][2]) + 90 - int(blastContigs[tag][4])

                        if b < 0 & c > len(g):
                                seq = t
                        elif b < 0 :
                                seq = t + g[c:d]
                        elif c > len(g):
                                seq = g[a:b] + t
                        else:
                                seq = g[a:b] + t + g[c:d]

                else:
                        blastContigs[tag].append("R")

                        if len(g) - (int(blastContigs[tag][1]) + int(blastContigs[tag][3]) - 1) < 250:
                                a = len(g)

                        else:
                                a = int(blastContigs[tag][1]) + int(blastContigs[tag][3]) - 1 + 250

                        if int(blastContigs[tag][2]) - int(blastContigs[tag][4]) + 90 - 1 <  250:
                                d = 0

                        else:
                                d = int(blastContigs[tag][2]) - int(blastContigs[tag][4]) - 90 - 1 - 250

                        c = int(blastContigs[tag][2]) + int(blastContigs[tag][4]) - 90 - 1
                        b = int(blastContigs[tag][1]) + int(blastContigs[tag][3]) - 1

                        first = list(Seq("".join(g[b:a])).reverse_complement())
                        second = list(Seq("".join(g[d:c])).reverse_complement())

                        if b < 0 & c > len(g):
                                seq = t
                        elif b < 0 :
                                seq = t + second
                        elif c > len(g):
                                seq = first + t
                        else:
                                seq = first + t + second

                tagSeq[tag] = "".join(seq)


                snpL.write(tag + "\t" + seltags[tag][0] + "\t" + blastContigs[tag][0] + "\t" + blastContigs[tag][5]  + "\t" + blastContigs[tag][1] + "\t" + blastContigs
[tag][2] + "\t" + tagSeq[tag] + "\n")

        snpL.close()
main()
