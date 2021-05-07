"""For this assignment, you will be reading in this file of apple genes and, based on these coding sequences, generated a codon usage bias table for this species.  If you are unfamiliar with this concept, the wiki page at the top will help get you up to speed.

Some things to keep in mind:

This is the first time you'll be using a full dataset.  Expect your numbers to get big
You will probably want to use the dictionary of amino acids and codons from Assignment 2 for this
We're looking for the frequency of the occurrence of a codon relative to other codons of the same amino acid.  As such, you will have multiple counts to keep track of.  Think about how you want to keep track of these counts, and what combinations of dictionaries, lists, and tuples are best suited to this task. Plan this BEFORE you start.
Just as important as it is to get these counts and generate a result, you must also consider how best to present these results. Your code should generate a human-readable file that shows your frequencies in a way that is verbose and makes sense.  Your success in presenting the results will be considered just as important as how you arrived to them.
"""

codons_dict = {'TTT' : 'Phe','TTC' : 'Phe', 'TTA' : 'Leu', 'TTG' : 'Leu', 'TCT' : 'Ser', 'TCC' : 'Ser', 'TCA' : 'Ser', 'TCG' : 'Ser',
	'TAT' : 'Tyr', 'TAC' : 'Tyr', 'TAA' : 'STOP', 'TAG' : 'STOP', 'TGT' : 'Cys', 'TGC' : 'Cys', 'TGA' : 'STOP', 'TGG' : 'Trp', 'CTT' : 'Leu',
	'CTC' : 'Leu', 'CTA' : 'Leu', 'CTG' : 'Leu', 'CCT' : 'Pro', 'CCC' : 'Pro', 'CCA' : 'Pro', 'CCA' : 'Pro', 'CCG' : 'Pro', 'CAT' : 'His', 'CAC' : 'His',
	'CAA' : 'Gln', 'CAG' : 'Gln', 'CGT' : 'Arg', 'CGC' : 'Arg', 'CGA' : 'Arg', 'CGG' : 'Arg', 'ATT' : 'Ile', 'ATC' : 'Ile', 'ATA' : 'Ile', 'ATG' : 'Met',
	'ACT' : 'Thr', 'ACC' : 'Thr', 'ACA' : 'Thr', 'ACG' : 'Thr', 'AAT' : 'Asn', 'AAC' : 'Asn', 'AAA' : 'Lys', 'AAG' : 'Lys', 'AGT' : 'Ser', 'AGC' : 'Ser',
	'AGA' : 'Arg', 'AGG' : 'Arg', 'GTT' : 'Val', 'GTC' : 'Val', 'GTA' : 'Val', 'GTG' : 'Val', 'GCT' : 'Ala', 'GCC' : 'Ala', 'GCA' : 'Ala', 'GCG' : 'Ala',
	'GAT' : 'Asp', 'GAC' : 'Asp', 'GAA' : 'Glu', 'GAG' : 'Glu', 'GGT' : 'Gly', 'GGC' : 'Gly', 'GGA' : 'Gly', 'GGG' : 'Gly'}

#aa_dict = {'Met':['ATG'], 'Phe':['TTT', 'TTC'], 'Leu':['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'Cys':['TGT', 'TGC'], 'Tyr':['TAC', 'TAT'], 'Trp':['TGG'], 
#	'Pro':['CCT', 'CCC', 'CCA', 'CCG'], 'His':['CAT', 'CAC'], 'Gln':['CAA', 'CAG'], 'Arg':['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'Ile':['ATT', 'ATC', 'ATA'],
#	'Thr':['ACT', 'ACC', 'ACA', 'ACG'], 'Asn':['AAT', 'AAC'], 'Lys':['AAA', 'AAG'], 'Ser':['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG'], 'Val':['GTT', 'GTC', 'GTA', 'GTG'],
#	'Ala':['GCT', 'GCC', 'GCA', 'GCG'], 'Asp':['GAT', 'GAC'], 'Glu':['GAA', 'GAG'], 'Gly':['GGT', 'GGC', 'GGA', 'GGG'], '*':['TAA','TAG','TGA']}	#"""Borrowing this, thanks for providing it"""

bias_codons_counter = {}

bias_amino_acids_counter = {}

#fh = open("/mnt/c/UNC_Fall_2019/BINF_6111/Lab_10_(apple_genome)/Mdomestica_SUBSET_testing.fa")	#I am testing my code on a smaller scale first. 
fh = open("/mnt/c/UNC_Fall_2019/BINF_6111/Lab_10_(apple_genome)/Mdomestica_491_v1.1.cds_primaryTranscriptOnly.fa")

d = {} 
header = ""
temp_empty_list = []

for line in fh:
	line = line.rstrip("\n")
	
	if (line.startswith(">")) and (line != header) and (header == ""):	#This will be use the first time a header is meet.
		header = line
		continue	
	
	if (line.startswith(">")) and (line != header) and (header != ""):	#This will be used for all headers after the first.
		a = "".join(temp_empty_list)
		
		d[header] = a
		header = line
		temp_empty_list = []
		a = None
		continue
	if (not line.startswith(">")) and (header != ""):
			temp_empty_list.append(line)
			
a = "".join(temp_empty_list)	#These two lines
d[header] = a					#are to add the last pair to the dictionary 

k = d.keys()
#print(len(k))

#for key in k:
#	print(key)

#print("\n\n")
#for key in k:
#	print(key,d[key] + "\n\n")
	
	
	
for h in k:
	value = d[h]	#The sequence
	for i in range(0,len(value)):
		if len(value[i:i+3]) == 3:	#len(value[(3*i):(3*i)+3])
			CODON = value[i:i+3]
			if "N" in CODON:
				continue
			else:
				bias_codons_counter[CODON] = bias_codons_counter.get(CODON,0) + 1
		else:
			continue
			

for i in bias_codons_counter.keys():
	print(i,bias_codons_counter[i])
	
for codon in bias_codons_counter.keys():
	
	if bias_amino_acids_counter.get(codons_dict[codon],None) == None:
		bias_amino_acids_counter[codons_dict[codon]] = [codon + " " + str(bias_codons_counter[codon])]
	else:
		bias_amino_acids_counter[codons_dict[codon]].append(codon + " " + str(bias_codons_counter[codon]))
		
for i in bias_amino_acids_counter.keys():
	print(i,bias_amino_acids_counter[i])
	
all = []
summary = ""
stats = bias_amino_acids_counter.keys()
for unit in stats:
	summary = summary + unit + ":" 
	values = bias_amino_acids_counter[unit]
	for x in values:
		y = x.rsplit(" ")
		#print(y)
		z = y[1]
		#print(z)
		z = int(z)
		all.append(z)
	aa_type_count = sum(all)
	for xx in values:
		yy = xx.rsplit(" ")
		zz = yy[1]
		zz = int(zz)
		summary = summary + str(yy[0]) + " has a frequency of " + str(round((zz/aa_type_count),4))  
		summary = summary + "\t"
	summary = summary + "\n"
	all = []
print(summary)

"""This format works"""
nf = open("FULL amino acid frequency output", 'w')
nf.write(summary)
nf.close()
  