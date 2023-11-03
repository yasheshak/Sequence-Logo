import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.image as mplimg
import numpy as np
import argparse
import sys
import time

parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--sequenceFile', default="path/to/Splice_Sequences.fasta ", type=str, action='store', help='input sequences file')
parser.add_argument('-A', '--Afile', default="-A path/to/A.png", type=str, action='store', help='input A file')
parser.add_argument('-C', '--Cfile', default="-C path/to/C.png", type=str, action='store', help='input C file')
parser.add_argument('-G', '--Gfile', default="-G path/to/C.png", type=str, action='store', help='input G file')
parser.add_argument('-T', '--Tfile', default="-T path/to/C.png", type=str, action='store', help='input T file')
parser.add_argument('-o', '--outFile', default='Kothari_Yashesha_BME163_Assignment_Week5.png', type=str, action='store', help='output file for figure')

args = parser.parse_args()

sequenceFile = args.sequenceFile
Afile = args.Afile
Cfile = args.Cfile
Gfile = args.Gfile
Tfile = args.Tfile
outFile = args.outFile

print(outFile, sequenceFile, Afile, Cfile, Gfile, Tfile )


#Figure and panel dimensions and positioning
figureWidth=5
figureHeight=2

plt.style.use('BME163.mplstyle')

plt.figure(figsize=(figureWidth,figureHeight))

panelWidth=1.5
panelHeight=0.5
panel1 = plt.axes([0.5/figureWidth,0.6/figureHeight,panelWidth/figureWidth,panelHeight/figureHeight])
panel2 = plt.axes([2.2/figureWidth,0.6/figureHeight,panelWidth/figureWidth,panelHeight/figureHeight])
           #[x values], [y values]
for panel in [panel1,panel2]:
    panel.set_xlim(-10,10)
    panel.set_ylim(0,2)
    panel.set_xlabel(' Distance to \nSplice Site')
    panel.axvline(0, color='black', linestyle='-', linewidth = 0.5) 

panel2.set_yticks([])
panel2.set_title("3'SS")
panel1.set_title("5'SS")
panel1.set_ylabel('Bits')



#Parse through the fasta file, code from given fasta_reader.py from lecture
start=time.time()
sequence=''
seqDict={}
for line in open('Splice_Sequences.fasta'):
    if line.startswith('>'):
        if sequence:
            seqDict[name]=sequence
            sequence=''
        name=line[1:].strip().split()[0]
    else:
        sequence+=line.strip()

if sequence:
    seqDict[name]=sequence



#Separate the current dictionary into a 5' dictionary and a 3'
splice_5={}
splice_3={}

for key in seqDict:
    if key.startswith("5'"):
        splice_5[key] = seqDict[key]
    elif key.startswith("3'"):
        splice_3[key] = seqDict[key]

#Now make functions for each calculation and call the functions
#final function for base height
#loop through for positions, can do it with just one plotting line instead of 4 like in lecture nucs[nuc] method


#Function to calculate relative frequncy of the bases
def relative_frequency(splice_dict):
    rel_freq_dict = {}
    for i in range(20):
        nuc_count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

        for j in splice_dict.values():
            if i < len(j):
                nuc_count[j[i]] += 1

        rel_freq_dict[i] = {nuc: count / len(splice_dict) for nuc, count in nuc_count.items()}

    return rel_freq_dict

rel_freq_5 = relative_frequency(splice_5)
rel_freq_3 = relative_frequency(splice_3)


# Varibales to hold the calculation for En's
en_5 = (1 / np.log(2)) * ((4 - 1) / (2 * len(splice_5)))
en_3 = (1 / np.log(2)) * ((4 - 1) / (2 * len(splice_3)))

#Function to calculate the shannon entropy/uncertainty of the bases
import math
def shannon(splice_dict):
    entropy_5 = {}
    entropy_3 = {}

    rel_freq_dict = relative_frequency(splice_dict)

    for i in range(20):
        # shannon for 5'
        if i <= 20:
            entropy = 0
            for nuc in ['A', 'C', 'G', 'T']:
                rel_freq = rel_freq_dict[i].get(nuc, 0)
                if rel_freq > 0:
                    entropy -= rel_freq * math.log2(rel_freq)
            entropy_5[i] = entropy

        # shannon for 3'
        if i >= 0:
            entropy = 0
            for nuc in ['A', 'C', 'G', 'T']:
                rel_freq = rel_freq_dict[i].get(nuc, 0)
                if rel_freq > 0:
                    entropy -= rel_freq * math.log2(rel_freq)
            entropy_3[20 - i - 1] = entropy

    return entropy_5, entropy_3

entropy_5 = shannon(splice_5)
entropy_3 = shannon(splice_3)




#Function to calculate the information content, can use the separate dictionaries made in each one as paramters
def information_content(entropy_5, entropy_3, en_5, en_3):
    ic_5 = {}
    ic_3 = {}

    for i in range(20):
        # ic for 5'
        entropy = entropy_5.get(i, 0)
        ic_5[i] = np.log2(4) - (entropy + en_5)

        # ic for 3'
        entropy = entropy_3.get(i, 0)
        ic_3[i] = np.log2(4) - (entropy + en_3)

    return ic_5, ic_3

ic_5, ic_3 = information_content(entropy_5[0], entropy_3[0], en_5, en_3)


#Function that calculates the base height 
def calculate_baseheight(rel_freq, ic):
    heights = []
    for i, rel_freq_i in rel_freq.items():
        infocontent_i = ic[i]

        height = [(base, rel_freq_i[base] * infocontent_i) for base in rel_freq_i.keys()]
        height.sort(key=lambda x: x[1])
        heights.append(height)
        
    return heights

bh_5 = calculate_baseheight(rel_freq_5, ic_5)
bh_3 = calculate_baseheight(rel_freq_3, ic_3)



#PLOTTING
#inc_x for x axis positions, moves over 1 each time
#prev_height for stack, builds on top of previous
#####used nucs[nuc] method instead of plotting all 4 because doing it separately was confusing for me
for panel in [panel,panel2]:
    panel.set_xlim(-10,10)
    panel.set_ylim(0,2)

A=mplimg.imread('A.png')
T=mplimg.imread('T.png')
C=mplimg.imread('C.png')
G=mplimg.imread('G.png')



#Plotting for 5' panel
nucs = {'A': A, 'T': T, 'C': C, 'G': G}

inc_x = panelWidth * 13.5 / len(bh_5)

for i in range(len(bh_5)):
    prev_height = 0
    
    for j in range(len(bh_5[i])): # Always 4
        nuc = bh_5[i][j][0]
        height = prev_height + bh_5[i][j][1] / (panelHeight *2)
        panel1.imshow(nucs[nuc], extent = [i-10,i+inc_x-10,prev_height,height ], aspect='auto', origin='upper')
        prev_height = height


inc_x = panelWidth * 13.5 / len(bh_3)
for i in range(len(bh_3)):
    prev_height = 0
    
    for j in range(len(bh_3[i])): # Always 4
        nuc = bh_3[i][j][0]
        height = prev_height + bh_3[i][j][1] / (panelHeight *2)
        panel2.imshow(nucs[nuc], extent = [i-10,i+inc_x-10,prev_height,height ], aspect='auto', origin='upper')
        prev_height = height



plt.savefig(outFile, dpi=600)