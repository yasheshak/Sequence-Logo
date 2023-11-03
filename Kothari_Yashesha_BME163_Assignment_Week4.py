import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import argparse
import random

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--identityFile', default="BME163_Input_Data_4.ident", type=str, action='store', help='input identity file')
parser.add_argument('-c', '--coverageFile', default="BME163_Input_Data_4.cov", type=str, action='store', help='input coverage file')
parser.add_argument('-o', '--outFile', default='Kothari_Yashesha_BME163_Assignment_Week4.png', type=str, action='store', help='output file for figure')

args = parser.parse_args()

identityFile = args.identityFile
coverageFile = args.coverageFile
outFile = args.outFile

print(outFile, identityFile, coverageFile)

iBlue=(44/255,86/255,134/255)
iOrange=(230/255,87/255,43/255)
iYellow=(248/255,174/255,51/255)
iGreen=(32/255,100/255,113/255)
colors = [iBlue, iGreen, iYellow, iOrange]

figureWidth=5
figureHeight=5

plt.style.use('BME163')

plt.figure(figsize=(figureWidth,figureHeight))

panelWidth=3
panelHeight=4

panel = plt.axes([0.5/figureWidth,0.75/figureHeight,panelWidth/figureWidth,panelHeight/figureHeight])
           #[x values], [y values]
panel.set_xticks(range(75,105,5))
panel.set_yticks([1,2,3,4])
panel.set_yticklabels(["1-3", "4-6", "7-9", ">=10"])
panel.set_xlim(75, 100)
panel.set_ylim(0.5, 4.5)


plt.xlabel('Identity (%)')
plt.ylabel('Subread Coverage')

# parse data
cov_dict = {}
with open("BME163_Input_Data_4.ident", "r") as covFile:
    for line in covFile:
            name, identitypercent = line.split()
            cov_dict[name] = float(identitypercent)

List1, List2, List3, List4 = [], [], [], []

with open("BME163_Input_Data_4.cov", "r") as identFile:
    for line in identFile:
        name, coverage = line.split() 
        if name in cov_dict.keys():
            if 1<=float(coverage)<=3 and len(List1) < 1000:
                List1.append(cov_dict[name])
            elif 4<=float(coverage)<=6 and len(List2) < 1000:
                List2.append(cov_dict[name])
            elif 7<=float(coverage)<=9 and len(List3) < 1000:
                List3.append(cov_dict[name])
            elif float(coverage)>=10 and len(List4) < 1000:
                List4.append(cov_dict[name])
        if len(List1) == 1000 and len(List2) == 1000 and len(List3) == 1000 and len(List4) == 1000:
            break

# swarm plot function - need to do this 4 times with different xValue list (list1, List2, etc)

markersize = 1
plottedPoints = []
notPlaced = []
minDist = markersize/30
mainWidth = 10
mainHeight = 8


for i, L in enumerate([List1, List2, List3, List4]):
    #print(i)
    ymin= i + 1
    ymax= i + 6
    xmin= 75
    xmax= 100
    
    yrange= ymax-ymin
    xrange = xmax-xmin
    
    shift = ((minDist/20)*xrange)/mainWidth
    
    plottedPoints = []
    notPlottedCount = 0
    notPlottedPos = []

    for x1 in L:
        yPos = ymin
        plotted=False
        multiplier = [1, -1]
      
        if len(plottedPoints)==0:
            plottedPoints.append((x1,yPos))
        
        else:
            for move in np.arange(0, 0.4, shift):
                for mult in multiplier:
                    y1 = yPos + move*mult
                    distList=[]
                    for coords2 in plottedPoints:
                        x2,y2 = coords2[0],coords2[1]
                        xdist=(np.abs(x1-x2)/xrange)*mainWidth
                        ydist=(np.abs(y1-y2)/yrange)*mainHeight
                        distance=(xdist**2+ydist**2)**0.5
                        distList.append(distance)
                    if min(distList)>minDist:
                        plotted=True
                        break
                if plotted:
                    break
            if plotted:
                plottedPoints.append((x1,y1))
            else:
                break

   # print("plotting", len(plottedPoints))
    #print(plottedPoints[10])
    for coords in plottedPoints:
        x,y = coords[0],coords[1]
        panel.plot(x,y,marker='o',mew=0, mfc=colors[i],ms=1)

    if notPlottedCount > 0:
        print(f"{notPlottedCount} points could not be plotted:")
        for pos in notPlottedPos:
            print(f"position {pos}")


# median1 = np.median(List1)
# panel.plot((median1, median1), (min(List1), max(List1)), color = 'red', linewidth = 0.75, ms = 0, mew = 0)
#panel.plot((median2, median2), (yPos of bottom of line, yPos of top of line), color = 'red', linewidth = 0.75, ms = 0, mew = 0)

median1 = np.median(List1)
median2 = np.median(List2)
median3 = np.median(List3)
median4 = np.median(List4)

panel.plot((median1, median1), (0.6, 1.4), color='red', linewidth=0.75, ms=0, mew=0)
panel.plot((median2, median2), (1.6, 2.4), color='red', linewidth=0.75, ms=0, mew=0)
panel.plot((median3, median3), (2.6, 3.4), color='red', linewidth=0.75, ms=0, mew=0)
panel.plot((median4, median4), (3.6,4.4), color='red', linewidth=0.75, ms=0, mew=0)

# # plt.savefig('pre.png', dpi=600)

panel.tick_params(bottom=True, labelbottom=True,\
                   left=True, labelleft=True, \
                   right=False, labelright=False,\
                   top=False, labeltop=False)

plt.savefig(outFile, dpi=600)



