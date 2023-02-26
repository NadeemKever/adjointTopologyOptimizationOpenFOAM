import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

lines = []
with open('objectiveVal.txt') as fileObject:
    
    for line in fileObject:
        lines.append(float(line.strip()))

stdDev = np.std(lines)
mean = np.mean(lines)
print("mean = ",mean,"\nStandard Deviation: ",stdDev,"\nmax = ", max(lines),"\nmin = ", min(lines))
print("Converged Value: ", lines[-1])

cutOff = int(0.7 * len(lines))
tail = lines[cutOff:]
stdDev2 = np.std(tail)
tailMean = np.mean(tail)
print("Mean of last 30% of  Values: ", tailMean)
print("Tail standard deviation: ", stdDev2)
mpl.rc("font",size = 20)
plt.figure(figsize=(25, 14))
plt.plot(lines,"-*")
plt.xlabel("Iteration")
plt.ylabel("Objective Value")
plt.ylim((tailMean - (1/8)*stdDev, tailMean+(1/8)*stdDev))
# plt.ylim((tailMean - 10*stdDev2, tailMean + 10*stdDev2))
plt.savefig("Objective.png")
