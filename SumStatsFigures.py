import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import Aim1

cosmic_summary = Aim1.readCSV("C:/Users/eelias13/Desktop/KinaseModulation/cosmic_summary_stats.csv")
gnomad_summary = Aim1.readCSV("C:/Users/eelias13/Desktop/KinaseModulation/gnomad_summary_stats.csv")

# set width of bar
barWidth = 0.25
fig = plt.subplots(figsize=(12, 8))

# set height of bar
COSMIC = [cosmic_summary.iloc[0][1], cosmic_summary.iloc[0][2], cosmic_summary.iloc[0][3], cosmic_summary.iloc[0][4], cosmic_summary.iloc[0][5], cosmic_summary.iloc[0][6], cosmic_summary.iloc[0][7]]
GNOMAD = [gnomad_summary.iloc[0][1], gnomad_summary.iloc[0][2], gnomad_summary.iloc[0][3], gnomad_summary.iloc[0][4], gnomad_summary.iloc[0][5], gnomad_summary.iloc[0][6], gnomad_summary.iloc[0][7]]

# Set position of bar on X axis
br1 = np.arange(len(COSMIC))
br2 = [x + barWidth for x in br1]

# Make the plot
plt.bar(br1, COSMIC, color='r', width=barWidth, edgecolor='grey', label='COSMIC')
plt.bar(br2, GNOMAD, color='g', width=barWidth, edgecolor='grey', label='GNOMAD')

# Adding Xticks
plt.xlabel('Mutation Type', fontweight='bold', fontsize=15)
plt.ylabel('Mutation Count', fontweight='bold', fontsize=15)
plt.xticks([r + barWidth for r in range(len(COSMIC))], ['Total', 'Total Phosphosite', 'Null', 'Phosophomimetic', 'ST to Y', 'Y to ST', 'Flanking Region'])

plt.legend()
plt.show()