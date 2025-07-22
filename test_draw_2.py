import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scipy as sp # type: ignore
import ast
import math
import pandas as pd
import csv
from scipy.integrate import quad # type: ignore
mpl.rcParams.update({
    "font.family": "serif",      
    "font.size": 20,             
    "axes.titlesize": 20,        
    "axes.labelsize": 20,        
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
})
y=1
S=1
epsilon = math.radians(23.4)
t_b= (-172/365)*2*math.pi
t_e= (193/365)*2*math.pi
epsilon_o = 1000 # epsilon_o in W/m^2
w = 2*math.pi/(24*3600)  # angular velocity in rad/s
C = (S*epsilon_o/w)
P_condition = math.atan(1/(math.tan(epsilon))) # P condition in radians
array = []
full_array = []
for n in range(0, 365):
    theta_s = (n-172)*2*math.pi/365
    array = []
    for l in range(-89,89):
        result = 0
        E_0 = lambda x: C*(math.cos(math.radians(l))*(math.sin(math.acos(-math.tan(epsilon)*math.tan(math.radians(l))*math.cos(x)))-math.sin(-math.acos(-math.tan(epsilon)*math.tan(math.radians(l))*math.cos(x))))+math.sin(math.radians(l))*math.cos(x)*math.sin(epsilon)*2*math.acos(-math.tan(epsilon)*math.tan(math.radians(l))*math.cos(x)))
        E_1 = lambda x: C*math.sin(math.radians(l))*math.sin(epsilon)*math.cos(x)*(2*math.pi)

        if ((0 <= math.radians(l) < P_condition) or (-P_condition < math.radians(l) <= 0)) == False:
            if math.radians(l) >= 0:
                t_1= -math.acos(-(1/(math.tan(math.radians(l))*math.tan(epsilon))))
                t_2= -math.acos((1/(math.tan(math.radians(l))*math.tan(epsilon))))
                t_3= math.acos((1/(math.tan(math.radians(l))*math.tan(epsilon))))
                t_4= math.acos(-(1/(math.tan(math.radians(l))*math.tan(epsilon))))
            else:
                t_1= -math.acos((1/(math.tan(math.radians(l))*math.tan(epsilon))))
                t_2= -math.acos(-(1/(math.tan(math.radians(l))*math.tan(epsilon))))
                t_3= math.acos(-(1/(math.tan(math.radians(l))*math.tan(epsilon))))
                t_4= math.acos((1/(math.tan(math.radians(l))*math.tan(epsilon))))
        if ((0 <= math.radians(l) < P_condition) or (-P_condition < math.radians(l) <= 0)):
            result = E_0(theta_s)
        elif (( t_1< theta_s < t_2 ) or (t_3 < theta_s < t_4)):
            result = E_0(theta_s)
        elif (l>0 and (t_2 < theta_s < t_3)) or (l<0 and (t_1 < theta_s < t_4) == False):
            result = E_1(theta_s)
        array.append(result)
    full_array.append(array)
calc = 0
for i in range(0, 365):         # Time or day (y-axis)
    calc += full_array[i][89]
print(calc)
csv_array = []

for i in range(0, 365):         # Time or day (y-axis)
    for j in range(0, 178):     # Latitude or longitude index (x-axis)
        csv_array.append([j - 89, i, round(full_array[i][j])])  # [x, y, z]

data = np.asarray(csv_array)
print(np.max(csv_array))
# Save CSV with header

with open("test_2.csv", mode="w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["x", "y", "z"])  # HEADER!
    writer.writerows(data)

# Read and pivot for plotting
df = pd.read_csv("test_2.csv")
pivot = df.pivot(index='y', columns='x', values='z')

plt.figure(figsize=(7, 7))
im = plt.imshow(pivot, origin='lower', aspect='auto', cmap='viridis', extent=[-90, 90, 0, 364])


cbar = plt.colorbar(im)
cbar.set_label('Solar Output (kWh)')


cbar.set_ticks([0, 6861416, 13722832, 20584248, 27445664, 34307081])
cbar.set_ticklabels(["0", "2.0", '4.0', '6.0',"8.0", "10.0"])

plt.xlabel('Latitude (°)')
plt.ylabel('Day of Year (n)')
plt.xticks(ticks=[-90, -60, -30, 0, 30, 60, 90])
plt.yticks(ticks=[0, 91, 182, 273, 365])
plt.title('Figure 10: Heatmap of daily solar output')
plt.tight_layout()
plt.show()