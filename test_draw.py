import matplotlib.pyplot as plt
import numpy as np
import scipy as sp # type: ignore
import ast
import math
from scipy.integrate import quad # type: ignore
y=1
S=1
epsilon = math.radians(23.4)
t_b= (-172/365)*2*math.pi
t_e= (193/365)*2*math.pi

epsilon_o = 1000 # epsilon_o in W/m^2
w = 2*math.pi/(24*3600)  # angular velocity in rad/s
C = (S*epsilon_o/w)
o = 0

P_condition = math.atan(1/(math.tan(epsilon))) # P condition in radians

def intE_iMod2(a,b,c):
    if a%2 == 0:
        return quad(E_0, b, c)[0]
    else:
        return quad(E_1, b, c)[0]
result_array = []


for k in range(-89,89):
    print(k)
    result = 0
    y = math.radians(k)
    E_0 = lambda x: (1/3600000)*C*(math.cos(y)*(math.sin(math.acos(-math.tan(epsilon)*math.tan(y)*math.cos(x)))-math.sin(-math.acos(-math.tan(epsilon)*math.tan(y)*math.cos(x))))+math.sin(y)*math.cos(x)*math.sin(epsilon)*2*math.acos(-math.tan(epsilon)*math.tan(y)*math.cos(x)))
    E_1 = lambda x: (1/3600000)*C*math.sin(y)*math.sin(epsilon)*math.cos(x)*(2*math.pi)
    
    if (0 <= y < P_condition) or (-P_condition < y <= 0):
        result = (365/(2*math.pi))*intE_iMod2(0, t_b, t_e)
        print(result)
        result_array.append(result)
    else:
        o = 0
        if y < 0: 
            o = 1
        if y > 0:
            t_1= -math.acos(-(1/(math.tan(y)*math.tan(epsilon))))
            t_2= -math.acos((1/(math.tan(y)*math.tan(epsilon))))
            t_3= math.acos((1/(math.tan(y)*math.tan(epsilon))))
            t_4= math.acos(-(1/(math.tan(y)*math.tan(epsilon))))
        else:
            t_1= -math.acos((1/(math.tan(y)*math.tan(epsilon))))
            t_2= -math.acos(-(1/(math.tan(y)*math.tan(epsilon))))
            t_3= math.acos(-(1/(math.tan(y)*math.tan(epsilon))))
            t_4= math.acos((1/(math.tan(y)*math.tan(epsilon))))
        R_n = np.sort([t_b, t_1, t_2, t_3, t_4, t_e])
        R_ndup = np.unique(R_n)
        keep_mask = []
        for val in R_ndup:
            if y >= 0 :
                if (val > t_b and val < t_e) and (val <= t_4 or val >= t_1):
                    keep_mask.append(True)
                else:
                    keep_mask.append(False)
            else:
        
                if (val >= t_b and val <= t_e) and (val <= t_2 or val >= t_3):
                    keep_mask.append(True)
                else:
                    keep_mask.append(False)

        R_n_filtered = R_ndup[keep_mask]
        
        for i in range(len(R_n_filtered)-1):
            if y < 0 and R_n_filtered[i] == t_2:
                pass
            else :  
                result += (365/(2*math.pi))*intE_iMod2(i+o, R_n_filtered[i], R_n_filtered[i+1])
                print("this is a calculation :", intE_iMod2(i+o, R_n_filtered[i], R_n_filtered[i+1]))
        print(result)
        print(y)
        result_array.append(result)


x_values = np.linspace(-90, 90, len(result_array))
tikz = []
for i in range(178):
    tikz.append([i-89,result_array[i]])
data = np.asarray(tikz)

# Using numpy.savetxt() to save the array 'data' into a CSV file named "test.csv"
# The 'delimiter=","' argument specifies the delimiter to use in the CSV file as a comma ","
np.savetxt("test.csv", data, delimiter=",") 
plt.plot(x_values, result_array, label='Result Array')
plt.title('Annual Solar Output based on Latitude')
plt.xlabel('Latitude (°)')
plt.ylabel('Solar Output (kWh)')
plt.xticks(np.arange(-90, 91, 30))
plt.grid(True)
plt.legend()
plt.savefig('test_plot.png')
plt.show()