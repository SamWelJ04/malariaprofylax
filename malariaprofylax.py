import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.linalg import norm
from scipy.integrate import quad

data = pd.read_excel("pk.xlsx")
time_points = np.array([0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 18.0 , 24.0, 48.0, 96.0])

def f(t, p):
    return p[0] * p[1]* p[2]/(p[1]-p[4])*(np.exp(-p[4]*t)-np.exp(-p[1]*t)) + p[0] * p[1] * p[3]/(p[1]-p[5])*(np.exp(-p[5]*t)-np.exp(-p[1]*t))
    
def Q(p, t, y) : # Omkastad ordning på p , t , y
    return norm(y - f(t, p))

F = 0.6
k_a = 0.9
A = 3
B = 4
l = 0.2 
mu = 0.1

p_start = np.array([F, k_a, A, B, l, mu])

def C(t):
    return p_opt[0] * p_opt[1]* p_opt[2]/(p_opt[1]-p_opt[4])*(np.exp(-p_opt[4]*t)-np.exp(-p_opt[1]*t)) + p_opt[0] * p_opt[1] * p_opt[3]/(p_opt[1]-p_opt[5])*(np.exp(-p_opt[5]*t)-np.exp(-p_opt[1]*t))

def tC(t):
    return t*C(t)

k_el_lst = []
AUC_lst = []
CL_lst = []
MRT_lst = []
V_ss_lst = []

xVar = np.linspace(0,100,1000)
for i in range(10):
    y = [f(t, p_start) for t in time_points]
    y_data = data["Conc"][10*i:10*(i+1)]
    p_opt = minimize(Q, p_start, (time_points , y_data)).x
    y_opt = [f(t, p_opt) for t in xVar]
    xVar = list(xVar)
    maxVal = max(y_opt)
    for j in range(y_opt.index(maxVal), len(y_opt)):
        if y_opt[j] <= maxVal/2:
            half_time = y_opt[j]
            break
    
    k_el = np.log(2)/half_time
    AUC, _ = quad(C, 0, 10000)
    dos = 150
    CL = dos/AUC
    g, _ = quad(tC, 0 , 10000)
    MRT = g/AUC
    V_ss = CL * MRT

    k_el_lst.append(k_el)
    AUC_lst.append(AUC)
    CL_lst.append(CL)
    MRT_lst.append(MRT)
    V_ss_lst.append(V_ss)
    
    #if i!=9:
        #plt.plot(xVar, y_opt, linestyle='dashed',color="lightblue")
    #else:
        #plt.plot(xVar, y_opt, linestyle='dashed',color="lightblue", label="Individuella dos-responskurvor")

y_data_all = data["Conc"]
time_points_all = data["Time"]
p_opt = minimize(Q, p_start, (time_points_all , y_data_all)).x

y_opt = [f(t, p_opt) for t in xVar]
#plt.plot(xVar, y_opt, linewidth=3, color='maroon', label="Skattad dos-responskurva") 
#plt.legend()
#plt.xlabel('Tid(h)')
#plt.ylabel('Koncentration(mg/L)')
#plt.show()

y_bad = data["Conc"][30:40]
p_bad = minimize(Q, p_start, (time_points , y_bad)).x
y_bad = [f(t, p_bad) for t in xVar]


doses = [150, 100, 100, 50, 50, 50, 50, 50, 50]
doses = [dose/150 for dose in doses]

'''
y_bad = data["Conc"][80:90]
p_bad = minimize(Q, p_start, (time_points , y_bad)).x
y_bad = [f(t, p_bad) for t in xVar]
y = y_bad
total_dose = [0 for i in range(1000)]
total_dose = doses[0]*np.array([x+y for x,y in zip(total_dose,y)])
for i in range(1,len(doses)):
    shifted_dose = doses[i]*np.array([0 for _ in range(120*i)] + list(y)[:-120*i])
    total_dose = [x+y for x,y in zip(total_dose,shifted_dose)]
plt.plot(xVar, total_dose, color="maroon", linewidth=1.5, label="Mest känslig")

p_opt = minimize(Q, p_start, (time_points_all , y_data_all)).x
y_opt = np.array([f(t, p_opt) for t in xVar])
y = y_opt
total_dose = [0 for i in range(1000)]
total_dose = doses[0]*np.array([x+y for x,y in zip(total_dose,y)])
for i in range(1,len(doses)):
    shifted_dose = doses[i]*np.array([0 for _ in range(120*i)] + list(y)[:-120*i])
    total_dose = [x+y for x,y in zip(total_dose,shifted_dose)]
plt.plot(xVar, total_dose, color="black", linewidth=1.5, label="Genomsnitt")

y_bad2 = data["Conc"][30:40]
p_bad2 = minimize(Q, p_start, (time_points , y_bad2)).x
y_bad2 = [f(t, p_bad2) for t in xVar]
y = y_bad2
total_dose = [0 for i in range(1000)]
total_dose = doses[0]*np.array([x+y for x,y in zip(total_dose,y)])
for i in range(1,len(doses)):
    shifted_dose = doses[i]*np.array([0 for _ in range(120*i)] + list(y)[:-120*i])
    total_dose = [x+y for x,y in zip(total_dose,shifted_dose)]
plt.plot(xVar, total_dose, color="darkgreen", linewidth=1.5, label="Minst känslig")

plt.fill_between(xVar, 0, 2, alpha=0.2, color="green")
plt.fill_between(xVar, 2, 3, alpha=0.3, color="yellow")
plt.fill_between(xVar, 3, 5, alpha=0.2, color="orange")
plt.fill_between(xVar, 5, 9, alpha=0.2, color="red")
plt.axhline(y=1, linestyle="dotted", color="black")
plt.xlabel('Tid(h)')
plt.ylabel('Koncentration(mg/L)')
plt.legend()

#plt.scatter(time_points_all, y_data_all)
plt.show()
'''

'''
x1 = []
y1 = []
x2 = []
y2 = []
x3 = []
y3 = []
x4 = []
y4= []
for s in range(len(data["Symptom"])):
    if data["Symptom"][s] == 0:
        y1.append(data["Symptom"][s])
        x1.append(data["Conc"][s])
    if data["Symptom"][s] == 1:
        y2.append(data["Symptom"][s])
        x2.append(data["Conc"][s])
    if data["Symptom"][s] == 2:
        y3.append(data["Symptom"][s])
        x3.append(data["Conc"][s])
    if data["Symptom"][s] == 3:
        y4.append(data["Symptom"][s])
        x4.append(data["Conc"][s])
plt.scatter(x1,y1,alpha=0.7,color="limegreen",edgecolors="black")
plt.scatter(x2,y2,alpha=0.7,color="yellow",edgecolors="black")
plt.scatter(x3,y3,alpha=0.7,color="orange",edgecolors="black")
plt.scatter(x4,y4,alpha=0.7,color="red",edgecolors="black")
plt.ylabel('Biverknings grad')
plt.xlabel('Koncentration(mg/L)')
plt.show()
'''


"""
plt.scBatter(symptom_sum, k_el_lst)
plt.title("Summa av biverknings grad mot $k_{el}$")
plt.ylabel('$k_{el}$')
plt.xlabel('Summa av biverknings grad')
plt.show()
plt.scatter(symptom_sum, AUC_lst)
plt.title("Summa av biverknings grad mot AUC")
plt.ylabel('AUC')
plt.xlabel('Summa av biverknings grad')
plt.show()
plt.scatter(symptom_sum, CL_lst)
plt.title('CL')
plt.show()
plt.scatter(symptom_sum, MRT_lst)
plt.title("MRT")
plt.show()
plt.scatter(symptom_sum, V_ss_lst)
plt.title("V_ss")
plt.show()
"""


print(f'k_el: {k_el_lst}')
print(f'AUC: {AUC_lst}')
print(f'CL: {CL_lst}')
print(f'MRT: {MRT_lst}')
print(f'V_ss: {V_ss_lst}')
######## BOXPLOTTAR

plt.figure()
fig, axs = plt.subplots(5)
medianprops = dict(linestyle='-', linewidth=2.5, color='firebrick')
axs[0].boxplot(k_el_lst, labels=["$k_{el}$"], vert=False, widths=1, medianprops=medianprops)
#axs[0].xlabel("k_el")
axs[1].boxplot(AUC_lst, widths=1, labels=["AUC"], vert=False, medianprops=medianprops)
axs[2].boxplot(CL_lst, widths=1, labels=["$CL$"], vert=False, medianprops=medianprops)
axs[3].boxplot(MRT_lst, widths=1, labels=["$MRT$"], vert=False, medianprops=medianprops)
axs[4].boxplot(V_ss_lst, widths=1, labels=["$V_{ss}$"], vert=False, medianprops=medianprops)
plt.tight_layout()
plt.show()


