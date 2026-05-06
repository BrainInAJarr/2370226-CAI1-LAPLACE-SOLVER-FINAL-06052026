import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.integrate import trapezoid
import scipy.constants

"""
TRAP ELECTRODES
    35,34,33,28,V_5
    V_5=80,90,100,110,120,130,140

B=40mt
Q=e
M=m_e

w_c=QB/M
w_z=(2QU_0/Md^2)**0.5
"""

idt_list=[50,79,57,20,39,88,90] #V_5 = 80,90,100,110,120,130,140

def numerical_frequency(ident,well_pot):

    file_path=f"C:/Users/joehu/Documents/Physics BSc (Hons) F300/PH-311 BSc PROJECT/CODE/SIMULATION DATA/{ident} PENNING TRAP, CYLINDRICAL NUMPY EXPORT.csv"

    with open(file_path, newline='') as csvfile:
        arr = list(csv.reader(csvfile))

    arr=arr[1:]

    volt_arr=max([float(i[2]) for i in arr])

    z=list()
    V_z=list()

    for i in arr:
        if float(i[0])==41.0:
            if 200<float(i[1])<500:
                z.append(float(i[1]))
                V_z.append(round(float(i[2]),6))
    
    volt_min=min(V_z)

    well_V=list()
    well_z=list()

    for i in range(0,len(V_z)):
        if volt_min<=V_z[i]<well_pot:
            well_V.append(V_z[i])
            well_z.append(i)

    # plt.scatter(well_z,well_V)
    # plt.title(f'{volt_arr} potential well\n{well_pot}V\n{min(well_z)}<=z<={max(well_z)}',fontsize=7,loc='left')
    # plt.savefig(f'{volt_arr} potential well {well_pot}V {min(well_z)} {max(well_z)}.png')
    # plt.show()

    well_z_int=np.linspace(min(well_z),max(well_z),100)
    well_V_int=np.interp(well_z_int,well_z,well_V)

    const_coefficient=2*scipy.constants.elementary_charge/scipy.constants.electron_mass #Coulombs per kg
    
    integrand = lambda j: (const_coefficient*(well_pot-j))**(-1/2)

    axial_bounce_w_calc = lambda f: 4*scipy.constants.pi/f

    well_V=[integrand(i) for i in well_V_int]
    well_z=[i/100 for i in well_z_int]

    V_integral_trapz=trapezoid(well_V,well_z)

    axial_bounce_trapz=axial_bounce_w_calc(V_integral_trapz)
    
    print(f'{well_pot}V {volt_arr}V, w={axial_bounce_trapz/1e6} Mrad/s')
    
    return [volt_arr, float(axial_bounce_trapz/1e6)]

a=np.linspace(32,33,2+4)
y_exp=[8.52,8.78,8.97,9.13,9.24,9.18,9.55]

#idt: 50,79,57, 20, 39, 88, 90
#vlt: 80,90,100,110,120,130,140

for j in a:
    print(f'\nWELL POTENTIAL: {j}V')
    pot_list=[numerical_frequency(i,j) for i in idt_list]

    y=[i[1] for i in pot_list]
    x=[i[0] for i in pot_list]

    plt.plot(x,y,label=f'Numerical, {j}V',linestyle='--')
    plt.scatter(x,y,marker='x')

plt.scatter(x,y_exp,label='beamline',marker='x',color='black',s=60)
plt.title(f'AXIAL BOUNCE FREQUENCIES',loc='left')
plt.xlabel('Gate voltage, V')
plt.ylabel('Axial bounce frequency, MRad/s')
plt.ylim(0,max(max(y),max(y_exp))+2)
plt.legend(loc='lower right',fontsize=7)
# plt.savefig('NUMERCICAL BOUNCE FREQUENCY PLOT.png')
plt.show()
