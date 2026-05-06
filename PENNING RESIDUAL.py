import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.special import i0

idt=np.random.randint(1,99) #file-handling number, no relevance to main program

def fin_diff_cylindrical(inp_arr,array_r,array_z):
    out_arr=np.full([array_r,array_z],0.0)
    for i in range(array_r):
        for j in range(array_z):
            if j==0 or i==array_r-1 or j==array_z-1:
                out_arr[i,j]=inp_arr[i,j] #maintains dirlichet boundary condition
            elif i==0:
                out_arr[i,j]=(inp_arr[i,j+1]+inp_arr[i,j-1]+2*inp_arr[i+1,j])/4 #maintains neumann boundary condition: inp_arr[i+1,j]=inp_arr[i-1,j]
            else:
                out_arr[i,j]=(inp_arr[i+1,j]+inp_arr[i-1,j]+inp_arr[i,j+1]+inp_arr[i,j-1])/4+(inp_arr[i+1,j]-inp_arr[i-1,j])/(8*(i+1/2))
    return out_arr

def converge_check(inp_arr,array_y,array_x,converge_val): #checks if values of array are sufficiently converged to approximate a solution
    comp_arr=fin_diff_cylindrical(inp_arr,array_y,array_x)
    a=round(np.max(np.abs(comp_arr-inp_arr)),9)
    if a<=converge_val:
        return True,a
    return False,a

def potential_sum(rad_n,zed_n):
    exp_V=v1*2/3
    for n in range(1,N):
        alpha_n=np.sin(np.pi*n*2/3)-np.sin(np.pi*n*4/3)
        exp_V+=(v1/np.pi)*alpha_n*(1/n)*(i0(2*np.pi*n*rad_n/Z_ZERO)/i0(2*np.pi*n*RADIUS_ZERO/Z_ZERO))*np.cos(2*np.pi*n*zed_n/Z_ZERO)
    return exp_V

v1=100.0
RADIUS_ZERO=20
Z_ZERO=3*RADIUS_ZERO
dim_r=RADIUS_ZERO
dim_z=Z_ZERO
N=300
GLOBAL_CMAP='jet'
z1=RADIUS_ZERO
z2=2*RADIUS_ZERO
array_r=RADIUS_ZERO
array_z=Z_ZERO
converge_val=1e-6 # upper limit for the numerical solution to be considered "solved", i.e. subsequent arrays do not differ by this amount

init_array=np.full([array_r,array_z], v1) # creates array and initialises all values as v1

init_array[0,0:z1]=v1
init_array[0,z1:z2]=0
init_array[0,z2:Z_ZERO]=v1
init_array=np.flipud(init_array)

arr=fin_diff_cylindrical(init_array,array_r,array_z)
i=1
t_a=time.time()

while True:
    chk,converge_val_n=converge_check(arr,array_r,array_z,converge_val)
    a=int(abs(np.log10(converge_val)))
    print(f'Step {i}, {100*round(converge_val/converge_val_n,a):.{a-2}f}% Convergence')
    if chk==True:
        t=round(time.time()-t_a, 10)
        plt.imshow(arr,cmap=GLOBAL_CMAP,origin='lower')
        plt.colorbar()
        plt.xlabel(f'z')
        plt.ylabel(f'r')
        plt.title(f'SIMPLE PENNING, NUMERICAL SOLUTION\n{RADIUS_ZERO}x{Z_ZERO}\n{v1}V -> GROUND\nConvergence of {converge_val} found in {i} steps, {t}s',loc='left',fontsize=7)
        # plt.savefig(f'{idt} NUMERICAL SIMPLE PENNING {v1}V {RADIUS_ZERO}.png')
        plt.show()
        break
    else:
        arr=fin_diff_cylindrical(arr,array_r,array_z)
    i+=1

arr2=np.full([dim_r,dim_z],0.0)

for i in range(0,dim_r):
    for j in range(0,dim_z):
        s=potential_sum(i,j)
        # print(i,j,s)
        arr2[i,j]=s

print(v1)
print(RADIUS_ZERO)

plt.imshow(arr2,cmap=GLOBAL_CMAP,origin='lower')
plt.xlabel(f'z')
plt.ylabel(f'r')
plt.colorbar()
plt.title(f'SIMPLE PENNING TRAP, ANALYTIC SOLUTION\nn={N}, {RADIUS_ZERO}x{Z_ZERO},\n{v1}V -> GROUND',loc='left',fontsize=7)
# plt.savefig(f'{idt} ANALYTIC SIMPLE PENNING {v1}V {RADIUS_ZERO}.png')
plt.show()

plt.imshow(arr2-arr,cmap=GLOBAL_CMAP,origin='lower')
plt.xlabel(f'z')
plt.ylabel(f'r')
plt.colorbar()
avg=np.average(arr2-arr)
plt.title(f'SIMPLE PENNING TRAP, RESIDUAL\nn={N}, {RADIUS_ZERO}x{Z_ZERO},\n{v1}V -> GROUND\nAverage: {avg}',loc='left',fontsize=7)
# plt.savefig(f'{idt} RESIDUAL SIMPLE PENNING {v1}V {RADIUS_ZERO}.png')
plt.show()
print(avg)

plt.imshow(np.abs(arr2-arr),cmap=GLOBAL_CMAP,origin='lower')
plt.xlabel(f'z')
plt.ylabel(f'r')
plt.colorbar()
avg=np.average(np.abs(arr2-arr))
plt.title(f'SIMPLE PENNING TRAP, ABSOLUTE RESIDUAL\nn={N}, {RADIUS_ZERO}x{Z_ZERO},\n{v1}V -> GROUND\nAverage: {avg}',loc='left',fontsize=7)
# plt.savefig(f'{idt} ABSOLUTE RESIDUAL SIMPLE PENNING {v1}V {RADIUS_ZERO}.png')
plt.show()
print(avg)

print(f'\n{idt}')
