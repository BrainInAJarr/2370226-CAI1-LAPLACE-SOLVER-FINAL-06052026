import numpy as np
import time
import matplotlib.pyplot as plt
np.set_printoptions(precision=9)

idt=np.random.randint(1,99) #file-handling number, no relevance to main program

N=1001
RAD_GLOBAL=50
dim=3+2*RAD_GLOBAL
elec_pot_1=0 #bottom electrode
elec_pot_2=1.0 #top electrode
converge_val=1e-5
GLOBAL_CMAP='jet'

def gen_circ(array_y,array_x,potential_1=1,potential_2=-1):
    elec_points=[]
    arr=np.full([array_y,array_x],0.0)
    x_mid=array_x//2
    y_mid=array_y//2
    for i in range(0,array_y):
        for j in range(0,array_x):
            q=np.sqrt((i-y_mid)**2+(j-x_mid)**2)
            if RAD_GLOBAL<q:
                    if i<y_mid:
                        arr[i,j]=potential_1
                    elif i>y_mid:
                        arr[i,j]=potential_2
                    elec_points.append([i,j])
            else:
                arr[i,j]=0
    return arr,elec_points

def fin_diff_circ(inp_arr,array_y,array_x):#executes one interation of finite difference method
    out_arr=np.empty([array_y,array_x])
    for i in range(array_y):
        for j in range(array_x):
            q=np.sqrt((i-array_y//2)**2+(j-array_x//2)**2)
            if i==0 or j==0 or i==array_y-1 or j==array_x-1 or RAD_GLOBAL<q<=RAD_GLOBAL+1:#maintains boundary condition(s)
                out_arr[i,j]=inp_arr[i,j]
            else:
                out_arr[i,j]=(inp_arr[i,j+1]+inp_arr[i,j-1]+inp_arr[i+1,j]+inp_arr[i-1,j])/4 #approximates second CARTESIAN derivative at [i,j], i.e. takes average of 4 nearest neighbours 
    return out_arr

def converge_check(inp_arr,array_y,array_x,converge_val): #checks if values of array are sufficiently converged to approximate a solution
    comp_arr=fin_diff_circ(inp_arr,array_y,array_x)
    a=round(np.max(np.abs(comp_arr-inp_arr)),9)
    if a<=converge_val:
        return True,a
    return False,a

def cleanup(inp_arr,array_y,array_x,bound_pot1,bound_pot2):
    out_arr=inp_arr
    y_mid=array_y//2
    x_mid=array_x//2
    for i in range(array_y):
        for j in range(array_x):
            q=np.sqrt((i-y_mid)**2+(j-x_mid)**2)
            if q>RAD_GLOBAL:
                out_arr[i,j]=np.nan
    return out_arr

def potential_sum(rad_n,theta_n,N):
    exp_V=elec_pot_2/2
    for n in range(1,N,2):
        exp_V+=(2*(elec_pot_2)/np.pi)*(1/n)*((rad_n/RAD_GLOBAL)**n)*np.sin(n*theta_n)
    return exp_V

def potential_avg(inp_arr):
    avg_sum=0
    avg_i=0
    for i in range(dim):
        for j in range(dim):
            if np.isnan(inp_arr[i,j])==False:
            # if np.isnan(inp_arr[i,j])==False and np.sqrt((i-dim//2)**2+(j-dim//2)**2)<RAD_GLOBAL-1:    
                avg_sum+=inp_arr[i,j]
                avg_i+=1
    return avg_sum/avg_i

init_array,elec_points=gen_circ(dim,dim,elec_pot_1,elec_pot_2)

arr=fin_diff_circ(init_array,dim,dim)
i=1
t_a=time.time()

while True:
    check,converge_val_n=converge_check(arr,dim,dim,converge_val)
    a=int(abs(np.log10(converge_val)))
    print(f'Step {i}, {100*round(converge_val/converge_val_n,a):.{a-2}f}% Convergence')
    if check==True:
        t=round(time.time()-t_a, 10)
        arr=cleanup(arr,dim,dim,elec_pot_1,elec_pot_2)
        plt.imshow(arr,cmap=GLOBAL_CMAP,origin='lower')
        plt.colorbar()
        plt.xlabel(f'x')
        plt.ylabel(f'y')
        plt.title(f'SEMICIRCULAR ELECTRODES, NUMERICAL SOLUTION\n{dim}x{dim}, Radius: {RAD_GLOBAL}\n{elec_pot_2}V -> GROUND\nConvergence of {converge_val} found in {i} steps, {t}s',loc='left',fontsize=7)
        # plt.savefig(f'{idt} SIMULATED SEMICIRCULAR ELECTRODES {elec_pot_2}V {RAD_GLOBAL}.png',dpi=600,transparent=True)
        plt.show()
        print(f'\n{t}')#s')
        break
    else:
        arr=fin_diff_circ(arr,dim,dim)
    i+=1

pot_arr=np.full([dim,dim], 0.0)
x_m=dim//2
y_m=dim//2
for i in range(0,dim):
    for j in range(0,dim):
        x_n=j-x_m
        y_n=i-y_m
        rad=((x_n)**2+(y_n)**2)**(1/2)
        the=np.atan2(y_n,x_n)
        if rad>RAD_GLOBAL:
            pot_arr[i,j]=np.nan
        else:
            pot_arr[i,j]=potential_sum(rad,the,N)

# print(elec_pot_2)
# print(RAD_GLOBAL)

plt.imshow(pot_arr,cmap=GLOBAL_CMAP,origin='lower')
plt.colorbar()
plt.xlabel(f'x')
plt.ylabel(f'y')
plt.title(f'SEMICIRCULAR ELECTRODES, ANALYTIC SOLUTION\nn={N}, Radius: {RAD_GLOBAL}, {dim}x{dim}\n{elec_pot_2}V -> GROUND',loc='left',fontsize=7)
# plt.savefig(f'{idt} ANALYTIC SEMICIRCULAR ELECTRODES {elec_pot_2}V {RAD_GLOBAL}.png',dpi=600,transparent=True)
plt.show()

plt.imshow(pot_arr-arr,cmap=GLOBAL_CMAP,origin='lower')
plt.colorbar()
plt.xlabel(f'x')
plt.ylabel(f'y')
avg=potential_avg(pot_arr-arr)
print(avg)
plt.title(f'SEMICIRCULAR ELECTRODES, RESIDUAL\n{dim}x{dim}, Radius: {RAD_GLOBAL}\n{elec_pot_2}V -> GROUND\nAverage: {avg}',loc='left',fontsize=7)
# plt.savefig(f'{idt} RESIDUAL SEMICIRCULAR ELECTRODES {elec_pot_2}V {RAD_GLOBAL}.png',dpi=600,transparent=True)
plt.show()

plt.imshow(np.abs(arr-pot_arr),cmap=GLOBAL_CMAP,origin='lower')
plt.colorbar()
plt.xlabel(f'x')
plt.ylabel(f'y')
avg=potential_avg(np.abs(arr-pot_arr))
print(avg)
plt.title(f'SEMICIRCULAR ELECTRODES, ABSOLUTE RESIDUAL\n{dim}x{dim}, Radius: {RAD_GLOBAL}\n{elec_pot_2}V -> GROUND\nAverage: {avg}',loc='left',fontsize=7)
# plt.savefig(f'{idt} ABSOLUTE RESIDUAL SEMICIRCULAR ELECTRODES {elec_pot_2}V {RAD_GLOBAL}.png')
plt.show()

plt.imshow(100*(np.abs((pot_arr-arr)/pot_arr)),cmap=GLOBAL_CMAP,origin='lower')
plt.xlabel(f'z')
plt.ylabel(f'r')
plt.colorbar()
avg=potential_avg(100*np.abs((pot_arr-arr)/pot_arr))
plt.title(f'SEMICIRCULAR ELECTRODES, RESIDUAL PERCENTAGE ERROR\n{dim}x{dim}, Radius: {RAD_GLOBAL}\n{elec_pot_2}V -> GROUND\nAverage: {avg}%',loc='left',fontsize=7)
# plt.savefig(f'{idt} RESIDUAL ERROR SEMICIRCULAR ELECTRODES {elec_pot_2}V {RAD_GLOBAL}.png',dpi=600,transparent=True)
plt.show()

print(f'{avg}')#%')
print(f'\n{RAD_GLOBAL}\n{idt}')
