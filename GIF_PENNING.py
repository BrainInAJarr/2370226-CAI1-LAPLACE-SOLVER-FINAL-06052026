import numpy as np
import time
import matplotlib.pyplot as plt
import os
import imageio.v2 as imageio
import shutil

np.set_printoptions(precision=9)

idt=np.random.randint(1,999) #file-handling number

v1=1.0
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
step=5

frame_directory=f'{idt} FRAME DIR'
os.mkdir(f'{frame_directory}')

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

init_array=np.full([array_r,array_z], v1) # creates array and initialises all values as v1

init_array=np.full([array_r,array_z], v1) # creates array and initialises all values as v1

init_array[0,0:z1]=v1
init_array[0,z1:z2]=0
init_array[0,z2:Z_ZERO]=v1
init_array=np.flipud(init_array)

plt.imshow(init_array,cmap=GLOBAL_CMAP,origin='lower')
plt.title(f'radius: {RADIUS_ZERO}\nframe {0//step} (i={0})\n0% convergence',fontsize=7,loc='left')
plt.colorbar()
plt.savefig(f'{frame_directory}/{0}')

arr=fin_diff_cylindrical(init_array,array_r,array_z)
i=1
t_a=time.time()

while True:
    check,converge_val_n=converge_check(arr,array_r,array_z,converge_val)    
    if check==True:
        t=round(time.time()-t_a, 10)
        plt.imshow(arr,cmap=GLOBAL_CMAP,origin='lower')
        plt.title(f'radius: {RADIUS_ZERO}\nframe {i//step} (i={i})\n{100*round(converge_val/converge_val_n,a):.{4}f}% convergence',fontsize=7,loc='left')
        plt.savefig(f'{frame_directory}/{i}.png')
        print(f'\n{t}s\n{idt}')
        break
    elif i%step==0:
        a=int(abs(np.log10(converge_val)))
        print(f'Step {i}, {100*round(converge_val/converge_val_n,a):.{a-2}f}% Convergence')
        plt.imshow(arr,cmap=GLOBAL_CMAP,origin='lower')
        plt.title(f'radius: {RADIUS_ZERO}\nframe {i//step} (i={i})\n{100*round(converge_val/converge_val_n,a):.{4}f}% convergence',fontsize=7,loc='left')
        plt.savefig(f'{frame_directory}/{i}.png')
    else:
        arr=fin_diff_cylindrical(arr,array_r,array_z)
    i+=1

f=os.listdir(frame_directory)
f2=[]
frames=[]

for i in f:
    f2.append(int(i.replace('.png','')))

for i in sorted(f2):
    j=os.path.join(frame_directory,str(i))
    print(f'j: {j}')
    frames.append(imageio.imread(f'{j}.png'))

duration=len(frames)
print(f'{duration} frames\n{idt}')

print(len(frames))

imageio.mimsave(f'SIMPLE PENNING SOLUTION {idt} {len(frames)} FRAMES.gif',frames,duration=0.001,loop=0)
shutil.rmtree(frame_directory)
