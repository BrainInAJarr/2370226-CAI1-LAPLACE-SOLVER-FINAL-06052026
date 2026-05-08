import numpy as np
import time
import matplotlib.pyplot as plt
import os
import imageio.v2 as imageio
import shutil

np.set_printoptions(precision=9)

idt=np.random.randint(1,999) #file-handling number

RAD_GLOBAL=35
dim=3+2*RAD_GLOBAL
elec_pot_1=0 #bottom electrode
elec_pot_2=1.0 #top electrode
converge_val=1e-5
GLOBAL_CMAP='Blues'
step=5

frame_directory=f'{idt} FRAME DIR'
os.mkdir(f'{frame_directory}')

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
                arr[i,j]=(potential_1+potential_2)/2
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

init_array,elec_points=gen_circ(dim,dim,elec_pot_1,elec_pot_2)

arr=fin_diff_circ(init_array,dim,dim)
plt.imshow(arr,cmap=GLOBAL_CMAP,origin='lower')
plt.title(f'radius: {RAD_GLOBAL}\nframe {0//step} (i={0})\n0% convergence',fontsize=7,loc='left')
plt.colorbar()
plt.savefig(f'{frame_directory}/{0}')

i=1

t_a=time.time()

while True:
    check,converge_val_n=converge_check(arr,dim,dim,converge_val)
    if check==True:
        t=round(time.time()-t_a, 10)
        plt.imshow(arr,cmap=GLOBAL_CMAP,origin='lower')
        plt.title(f'radius: {RAD_GLOBAL}\nframe {i//step} (i={i})\n{100*round(converge_val/converge_val_n,a):.{4}f}% convergence',fontsize=7,loc='left')
        plt.savefig(f'{frame_directory}/{i}.png')
        print(f'\n{t}s\n{idt}')
        break
    elif i%step==0:
        a=int(abs(np.log10(converge_val)))
        print(f'Step {i}, {100*round(converge_val/converge_val_n,a):.{a-2}f}% Convergence')    
        plt.imshow(arr,cmap=GLOBAL_CMAP,origin='lower')
        plt.title(f'radius: {RAD_GLOBAL}\nframe {i//step} (i={i})\n{100*round(converge_val/converge_val_n,a):.{4}f}% convergence',fontsize=7,loc='left')
        plt.savefig(f'{frame_directory}/{i}.png')
    else:
        arr=fin_diff_circ(arr,dim,dim)
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

imageio.mimsave(f'SEMICIRCULAR NUMERICAL SOLUTION {idt} {len(frames)} FRAMES.gif',frames,duration=0.001,loop=0)
shutil.rmtree(frame_directory)
