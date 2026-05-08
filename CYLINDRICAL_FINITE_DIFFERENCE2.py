import numpy as np
import time
import matplotlib.pyplot as plt



idt=np.random.randint(1,99) #file-handling number, no relevance to main program



def fin_diff_cylindrical(inp_arr,array_r,array_z):
    out_arr=np.empty([array_r,array_z])
    for i in range(array_r):
        for j in range(array_z):
            if j==0 or i==array_r-1 or j==array_z-1:
                out_arr[i,j]=inp_arr[i,j] #maintains dirlichet boundary condition
            elif i==0:
                out_arr[i,j]=(inp_arr[i,j+1]+inp_arr[i,j-1]+2*inp_arr[i+1,j])/4 #maintains neumann boundary condition
            else:
                out_arr[i,j]=(inp_arr[i+1,j]+inp_arr[i-1,j]+inp_arr[i,j+1]+inp_arr[i,j-1])/4+(inp_arr[i+1,j]-inp_arr[i-1,j])/(8*(i+1/2))
    return out_arr

def converge_check(inp_arr,array_y,array_x,converge_val): #checks if values of array are sufficiently converged to approximate a solution
    comp_arr=fin_diff_cylindrical(inp_arr,array_y,array_x)
    a=round(np.max(np.abs(comp_arr-inp_arr)),9)
    if a<=converge_val:
        return True,a
    return False,a

def v_plot(inp_arr,y_axis,x_axis,selection='x'):
    y=list()
    V_y=list()
    if selection=='y':
        for i in range(0,y_axis):
            V_y.append(round(inp_arr[i,x_axis//2],10))
            y.append(i)
        plt.title(f'Potential at horizontal midpoint of electrode(s)\n{title}\n{y_axis}x{x_axis}',fontsize=8,loc='left')
    else:
        for i in range(0,x_axis):
            V_y.append(round(inp_arr[y_axis//2,i],10))
            y.append(i)
        plt.title(f'Potential at vertical midpoint of electrode(s)\n{title}\n{y_axis}x{x_axis}',fontsize=8,loc='left')
    plt.xlabel(f'{selection}, m')
    plt.ylabel(f'V({selection}), C/m')
    plt.plot(y,V_y)

def upsize(inp_arr,inp_y,inp_x,scale_factor):
    out_arr=np.zeros([scale_factor*inp_y,scale_factor*inp_x])
    for i in range(inp_y):
        for j in range(inp_x):
            out_arr[scale_factor*i:scale_factor*i+scale_factor,scale_factor*j:scale_factor*j+scale_factor]=inp_arr[i,j]
    return out_arr



array_z=2+50
array_r=2+50
bound_pot_1=1 #left electrode
bound_pot_2=-1 #right electrode
mid_point='x'
upsize_val=20
converge_val=1e-8 #upper limit for the numerical solution to be considered "solved", i.e. subsequent arrays do not differ by this amount
init_array=np.full([array_r,array_z], 0) #creates array and initialises all values as 0
title='placeholder'



init_array[0:,0]=bound_pot_1
init_array[0:,array_z-1]=bound_pot_2
title='PARALLEL ELECTRODES, CYLINDRICAL'



plt.imshow(init_array,cmap='jet',origin='lower')
plt.colorbar()
plt.title(f'{title} {array_z} by {array_r} INITIAL',loc='left',fontsize=7)
plt.imsave(f'{idt} {title} {array_z} by {array_r} INITIAL RAW.png',init_array,cmap='jet')
resize_init_array=upsize(init_array,array_r,array_z,upsize_val)
plt.imsave(f'{idt} {title} {array_z} by {array_r} INITIAL RESIZE.png',resize_init_array,cmap='jet')
resize_init_array=upsize(init_array,array_r,array_z,upsize_val)
plt.savefig(f'{idt} {title} {array_z} by {array_r} INITIAL.png')
plt.show()



arr=fin_diff_cylindrical(init_array,array_r,array_z)
i=1
t_a=time.time()



while True:
    chk,converge_val_n=converge_check(arr,array_r,array_z,converge_val)
    a=int(abs(np.log10(converge_val)))
    print(f'Step {i}, {100*round(converge_val/converge_val_n,a):.6f}% Convergence')
    if chk==True:
        t=round(time.time()-t_a, 10)
        v_plot(arr,array_r,array_z,mid_point)
        plt.savefig(f"{idt} {title} POTENTIAL GRAPH {array_r} {array_z}.png")
        plt.show()
        plt.imshow(arr,cmap='jet',origin='lower')
        plt.colorbar()
        arr=np.flipud(arr)
        plt.title(f'{title}\n{array_r}x{array_z} grid, {bound_pot_1}->{bound_pot_2}V\nConvergence of {converge_val} found in {i} steps\n{t}s',loc='left',fontsize=7)
        plt.imsave(f'{idt} {title} {array_z} by {array_r} {i} steps {bound_pot_1} {bound_pot_2}V {converge_val} {t}s RAW.png',arr,cmap='jet')
        resize_arr=upsize(arr,array_r,array_z,upsize_val)
        plt.imsave(f'{idt} {title} {array_z} by {array_r} {i} steps {bound_pot_1} {bound_pot_2}V {converge_val} {t}s RESIZE.png',resize_arr,cmap='jet')
        plt.savefig(f'{idt} {title} {array_z} by {array_r} {i} steps {bound_pot_1} {bound_pot_2}V {converge_val} {t}s.png')
        plt.show()
        print(f'\n{t}s\n{idt}')
        break
    else:
        arr=fin_diff_cylindrical(arr,array_r,array_z)
    i+=1
