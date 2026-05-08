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
                out_arr[i,j]=(inp_arr[i,j+1]+inp_arr[i,j-1]+2*inp_arr[i+1,j])/4 #maintains neumann boundary condition: inp_arr[i-1,j]=inp_arr[i+1,j]
            else:
                out_arr[i,j]=(inp_arr[i+1,j]+inp_arr[i-1,j]+inp_arr[i,j+1]+inp_arr[i,j-1])/4+(inp_arr[i+1,j]-inp_arr[i-1,j])/(8*(i+1/2))
    return out_arr

def converge_check(inp_arr,array_y,array_x,converge_val): #checks if values of array are sufficiently converged to approximate a solution
    comp_arr=fin_diff_cylindrical(inp_arr,array_y,array_x)
    a=round(np.max(np.abs(comp_arr-inp_arr)),9)
    if a<=converge_val:
        return True,a
    return False,a

def v_plot2(inp_arr,y_axis,x_axis,position):
    V_y=[]
    y=[]
    for i in range(0,x_axis):
        V_y.append(round(inp_arr[position,i],10))
        y.append(i)
    plt.title(f'Potential at vertical midpoint of electrode(s)\n{title}\nr={position}\n{y_axis}x{x_axis} {resolution} resolution',fontsize=8,loc='left')
    plt.xlabel(f'z, mm')
    plt.ylabel(f'V(z), C/m')
    plt.plot(y,V_y)
    # plt.savefig(f"{idt} {title} POTENTIAL GRAPH {position} {array_r} {array_z}.png")
    plt.show()

def upsize(inp_arr,inp_y,inp_x,scale_factor):
    out_arr=np.zeros([scale_factor*inp_y,scale_factor*inp_x])
    for i in range(inp_y):
        for j in range(inp_x):
            out_arr[scale_factor*i:scale_factor*i+scale_factor,scale_factor*j:scale_factor*j+scale_factor]=inp_arr[i,j]
    return out_arr



elec_len=100 #=50mm
end_len=1+120 #=60mm
array_r=1+41 #=20.5mm
bound_pot_lst=[35,34,33,28,140]
end_pot=0



elec_num=len(bound_pot_lst)
N=elec_num*elec_len
array_z=elec_num*elec_len+2*end_len



mid_point='x'
upsize_val=20
converge_val=1e-5 # upper limit for the numerical solution to be considered "solved", i.e. subsequent arrays do not differ by this amount
init_array=np.full([array_r,array_z], 0.0) # creates array and initialises all values as 0
title='PENNING TRAP, CYLINDRICAL'
resolution='0.5mm'



init_array=np.full([array_r,N],0.0)
for i in range(elec_num):
    init_array[0:,i*N//elec_num:(i+1)*N//elec_num]=bound_pot_lst[i]
init_array=np.fliplr(np.append(np.fliplr(np.append(init_array,np.full([array_r,end_len],end_pot),axis=1)),np.full([array_r,end_len],end_pot),axis=1))
init_array[1:,1:array_z-1]=0
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
        v_plot2(arr,array_r,array_z,0)
        plt.imshow(arr,cmap='jet',origin='lower')
        plt.colorbar()
        arr=np.flipud(arr)
        plt.title(f'{title}\n{array_r}x{array_z} grid {resolution} resolution\n{bound_pot_lst}V\nConvergence of {converge_val} found in {i} steps\n{t}s',loc='left',fontsize=7)
        # plt.imsave(f'{idt} {title} {array_z} by {array_r} {i} steps {converge_val} {t}s RAW.png',arr,cmap='jet')
        resize_arr=upsize(arr,array_r,array_z,upsize_val)
        plt.imsave(f'{idt} {title} {array_z} by {array_r} {i} steps {converge_val} {t}s RESIZE.png',resize_arr,cmap='jet')
        # plt.savefig(f'{idt} {title} {array_z} by {array_r} {i} steps {converge_val} {t}s.png')
        plt.show()
        print(f'\n{t}s\n{idt}')
        break
    else:
        arr=fin_diff_cylindrical(arr,array_r,array_z)
    i+=1



arr2=np.concatenate([arr,np.flipud(arr)], axis=0)
plt.imshow(arr2,cmap='jet')
plt.title(f'{title} MIRRORED\n{2*array_r}x{array_z} grid {resolution} resolution\n{bound_pot_lst}V\nConvergence of {converge_val} found in {i} steps\n{t}s',loc='left',fontsize=7)
plt.colorbar()
# plt.imsave(f'{idt} MIRRORED {title} {array_z} by {2*array_r} {i} steps {converge_val} {t}s RAW.png',arr2,cmap='jet')
plt.show()
resize_arr2=upsize(arr2,2*array_r,array_z,upsize_val)
plt.imshow(resize_arr2)
# plt.imsave(f'{idt} MIRRORED {title} {array_z} by {2*array_r} {i} steps {converge_val} {t}s RESIZE.png',resize_arr2,cmap='jet')



exp_lst=[]
for i in range(array_r):
    for j in range(array_z):
        exp_lst.append([i,j,arr[i,j]])
# np.savetxt(f'{idt} {title} NUMPY EXPORT.csv', exp_lst,delimiter=',',fmt='%s',header='r,z,V')
