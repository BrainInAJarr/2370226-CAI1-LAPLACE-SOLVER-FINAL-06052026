import numpy as np
import time
import matplotlib.pyplot as plt
np.set_printoptions(precision=9)



idt=np.random.randint(1,99) #file-handling number, no relevance to main program



n=19
array_x=n
array_y=n
elec_pot_1=0 #left electrode
elec_pot_2=1 #right electrode
converge_val=1e-8
upsize_val=20



def gen_circ(array_y,array_x,potential_1,potential_2):
    elec_points=[]
    # r=min(array_x,array_y)//2
    # r=10
    # r=n//2
    r=9
    arr=np.full([array_y,array_x],abs((elec_pot_1-elec_pot_2)/2))
    x_mid=array_x//2
    y_mid=array_y//2
    for i in range(0,array_y):
        for j in range(0,array_x):
            q=np.sqrt((i-y_mid)**2+(j-x_mid)**2)
            if r-1<=q<r+1:
                    # if i<y_mid:
                    if j<x_mid:
                        arr[i,j]=potential_1
                    # elif i>y_mid:
                    elif j>x_mid:
                        arr[i,j]=potential_2
                    elec_points.append([i,j])
            elif q>r:
                if j<x_mid:
                    arr[i,j]=elec_pot_1
                elif j>x_mid:
                    arr[i,j]=elec_pot_2
            # else:
            #     arr[i,j]=abs((elec_pot_2+elec_pot_1)/2)
    return arr,elec_points,r

def upsize(inp_arr,inp_y,inp_x,scale_factor):
    out_arr=np.zeros([scale_factor*inp_y,scale_factor*inp_x])
    for i in range(inp_y):
        for j in range(inp_x):
            out_arr[scale_factor*i:scale_factor*i+scale_factor,scale_factor*j:scale_factor*j+scale_factor]=inp_arr[i,j]
    return out_arr

def fin_diff_circ(inp_arr,array_y,array_x):#executes one interation of finite difference method
    out_arr=np.empty([array_y,array_x])
    for i in range(array_y):
        for j in range(array_x):
            # r=min(array_x,array_y)//2.2
            r=n//2
            q=np.sqrt((i-array_y//2)**2+(j-array_x//2)**2)
            if i==0 or j==0 or i==array_y-1 or j==array_x-1:# or r<q<=r+1:#maintains boundary condition(s)
                out_arr[i,j]=inp_arr[i,j]
            elif q>r:
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

def cleanup(inp_arr,array_y,array_x,bound_pot1,bound_pot2):
    radius=min(array_x,array_y)//2.2
    out_arr=inp_arr
    y_mid=array_y//2
    x_mid=array_x/2
    for i in range(array_y):
        for j in range(array_x):
            q=np.sqrt((i-y_mid)**2+(j-x_mid)**2)
            if q>radius:
                # if i>=y_mid:
                if j>=x_mid:
                    out_arr[i,j]=bound_pot2
                else:
                    out_arr[i,j]=bound_pot1
    return out_arr

# def csv_export(inp_arr,array_y,array_x,title=''):
#     exp_lst=[]
#     for i in range(array_x):
#         for j in range(array_y):
#             exp_lst.append([i,j,inp_arr[j,i]])
#     np.savetxt(f'{idt} {title} NUMPY EXPORT.csv', exp_lst,delimiter=',',fmt='%s',header='x,y,V')



init_array,elec_points,elec_radius=gen_circ(array_y,array_x,elec_pot_1,elec_pot_2)
title='CONDUCTIVE PAPER SIMULATION'
mid_point='x'



plt.imshow(init_array,cmap='jet',origin='lower')
plt.colorbar()
plt.title(f'{title} {array_x} by {array_y}\nRadius: {elec_radius}\nINITIAL',loc='left',fontsize=7)
# plt.imsave(f'{idt} {title} {array_x} by {array_y} INITIAL RAW.png',init_array,cmap='jet')
# resize_init_array=upsize(init_array,array_y,array_x,upsize_val)
# plt.imsave(f'{idt} {title} {array_x} by {array_y} INITIAL RESIZE.png',resize_init_array,cmap='jet')
# resize_init_array=upsize(init_array,array_y,array_x,upsize_val)
# plt.savefig(f'{idt} {title} {array_x} by {array_y} INITIAL.png')
plt.show()



arr=fin_diff_circ(init_array,array_y,array_x)
i=1
t_a=time.time()



while True:
    check,converge_val_n=converge_check(arr,array_y,array_x,converge_val)
    a=int(abs(np.log10(converge_val)))
    print(f'Step {i} {100*round(converge_val/converge_val_n,a):.5f}% Convergence')
    if check==True:
        t=round(time.time()-t_a, 10)
        plt.imshow(arr,cmap='jet',origin='lower')
        plt.colorbar()
        plt.title(f'{title}\n{array_y}x{array_x} grid, radius:{elec_radius}, {elec_pot_1}->{elec_pot_2}V\nConvergence of {converge_val} found in {i} steps\n{t}s',loc='left',fontsize=7)
        # resize_arr=upsize(arr,array_y,array_x,upsize_val)
        # plt.imsave(f'{idt} {title} {array_x} by {array_y} {i} steps {elec_pot_1}V {elec_pot_2}V {converge_val} {t}s.png',resize_arr,cmap='jet')
        plt.show()
        print(np.size(arr))
        # csv_export(arr,array_y,array_x,f'{title} {array_y} {array_x}')
        arr=np.delete(arr,n-1,axis=1)
        arr=np.delete(arr,n-1,axis=0)
        arr=np.delete(arr,0,axis=1)
        arr=np.delete(arr,0,axis=0)
        arr=np.delete(arr,n-3,axis=1)
        arr=np.delete(arr,n-3,axis=0)
        arr=np.delete(arr,0,axis=1)
        arr=np.delete(arr,0,axis=0)
        print(np.size(arr))
        plt.imshow(arr,cmap='jet',origin='lower')
        plt.colorbar()
        plt.title(f'{title}CUT\n{array_y-4}x{array_x-4} grid, radius:{elec_radius}, {elec_pot_1}->{elec_pot_2}V\nConvergence of {converge_val} found in {i} steps\n{t}s',loc='left',fontsize=7)
        # resize_arr=upsize(arr,array_y-4,array_x-4,upsize_val)
        # plt.imsave(f'{idt} CUT {title} {array_x-4} by {array_y-4} {i} steps {elec_pot_1}V {elec_pot_2}V {converge_val} {t}s.png',resize_arr,cmap='jet')
        plt.show()
        # csv_export(arr,array_y-4,array_x-4,f'{title} CUT {array_y-4} {array_x-4}')
        print(f'\n{t}s\n{idt}')
        break
    else:
        arr=fin_diff_circ(arr,array_y,array_x)
    i+=1
