import numpy as np
import time
import matplotlib.pyplot as plt



idt=np.random.randint(1,99) #file-handling number, no relevance to main program



def set_boundary():
    while True:
        dimension_y=int(input('INPUT Y:\n'))
        dimension_x=int(input('INPUT X:\n'))
        out_arr=np.full([dimension_y,dimension_x], 0)
        selection=str(input('SELECT BOUNDARY CONDITION:\n1)Electrode at end\n2)Parallel plates\n3)Circular electrodes\n'))
        if selection=='1':
            out_arr[0:,dimension_x-1]=bound_pot
            return out_arr,dimension_y,dimension_x
        elif selection=='2':
            out_arr[0,0:]=bound_pot
            out_arr[dimension_y-1,0:]=-bound_pot
            return out_arr,dimension_y,dimension_x
        elif selection=='3':
            return gen_circ(dimension_y,True),dimension_y,dimension_x
        else:
            selection=input('INVALID SELECTION\n1)Electrode at end\n2)Parallel plates\n3)Circular electrodes\n')

def gen_circ(array_y,array_x,potential_1=1,potential_2=-1,full=False):
    r=min(array_x,array_y)//2.5
    arr=np.full([array_y,array_x],0)
    x_mid=array_x//2
    y_mid=array_y//2
    for i in range(0,array_y):
        for j in range(0,array_x):
            q=np.sqrt((i-y_mid)**2+(j-x_mid)**2)
            if full==True:
                if r<q+2:
                    if i<y_mid:
                        arr[i,j]=potential_1
                    else:
                        arr[i,j]=potential_2
                else:
                    arr[i,j]=0
            else:
                if r<q<=r+1:
                    if i>=y_mid:
                        arr[i,j]=potential_1
                    else:
                        arr[i,j]=potential_2
                else:
                    arr[i,j]=0
    return arr

def fin_diff(inp_arr,array_y,array_x,bound_pot): #executes one interation of finite difference method
    out_arr=np.empty([array_y,array_x])
    for i in range(array_y):
        for j in range(array_x):
            if i==0 or j==0 or i==array_y-1 or j==array_x-1 or inp_arr[i,j]==bound_pot:# or inp_arr[i,j]==-bound_pot: #maintains boundary condition(s) and avoids IndexError
                out_arr[i,j]=inp_arr[i,j] 
            else:
                out_arr[i,j]=(inp_arr[i,j+1]+inp_arr[i,j-1]+inp_arr[i+1,j]+inp_arr[i-1,j])/4 #approximates second CARTESIAN derivative at [i,j], i.e. takes average of 4 nearest neighbours 
    return out_arr

def converge_check(inp_arr,array_y,array_x,converge_val): #checks if values of array are sufficiently converged to approximate a solution
    comp_arr=fin_diff(inp_arr,array_y,array_x,bound_pot)
    # a=round(np.max(np.abs(comp_arr-inp_arr)),10)
    a=np.max(np.abs(inp_arr-comp_arr))
    if a<converge_val:
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

def cleanup(inp_arr,array_y,array_x,bound_pot1,bound_pot2):
    radius=min(array_x,array_y)//2.5
    out_arr=inp_arr
    y_mid=array_y//2
    x_mid=array_x/2
    for i in range(array_y):
        for j in range(array_x):
            q=np.sqrt((i-y_mid)**2+(j-x_mid)**2)
            if q>radius+2:
                if i>=y_mid:
                    out_arr[i,j]=bound_pot2
                else:
                    out_arr[i,j]=bound_pot1
    return out_arr



# array_x=165
# array_y=116
array_x=60
array_y=40
bound_pot=1
bound_pot_1=10
bound_pot_2=-10
mid_point='y'
upsize_val=10
converge_val=1e-3 #upper limit for the numerical solution to be considered "solved", i.e. subsequent arrays do not differ by this amount
init_array=np.full([array_y,array_x], 0.0) #creates array and initialises all values as 0
title='placeholder'
# GLOBAL_CMAP='RdYlBu'
# GLOBAL_CMAP='Blues'
GLOBAL_CMAP='winter'


# init_array=gen_circ(array_y,array_x,bound_pot_1,bound_pot_2,False) #creates circular electrode with two seperate potentials
# title='CIRCULAR ELECTRODES, CARTESIAN'
# mid_point='y'

# init_array[0:,array_x-1]=bound_pot #sets right edge of array to 1
# title='RIGHT ELECTRODE, CARTESIAN'
# mid_point='x'

# init_array[0:,0]=bound_pot #sets left edge of array to 1
# title='LEFT ELECTRODE, CARTESIAN'
# mid_point='x'

# init_array[array_y//6:5*array_y//6+1,array_x//6]=-bound_pot #similar to above, but smaller and moved closer to centre
# title='SMALL LEFT ELECTRODE, CARTESIAN'
# mid_point='x'

init_array[array_y//6:5*array_y//6+1,5*array_x//6]=bound_pot #similar to above, but smaller and moved closer to centre
title='SMALL RIGHT ELECTRODE, CARTESIAN'
mid_point='x'

# init_array[0,0:]=bound_pot #parallel plate electrodes
# init_array[array_y-1,0:]=-bound_pot
# title='PARALLEL ELECTRODES, CARTESIAN'
# mid_point='y'

# init_array[array_y//2:array_y//2+1,array_x//2:array_x//2+2]=bound_pot #small electrode in centre
# init_array[array_y//2+1:array_y//2+2,array_x//2:array_x//2+2]=-bound_pot
# title='CENTRE ELECTRODE, CARTESIAN'
# mid_point='y'

# init_array[array_y//6:5*array_y//6,array_x//2+1]=bound_pot #long electrode in centre
# init_array[array_y//6:5*array_y//6,array_x//2]=bound_pot
# title='TALL CENTRE ELECTRODE, CARTESIAN'
# mid_point='x'



# plt.imshow(init_array,cmap='jet',origin='lower')
# plt.colorbar()
# plt.title(f'{title} {array_x} by {array_y} INITIAL',loc='left',fontsize=7)
# plt.imsave(f'{idt} {title} {array_x} by {array_y} INITIAL RAW.png',init_array,cmap='jet')
# resize_init_array=upsize(init_array,array_y,array_x,upsize_val)
# plt.imsave(f'{idt} {title} {array_x} by {array_y} INITIAL RESIZE.png',resize_init_array,cmap='jet')
# resize_init_array=upsize(init_array,array_y,array_x,upsize_val)
# plt.savefig(f'{idt} {title} {array_x} by {array_y} INITIAL.png')
# plt.show()



arr=fin_diff(init_array,array_y,array_x,bound_pot)
i=1
t_a=time.time()



while True:
    chk,converge_val_n=converge_check(arr,array_y,array_x,converge_val)
    a=int(abs(np.log10(converge_val)))
    print(f'Step {i}, {100*round(converge_val/converge_val_n,a):.6f}% Convergence')
    if chk==True:
        t=round(time.time()-t_a, 10)
        v_plot(arr,array_y,array_x,mid_point)
        plt.savefig(f"{idt} {title} POTENTIAL GRAPH {array_y} {array_x}.png")
        plt.show()
        plt.imshow(arr,cmap=GLOBAL_CMAP,origin='lower')
        plt.colorbar()
        plt.title(f'{title}\n{array_y}x{array_x} grid, {bound_pot}V\nConvergence of {converge_val} found in {i} steps\n{t}s',loc='left',fontsize=7)
        plt.imsave(f'{idt} {title} {array_x} by {array_y} {i} steps {bound_pot}V {converge_val} {t}s RAW.png',arr,cmap='jet')
        resize_arr=upsize(arr,array_y,array_x,upsize_val)
        plt.imsave(f'{idt} {title} {array_x} by {array_y} {i} steps {bound_pot}V {converge_val} {t}s RESIZE.png',resize_arr,cmap=GLOBAL_CMAP)
        plt.savefig(f'{idt} {title} {array_x} by {array_y} {i} steps {bound_pot}V {converge_val} {t}s.png')
        plt.show()
        print(f'\n{t}s\n{idt}')
        break
    else:
        arr=fin_diff(arr,array_y,array_x,bound_pot)
    i+=1
