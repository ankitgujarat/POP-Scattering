#!bin/python
import numpy as np
import matplotlib.pyplot as plt

def plot(t):
    ir_mode = [0.0001,0.0000,0.0000,0.0000,0.0000,0.0000,0.3383,11.2142,0.0000,0.0000,0.7520,16.9140,1.8258,0.0000,0.0000,0.0000,4.4220,0.0000,21.5583,0.0000,0.0000,28.9572,0.0000,18.3179,0.0000,14.6069,0.0000,3.2368,0.0000,32.0120]    
    x = np.linspace(1,30,30)
    plt.close('all')
    
    f,(ax1,ax2) = plt.subplots(2,sharex=True)
    #print x
    for i in range(0,len(t)):
        ax1.bar(x[i],t[i],align='center',alpha=0.5)
        ax2.bar(x[i],ir_mode[i],align='center',alpha=0.5)
    ax1.set_title('Dipole Strength & IR activity for each Mode')
    ax2.set_xlim([0,32])
    ax1.grid()
    ax2.grid()
    plt.show()
    return 0

def func1(x):
    l=0
    t=""
    temp=[]
    disp=[[],[],[]]
    new = []
    x=x.strip()
    while l<len(x):
        if x[l]!=chr(32) and x[l]!=chr(40) and x[l]!=chr(41):
            t+=x[l]
        if x[l]==chr(32):
            temp.append(l)
        l=l+1
    for i in range (0,len(temp)-1):
        if temp[i+1]-temp[i]>1:
            new.append(temp[i])
            new.append(temp[i+1])
    disp[0].append(float(x[new[0]:new[1]]))
    disp[1].append(float(x[new[4]:new[5]]))
    disp[2].append(float(x[new[8]:new[9]]))
    return disp
            

#The above function first removes all the spaces and brackets and appends the entire line into a string
#Then the location of the spaces is calculated and then substracted to find the length of the actual string
#Depending on the presence of the negative sign and the magnitude, the string lengths are varied
#The string is updated every time a value is entered into the displacement matrix

def dielectric(r):
    f=[]
    g=0
    diele=[[],[],[]]
    x=r[0].strip()
    while g < len(x):
        if x[g]==chr(32):
            f.append(g)
        g=g+1
    val1 = fin_max_diff(f)
    for i in range (0,3):
        diele[i].append(float(r[i].strip()[0:f[0]]))
        diele[i].append(float(r[i].strip()[int(val1[0]):int(val1[1])]))
        diele[i].append(float(r[i].strip()[f[len(f)-1]:]))
    for i in range (0,3):
        print diele[i][0],diele[i][1],diele[i][2]
        

def fin_max_diff(f):
    k=[]
    for i in range(0,len(f)-1):
        if f[i+1]-f[i]>=2:
            k.append(f[i]+1)
            k.append(f[i+1]+1)
    return k

def eff_chg1(q,atm):
    temp_q=""
    atoms = []
    val1 = []
    n=0
    zeff = []
    temp=[[],[],[]]
    for i in range (0,len(q)):
        temp_q=q[i].strip()
        if temp_q[0]==chr(97):
            atoms.append(i)
    
    for j in range (0,len(atoms)):
        q.pop(atoms[j]-j)
    
    for i in range (0,len(q)):
        w = []
        temp_q = q[i].strip()
        
        for j in range(0,len(temp_q)):
            if temp_q[j]==chr(32):
                w.append(j)
        val1 = fin_max_diff(w)
        temp[n].append(float(temp_q[0:w[0]]))
        temp[n].append(float(temp_q[val1[0]:val1[1]]))
        temp[n].append(float(temp_q[w[len(w)-1]:]))
        n = n + 1
        if n > 2:
            n = 0
            zeff.append(temp)
            temp = [[],[],[]]
    return zeff
    

def ir_strength(dis,chg):
    dist=dis
    stgt = []
    x_dir=[]
    y_dir=[]
    z_dir=[]
    for mode in range (0,30):
        ir=np.array([[0],[0],[0]])
        temp=[]
        x=[]
        y=[]
        for j in range (0,10):
            temp.append(dist[j])
        for i in range (0,10):
            x=np.array(temp[i])
            y=np.array(chg[i])
            ir=ir+y.dot(x)
        x_dir.append(ir[0])
        y_dir.append(ir[1])
        z_dir.append(ir[2]) 
        stgt.append(np.linalg.norm(ir))
        for k in range (0,10):
            dist.pop(0) 
    plot(x_dir)
    plot(y_dir)
    plot(z_dir)
    return stgt
        
        

disp_pat=open("G:\Charge Transport Theory\Calculations\Ga2O3\K_D\dynmat.out","r")
eff_ch=open("G:\Charge Transport Theory\Calculations\Ga2O3\K_D\Ga2O3.dyn","r")
a=[]
b=[]
c=[]
d=[]
z=[]
p=[]
dpat=[]
dpat1 = {}
j=0
w=0

#**************************************************************************************************************************************
#Code for Dielectric Tensor and Effective Charge
#**************************************************************************************************************************************
line=input("Enter the line number where the word Dielectric Tensor is mentioned")
eff=input("Enter the line number of the word effective charge")
atm=input("Enter the number of atoms per unit cell")
for i in eff_ch:
    d.append(i)
dielectric(d[line+1:line+5])
effective_charge=eff_chg1(d[eff+1:eff+41],atm)
print "********************************************************************"
for p in range(0,10):
    for q in range(0,3):
        print effective_charge[p][q][0],effective_charge[p][q][1],effective_charge[p][q][2]



#***********************************************************************************************************************************
#Code for displacement pattern
#***********************************************************************************************************************************
for i in disp_pat:
    a.append(len(i))
    c.append(i)
#print a
#The above code reads lines from the files and stores it in an array "c". It also stores the length of each line
#print "**************************************************"
val=max(a)
#print val
while j<len(a):
    if a[j]==val:
        b.append(j)
    j=j+1
#print b
#The above code finds the index of the lines with max length and stores it in an array "b"
#print "***************************************************"
#Testing for the IR strength of the first mode... The code will later be modified for other modes too
count = 0
while w<len(b):
    x=c[b[w]]
    displacement=func1(x)
    dpat.append(displacement)
    w=w+1
print dpat

#s = ir_strength(dpat,effective_charge)
#plot(s)


#The above code picks lines from the array c with the max length and calls the function