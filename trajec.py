#code for simulation with a jump kernel using levy distribution and produces the file for entry and exit times into coalascent zone. We are using a lattice and spherical pbc. 
#file format entrytime[space]exittime
from scipy.stats import levy_stable
import random
import math
import re
def jump(x,y):#takes an input of coordinate and gives a new coordinate after jump
    a1=levy_stable.rvs(2.0, 1.0,loc=5, scale=1, size=1 )#rvs(alpha,beta,loc,scale,size). aplha gives the power of the decay of tail, beta gives the skew(with beta zero the distribution is symmetric,with beta=1.0 we wont have any negative values and curve is skewed to right.loc  
    jumplength=a1[0]#scale allows us to scale the distribution and size gives us number of elements we are drawing from the levy stable distribution. loc is set at 5 as for small jump sizes one has to deal with assymetryties because of forcing everything on the lattice.
    a2=random.random()*2*math.pi#uniformly chooses an angle between [0,2pi). Uses mt(19337) pseudo random generator
    #We want the 4 corners of the lattice at the center of which the jumping point falls. Then we choose which of the 4 available lattice sites the particle jumps to by finding the probabilities of jump to each of the sites by calculating r**-alpha for all 4 sites. r is the distance of the site
    #from the original jumping point.
    x0=int(jumplength*math.cos(a2))#first lattice point
    y0=int(jumplength*math.sin(a2))
    r=math.sqrt(x0**2+y0**2)
    if r>0:
        prob0=r**(-2)
    else:
        prob0=0    
    x1=int(jumplength*math.cos(a2))#second lattice point
    y1=int(jumplength*math.sin(a2))+1
    r=math.sqrt(x1**2+y1**2)
    if r>0:
        prob1=r**(-2)
    else:
        prob1=0 
    x2=int(jumplength*math.cos(a2))+1#third lattice point
    y2=int(jumplength*math.sin(a2))
    r=math.sqrt(x2**2+y2**2)
    if r>0:
        prob2=r**(-2)
    else:
        prob2=0 
    x3=int(jumplength*math.cos(a2))#fourth lattice point
    y3=int(jumplength*math.sin(a2))
    r=math.sqrt(x3**2+y3**2)
    if r>0:
        prob3=r**(-2)
    else:
        prob3=0 
    sumprob=prob0+prob1+prob2+prob3
    
    #choose the lattice point where the jump ends
    a3=random.random()
    
    if a3<(prob0/sumprob):
        x=x+x0
        y=y+y0
    elif a3<((prob0+prob1)/sumprob):
        x=x+x1
        y=y+y1
    elif a3<((prob0+prob1+prob2)/sumprob):
        x=x+x2
        y=y+y2
    else:
        x=x+x3
        y=y+y3
        
    #check for pbc. Circle of radius 10000000
    r=math.sqrt(x**2+y**2)
    if r<10000000:
        x=x 
        y=y
    else:
        r1=10000000-r%10000000 #if the particle is x outside of the circle then under pbc it will be x from the boundary and reappear at the diametrically opposite point
        x=int(r1*math.cos(a2+math.pi)) 
        y=int(r1*math.sin(a2+math.pi))   
        
    return x,y
outr=open('out.txt','w')    
for replicates in range(5000):
    coalascetime=[0,0]#keeps tabs of entry and exit times
    xs=1.0 #(starting point)
    ys=0.0 #(starting point)
    rs=xs**2+ys**2
    coal_state=0#gives whether particle is in coalascent zone. We are using a circle of radius 10 as our coalascent zone
    if rs<100:
        coal_state=1
    else:
        coal_state=0
        
        
    for times in range(100000):
        x,y=jump(xs,ys)
        r=x**2+y**2
        if r<100:
            if coal_state==0:
                coalascetime[0]=times
                coal_state=1
            if coal_state==1:
                coal_state=1    
   
        if r>100:
            if coal_state==0:
                coal_state=0
            if coal_state==1:
                coalascetime[1]=times
                coal_state=0
                outr.write('%f %s\n'%(coalascetime[0],coalascetime[1]))
        
        xs=x
        ys=y
                        
outr.close()            
            
                                                                                                          
