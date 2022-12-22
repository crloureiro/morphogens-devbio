from __future__ import print_function
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
from IPython.core.display import HTML
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import numpy as np
from math import *
from Diffusion1D import *

# Return the probability of gene activation
def activation(b,c,n,K): 
    return b*c**n/(c**n+K**n)

# Return the
def limit_tissue(c0,k,D,s1,s2):
    l1=-log(s1/c0)*sqrt(D/k)
    l2=-log(s2/c0)*sqrt(D/k)
    return l1,l2

def Cmorphogen_at_activation(p,K,n):
    return K/(((1/p)-1)**(1./n))

def Cmorphogen(c0,x,k,D):
    l=sqrt(D/k)
    return c0*exp(-x/l)

def limit_tissue(c0,k,D,s1,s2):
    l1=-log(s1/c0)*sqrt(D/k)
    l2=-log(s2/c0)*sqrt(D/k)
    return l1,l2

def plot_french_flag(c0=(5,10,0.01), D=(0.1,100), k=(0.1,1,0.01)):
    x=np.arange(0,50,0.5)
    c= [Cmorphogen(c0,x[i],k,D) for i in range(len(x))]
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    L1,L2=limit_tissue(c0,k,D,4,0.5)
    rect = plt.Rectangle((0, -1.5), L1, 1.5, color='blue', alpha=0.8)
    rect2 = plt.Rectangle((-12, -1.5),12, 1.5, color='black', alpha=0.5)
    rect3 = plt.Rectangle((L1, -1.5), L2-L1, 1.5, color='gray', alpha=0.1)
    rect4 = plt.Rectangle((L2, -1.5), 50-L2, 1.5, color='red', alpha=0.8)
    plt.axes(ax).spines['left'].set_position(('data',0))
    #ax.set_xlim(10,30)
    plt.text(-9,-1.2,'Source \n cells', color='white', fontsize=12)
    plt.text(L1/2,-1,'A', color='white', fontsize=14)
    plt.text((L2-L1)/2+L1,-1,'B', color='black', fontsize=14)
    plt.text(L2+(50-L2)/2,-1,'C', color='white', fontsize=14)
    plt.text(-12,4,'Threshold 1', color='blue', fontsize=12)
    plt.text(-12,0.5,'Threshold 2', color='red', fontsize=12)
    ax.add_patch(rect)
    ax.add_patch(rect2)
    ax.add_patch(rect3)
    ax.add_patch(rect4)
    plt.plot([L1, L1], [-1, 4], lw=1, color='blue', )
    plt.plot([0, L1], [4, 4], lw=1, color='blue', )
    plt.plot([L2, L2], [-1, 0.5], lw=1, color='red')
    plt.plot([0, L2], [0.5, 0.5], lw=1, color='red', )
    plt.xlabel('Distance to the source x', fontsize=12)
    plt.ylim((-1.5,10))
    plt.xlim((-12,50))
    #plt.ylabel('Morphogen concentration', fontsize=12)
    plt.plot(x, c, label='Morphogen concentration')
    plt.legend()
    ax.grid(color='gray',alpha=0.3)

# Return the probability of gene activation
def activation(b,c,n,K): 
    return b*c**n/(c**n+K**n)

# Return the 
def limit_tissue(c0,k,D,s1,s2): 
    l1=-log(s1/c0)*sqrt(D/k)
    l2=-log(s2/c0)*sqrt(D/k)
    return l1,l2

# Return the
def Cmorphogen_at_activation(p,K,n):
    return K/(((1/p)-1)**(1./n))

# Plot morphogen distribution
def plot_activation(K1=(0.1,0.5,0.01),K2=(0.01,0.1,0.01)):
    c0=1 # Morphogen concentration at the source
    k=0.5 # Morphogen degradation rate
    D=20 # Morphogen diffusion
    p=0.9
    x=np.arange(0,30,0.5) # Set the list of distances x
    c=[Cmorphogen(c0,x[i],k,D) for i in range(len(x))] # Calculate concentrations of morphogen
    e=[activation(1,c[i],4,K1) for i in range(len(x))] # Calculate concentrations of effector A
    f=[activation(1,c[i],4,K2) for i in range(len(x))] # Calculate concentrations of effector B
    
    C1,C2=Cmorphogen_at_activation(p,K1,4),Cmorphogen_at_activation(p,K2,4)
    L1,L2=limit_tissue(c0,k,D,C1,C2) # Determine the distance x corresponding to activation thresholds
    
    # Create the figure
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    # Plot the concentrations of morphogens and effectors
    plt.plot(x, c, label='Morphogen concentration (C)')
    plt.plot(x, e, 'g--', label='Probability of effector A activation', color='blue')
    plt.plot(x, f, 'g--',label='Probability of effector B activation', color='dimgray')
    plt.plot([0, 30], [p, p], 'g--', lw=1, color='black' )
    plt.text(-10,p,'Activation condition p', color='black', fontsize=12) #Threshold 2, c=K2
    
    # Draw Source cell domain
    rect = plt.Rectangle((-9, -0.1),9, 0.1, color='black', alpha=0.5) # Source cells
    plt.text(-8,-0.08,'Source cells', color='white', fontsize=14) # Source cells
    ax.add_patch(rect)
    
    if L1<0 and L2>0: # Morphogen concentration below activation threshold for effector A
        
        plt.plot([0, L2], [C2, C2], lw=1, color='dimgray' )
        plt.plot([L2, L2], [0, C2], lw=1, color='dimgray')
        plt.plot([L2, L2], [0, 1.2], 'g--', lw=1, color='dimgray' )
        rect3 = plt.Rectangle((0, -0.1), L2, 0.1, color='gray', alpha=0.1) # Receiving cells, only B activated
        rect4 = plt.Rectangle((L2, -0.1), 30-L2, 0.1, color='red', alpha=0.8) # Receiving cells, none activated
        ax.add_patch(rect3)
        ax.add_patch(rect4)
        plt.text((L2-L1)/2+L1,-0.08,'B', color='black', fontsize=14) # domain B
        plt.text(L2+(30-L2)/2,-0.08,'C', color='white', fontsize=14) # domain C
        plt.text(-8,C2,'C threshold for\nB activation', color='dimgray', fontsize=12) # Threshold 2, c=K2
        plt.text(10,0.6,'Warning, initial morphogen\nconcentration below A\nactivation threshold!', color='black', fontsize=14) 
    
    elif L2<0: # Morphogen concentration below activation thresold for effector A and B

        rect4 = plt.Rectangle((0, -0.1), 30, 0.1, color='red', alpha=0.8) # Receiving cells, none activated
        ax.add_patch(rect4)
        plt.text(L2+(30-L2)/2,-0.08,'C', color='white', fontsize=14) # domain C
        plt.text(10,0.6,'Warning, initial morphogen concentration \n below A&B activation threshold!', color='black', fontsize=14) 
       
    else:
    # Plot the morphogen concentration corresponding to the activation threshold (probability of activation condition)
        plt.plot([L1, L1], [0, C1], lw=1, color='blue', )
        plt.plot([0, L1], [C1, C1], lw=1, color='blue', )
        plt.plot([L2, L2], [0, C2], lw=1, color='dimgray')
        plt.plot([0, L2], [C2, C2], lw=1, color='dimgray', )
        plt.plot([L1, L1], [0, 1.2], 'g--', lw=1, color='blue' )
        plt.plot([L2, L2], [0, 1.2], 'g--', lw=1, color='dimgray' )
    
    # Draw the different tissue domain A,B,C
        rect2 = plt.Rectangle((0, -0.1), L1, 0.1, color='blue', alpha=0.8) # Receiving cells, A and B activated
        rect3 = plt.Rectangle((L1, -0.1), L2-L1, 0.1, color='gray', alpha=0.1) # Receiving cells, only B activated
        rect4 = plt.Rectangle((L2, -0.1), 30-L2, 0.1, color='red', alpha=0.8) # Receiving cells, none activated
        ax.add_patch(rect2)
        ax.add_patch(rect3)
        ax.add_patch(rect4)
    
    # Annotations of the graph
        plt.text(L1/2,-0.08,'A', color='white', fontsize=14) # domain A
        plt.text((L2-L1)/2+L1,-0.08,'B', color='black', fontsize=14) # domain B
        plt.text(L2+(30-L2)/2,-0.08,'C', color='white', fontsize=14) # domain C
        plt.text(-6.5,C1-0.02,'Activation\nthreshold A', color='blue', fontsize=12) # Threshold 1, c=K1
        plt.text(-6.5,C2-0.02,'Activation\nthreshold B', color='dimgray', fontsize=12) # Threshold 2, c=K2
    
    # Configuration of the axis
    plt.axes(ax).spines['left'].set_position(('data',0))
    plt.ylim((-0.1,1.2))
    plt.xlim((-10,30))
    plt.xlabel('Distance to the source x', fontsize=12)
    plt.legend()
    plt.show()

    
# Return the probability of receptor activation
def activation_receptor(b,c,n,K,r):
    return r*c**n/(c**n+K**n)


# Plot competence
def plot_activation_competence(Rlim=(1,40,1)):
    c0=1
    k=0.2
    D=50
    K1=0.3
    K2=0.14
    p=0.9
    xmax=40
    
    x=np.arange(0,xmax,0.5) # Set the list of distances x
    r1=[1 for i in range (0,2*Rlim)]
    r2=[0 for i in range (2*Rlim,len(x))]
    R=r1+r2
    c= [Cmorphogen(c0,x[i],k,D) for i in range(len(x))] # Calculate concentrations of morphogen
    e=[activation_receptor(1,c[i],4,K1,R[i]) for i in range(len(x))] # Calculate concentrations of effector A
    f=[activation_receptor(1,c[i],4,K2,R[i]) for i in range(len(x))] # Calculate concentrations of effector B
    e2=[activation_receptor(1,c[i],4,K1,1) for i in range(len(x))] # Calculate concentrations of effector A
    f2=[activation_receptor(1,c[i],4,K2,1) for i in range(len(x))] # Calculate concentrations of effector B
    C1,C2=Cmorphogen_at_activation(p,K1,4),Cmorphogen_at_activation(p,K2,4)
    L1,L2=limit_tissue(c0,k,D,C1,C2) # Determine the distance x corresponding to activation thresholds
    
        
    # Create the figure
    fig, (ax2, ax1) = plt.subplots(1, 2, figsize=(18,6), gridspec_kw={'hspace': 0, 'wspace': 0.1})
        
    # Plot the morphogen concentration and effector response for competent cells#
    
    ax2.set_xlabel("Distance X")
    ax2.set_title("Equal cell competence", fontsize=14)
    
    ax2.plot(x, c, label='Morphogen concentration')
    ax2.plot(x, e2, 'g--', label='Probability of effector A activation', color='blue')
    ax2.plot(x, f2, 'g--', label='Probability of effector B activation', color='dimgray')
    
    rectb = plt.Rectangle((-9, -0.1),9, 0.1, color='black', alpha=0.5) # Source cells
    rect4 = plt.Rectangle((0, -0.1), L1, 0.1, color='blue', alpha=0.8) # Receiving cells, A and B activated
    rect5 = plt.Rectangle((L1, -0.1), L2-L1, 0.1, color='gray', alpha=0.1) # Receiving cells, only B activated
    rect6 =plt.Rectangle((L2, -0.1), xmax-L2, 0.1, color='red', alpha=0.8)
    ax2.add_patch(rectb)
    ax2.add_patch(rect4)
    ax2.add_patch(rect5)
    ax2.add_patch(rect6)
    ax2.plot([L1, L1], [0, C1], lw=1, color='blue', )
    ax2.plot([0, L1], [C1, C1], lw=1, color='blue', )
    ax2.plot([L2, L2], [0, C2], lw=1, color='dimgray', )
    ax2.plot([0, L2], [C2, C2], lw=1, color='dimgray', )
    ax2.text(-9,-0.08,'Source cells', color='white', fontsize=12) # source cells
    ax2.text(L1/2,-0.08,'A', color='white', fontsize=14) # domain A
    ax2.text((xmax-L2)/2+L2,-0.08,'C', color='black', fontsize=14) # domain C
    ax2.text((L2-L1)/2+L1,-0.08,'B', color='black', fontsize=14) # domain B
    ax2.text(-10,C1-0.02,'Activation\nthreshold A', color='blue', fontsize=12) # Threshold 1
    ax2.text(-10,C2-0.02,'Activation\nthreshold B', color='dimgray', fontsize=12)

    
    # Plot the morphogen concentration and effector response for inequally competent cells##
    rect = plt.Rectangle((-9, -0.1),9, 0.1, color='black', alpha=0.5) #Source cells
    ax1.add_patch(rect)
    ax1.plot([L1, L1], [0, C1], lw=1, color='blue', )
    ax1.plot([0, L1], [C1, C1], lw=1, color='blue', )
    ax1.plot([L2, L2], [0, C2], lw=1, color='dimgray', )
    ax1.plot([0, L2], [C2, C2], lw=1, color='dimgray', )
    ax1.text(-10,C1-0.02,'Activation\nthreshold A', color='blue', fontsize=12) #Threshold 1
    ax1.text(-10,C2-0.02,'Activation\nthreshold B', color='dimgray', fontsize=12)
    if L2>Rlim and L1<Rlim:
        L2=Rlim
        rect2 = plt.Rectangle((0, -0.1), L1, 0.1, color='blue', alpha=0.8) #A
        rectc =plt.Rectangle((L2, -0.1), xmax-L2, 0.1, color='red', alpha=0.8)
        rect3 = plt.Rectangle((L1, -0.1), L2-L1, 0.1, color='gray', alpha=0.1)
        ax1.text(L1/2,-0.08,'A', color='white', fontsize=14) # domain A
        ax1.text((L2-L1)/2+L1,-0.08,'B', color='black', fontsize=14) # domain B
        ax1.text((xmax-L2)/2+L2,-0.08,'C', color='black', fontsize=14) # domain C
        ax1.add_patch(rect2)
        ax1.add_patch(rect3)
        ax1.add_patch(rectc)
    elif L1>Rlim:
        L1=Rlim
        rectc =plt.Rectangle((L1, -0.1), xmax-L1, 0.1, color='red', alpha=0.8) # C
        rect2 = plt.Rectangle((0, -0.1), L1, 0.1, color='blue', alpha=0.8) # A
        ax1.text((xmax-L1)/2+L1,-0.08,'C', color='black', fontsize=14) # domain C
        ax1.text(L1/2,-0.08,'A', color='white', fontsize=14) # domain A
        ax1.add_patch(rect2)
        ax1.add_patch(rectc)
    else:
        rect = plt.Rectangle((-9, -0.1),9, 0.1, color='black', alpha=0.5) #Source cells
        rect2 = plt.Rectangle((0, -0.1), L1, 0.1, color='blue', alpha=0.8) #Receiving cells, A and B activated
        rect3 = plt.Rectangle((L1, -0.1), L2-L1, 0.1, color='gray', alpha=0.1) #Receiving cells, only B activated
        rectc =plt.Rectangle((L2, -0.1), xmax-L2, 0.1, color='red', alpha=0.8)
        ax1.text(L1/2,-0.08,'A', color='white', fontsize=14) #domain A
        ax1.text((L2-L1)/2+L1,-0.08,'B', color='black', fontsize=14) #domain B
        ax1.text((xmax-L2)/2+L2,-0.08,'C', color='black', fontsize=14) #domain C
        ax1.add_patch(rect2)
        ax1.add_patch(rect3)
        ax1.add_patch(rectc)
    
    ax1.set_xlabel("Distance X")    
    ax1.set_title("Differential cell competence (ON/OFF activation switch)",fontsize=14)
    
    ax1.plot(x, c, label='Morphogen concentration')
    ax1.plot(x, e, 'g--', label='Probability of effector A activation', color='blue')
    ax1.plot(x, f, 'g--', label='Probability of effector B activation', color='dimgray')

    #ax1.plot([L1, L1], [0, C1], lw=1, color='blue', )
    #ax1.plot([0, L1], [C1, C1], lw=1, color='blue', )
    #plt.text(L1/2,-0.08,'A', color='white', fontsize=14) #domain A
    ax1.text(-9,-0.08,'Source cells', color='white', fontsize=14) #source cells
    #ax1.text(-9,C1-0.02,'Activation\nthreshold A', color='blue', fontsize=12) #Threshold 1, c=K1
    ax1.plot([Rlim, Rlim], [0, 1.4], color='black', lw=2 )
    ax1.text(Rlim/2-4,1.05,'Competent\n   cells', color='black', fontsize=14)
    ax1.text((xmax-Rlim)/2-4+Rlim,1.05,'Uncompetent\n     cells', color='black', fontsize=14)
    
    # Configuration of the axis
    plt.axes(ax1).spines['left'].set_position(('data',0))
    plt.axes(ax2).spines['left'].set_position(('data',0))
    ax1.set_ylim((-0.1,1.4))
    ax1.set_xlim((-9.5,xmax))
    ax2.set_ylim((-0.1,1.4))
    ax2.set_xlim((-9.5,xmax))
    #ax1.legend()
    ax2.legend()
    plt.show()

    
# Construct the gradient
def construct_gradient(C0=(1,10),D=(0.1,1)):   
    N=100
    U=np.zeros(N)
    S=np.zeros(N)
    x=np.arange(N)*10.0/N
    coef = [[1,D]] # First coeff is initial concentration, second one diffusion vector
    t=0
    diffusion = Diffusion1D(N)
    diffusion.config(0.0001,coef,"dirichlet",C0,"dirichlet",0,S)
    [U1,t]=diffusion.iterations(U,t,0.001)
    diffusion.config(0.001,coef,"dirichlet",C0,"dirichlet",0,S)
    [U2,t]=diffusion.iterations(U1,t,0.01)
    diffusion.config(0.01,coef,"dirichlet",C0,"dirichlet",0,S)
    [U3,t]=diffusion.iterations(U2,t,0.1)
    diffusion.config(0.1,coef,"dirichlet",C0,"dirichlet",0,S)
    [U4,t]=diffusion.iterations(U3,t,1)
    figure(figsize=(8,6))
    plot(x,U1,label="t=0.001")
    plot(x,U2,label="t=0.01")
    plot(x,U3,label="t=0.1")
    plot(x,U4,label="t=1")
    legend(loc="upper right")
    xlabel("Distance from the source")
    ylabel("Concentration of the Morphogen")
    axis([0,10,0,10])
    grid()