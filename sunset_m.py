import  numpy as np
import scipy as sp
from scipy import optimize
import matplotlib.pyplot as plt
import pandas as pd

import time
import decimal
from decimal import Decimal

decimal.getcontext().prec = 256

lb=np.linspace(10, 3990, num=200, endpoint = True)
h=6.62607015e-34
c=299792458
k=1.3806488e-23
Ts=5.75e3
Tbt=3e3
Tbl=22e3
Tac=28e3
n=1.00029
N=2.504e25
Rt= 6.378e6
Sz=8000
So=np.sqrt((Rt+Sz)**2-Rt**2)
Angoli=np.linspace(0, 90, num=91, endpoint = True)

def Beta():
    B=(8*np.pi**3)/(3*(lb*10**(-9))**4*N)*(n**2-1)**2
    return B

def Beta1(L):
    B=(8*np.pi**3)/(3*(L*10**(-9))**4*N)*(n**2-1)**2
    return B

def St(theta):
    theta=theta*np.pi/180
    return np.sqrt((Rt*np.cos(theta))**2+2*Rt*Sz+Sz**2)-Rt*np.cos(theta)

def DEN(L,T):
    E=np.expm1((h*c)/(L*k*T*10**(-9)))
    D=(2*c)/((L*10**(-9))**4*(E))
    return D

def Phot(T):
    xr=np.random.uniform(10, 3990, 100000)
    yr=np.random.random(100000)
    hitmask=yr<(DEN(xr,T)/DEN(2898/T*1000,T))

    n=np.histogram(xr[hitmask], bins=200, range=[0,4000])
           
    return np.array(n[0], dtype=int)

def Obs(S, D):
    B=Beta()
    ex=np.exp(-B*S)
    res=D*ex
    return res

def zenith(P):
    
    F=Obs(Sz, P)
    return F

def orizzonte(P):

    F=Obs(So, P)
    return F

def angolato(P, theta):
    Sa=St(theta)
    return Obs(Sa, P)

def graph_fl(N):

    mcint=np.empty(0)

    for t in Angoli:
        mcint=np.append(mcint, [(4000*10**(-9))/100000*np.sum(N*np.exp(-Beta()*St(t)))/((4000*10**(-9))/100000*np.sum(N*np.exp(-Beta()*Sz)))])
        
    plt.plot(Angoli, mcint, 'o-')
    plt.xlabel("Angolo rispetto allo zenith ( °)")
    plt.ylabel("Flusso relativo di fotoni per unità di superficie")
    plt.title("Flusso relativo integrato di fotoni in dipendenza dell'angolo di incidenza")
    plt.grid(visible=True)
    plt.show()

def graph(na, ze, oz, st):
    
    plt.bar(lb, height=na, width=20, color='dodgerblue', alpha=0.6, label='Nessun assorbimento', linewidth=0.5, edgecolor='dimgray')
    plt.bar(lb, height=ze, width=20, color='orange', alpha=0.6, label='Zenith', linewidth=0.5, edgecolor='dimgray')
    plt.bar(lb, height=oz, width=20, color='crimson', alpha=0.6, label='Tramonto', linewidth=0.5, edgecolor='dimgray')

    plt.title('Grafico della distribuzione di fotoni per unità di superficie della stella '+ st)
    plt.legend()
    plt.xlabel("Lunghezza d'onda (nm)")
    plt.ylabel("Numero di fotoni per unità di superficie")
    plt.show()

def graph2(L, N):
    plt.bar(L, height=N, width=20, color='crimson', alpha=0.6, label='Stella X, 28°', linewidth=0.5, edgecolor='dimgray')
    plt.title('Grafico della distribuzione di fotoni per unità di superficie della stella sconosciuta')
    plt.legend()
    plt.xlabel("Lunghezza d'onda (nm)")
    plt.ylabel("Numero di fotoni per unità di superficie")
    plt.show()

def exp_large_float(x):
    d=[Decimal(i) for i in x]
    D=np.array(d)
    return np.exp(D)

def Temp(L, T, A):

    #Per questa funzione è necessario passare ai Decimal durante lo svolgimento per via delle dimensioni che l'esponenziale può arrivare ad avere durante i tentativi di calcolo del fit

    B=Beta1(L)
    b=[Decimal(i) for i in B]
    B=np.array(b)
    
    l=[Decimal(i) for i in L]
    L=np.array(l)*Decimal(10**(-9))
    
    h=Decimal(6.62607015e-34)
    c=Decimal(299792458)
    k=Decimal(1.3806488e-23)
    
    S=Decimal(St(28))
    
    ebs=np.exp(-B*S)
    gamma=(h*c)/k

    T=Decimal(T)

    A=Decimal(A)
    
    pot=gamma/(L*T)
    E=exp_large_float(pot)-1
    N=np.array(A*2*c*ebs/(L**4*E), dtype=np.float64)
    T=float(T)
    A=float(A)
    return N

def find_T(L, N):

    timestart=time.time()
    pnames=['T', 'A']
    params, params_covariance = optimize.curve_fit(Temp, L, N, p0=[4e3, 1e-25], maxfev=1000, bounds=((1, 3e-30),(3e4, 1e-8)))
    timestop=time.time()
    Tem=params[0]
    Temerr=params_covariance.diagonal()[0]
    return (Tem, Temerr)

def star(choice):
    match choice:
        case 'sole':
            
            print("Hai selezionato il Sole.")

            ph=Phot(Ts)
            
            na=ph
            ze=zenith(ph)
            oz=orizzonte(ph)

            graph_fl(ph)
            
            graph(na, ze, oz, 'Sole.')

            return 0
        
        case 'bet':
            
            print("Hai selezionato Betelgeuse.")

            ph=Phot(Tbt)
            
            na=ph
            ze=zenith(ph)
            oz=orizzonte(ph)
               
            graph_fl(ph)
            
            graph(na, ze, oz, 'Betelgeuse.')

            return 0

        case 'bel':
            
            print("Hai selezionato Bellatrix.")

            ph=Phot(Tbl)
            
            na=ph
            ze=zenith(ph)
            oz=orizzonte(ph)

            graph_fl(ph)
            
            graph(na, ze, oz, 'Bellatrix.')
            
            return 0
        
        case 'alpha':
            
            print("Hai selezionato Alpha Crucis.")

            ph=Phot(Tac)
            
            na=ph
            ze=zenith(ph)
            oz=orizzonte(ph)
               
            graph_fl(ph)
            
            graph(na, ze, oz, 'Alpha Crucis.')
            
            return 0
        
        case 'q':
            
            return 2
        
        case 'm':
            
            return 1
        
        case _:
            
            print("Errore, per favore inserire una delle opzioni indicate.")
            
            return 0
