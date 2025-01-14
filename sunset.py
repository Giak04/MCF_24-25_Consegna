import numpy as np
import scipy as sp
from scipy import optimize
import matplotlib.pyplot as plt
import pandas as pd
import sunset_m as sns

h=6.62607015e-34
c=299792458
k=1.3806488e-23
Ts=5.75e3

lb=np.linspace(10, 3990, num=200, endpoint = True)
Ts=5.75e3
Tbt=3e3
Tbl=22e3
Tac=28e3

on=1

while on == 1:

    cat=input("Digitare '1' se si desidera studiare la distribuzione di fotoni di una stella nota.\nDigitare '2'se si desidera determinare la temperatura di una stella ignota data la distribuzione di fotoni. \nDigitare '0' se si desidera terminare l'esecuzione.\n")
    
    if cat == '0':
        on=0
        
    elif cat == '1':
        print("\nStudio della distribuzione di fotoni di una stella nota.")
        control=0
        
        while control == 0:
            print("Selezionare quale stella usare per lo studio:")
            choice=input("Digitare 'sole' per studiare il Sole.\nDigitare 'bet' per studiare Betelgeuse.\nDigitare 'bel' per studiare Bellatrix.\nDigitare 'alpha' per studiare Alpha Crucis.\nDigitare 'm' per tornare al menù di selezione iniziale.\nDigitare 'q' per terminare il programma.\n")
            
            control=sns.star(choice)
                 
            if control == 2:
                on=0
               
    elif cat == '2':
        
        print('Ricavare la temperatura di una stella ignota a partire dalla distribuzione di fotoni')
        
        print('Lettura dei dati in corso...')
        
        Tab=pd.read_csv('observed_starX.csv')
        
        print('Dati acquisiti.\n', Tab, "\n L'angolo di inclinazione rispetto allo zenith è di 28°")

        la=np.array(Tab['lambda (nm)'])
        Nobs=np.array(Tab['photons'])

        sns.graph2(la,Nobs)
       
        print('Calcolo della temperatura in corso')
        
        Tem, Temerr = sns.find_T(la, Nobs)
        
        print('La temperatura della stella si stima essere di ', Tem,'+-', Temerr, ' K.')

    else:
        
        print('\nErrore. Si prega di inserire una risposta tra quelle elencate.\n')

print("Programma terminato. Arrivederci.")
