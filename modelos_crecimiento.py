#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Exportamos todo lo necesario
from flask import Flask, render_template, request, redirect, jsonify
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import io
import base64


# In[9]:


app = Flask(__name__)

#Damos dos opciones al usurio, una si conoce directamente el volumen del tumor y otra si solo conoce largo y ancho
@app.route('/')
def index():
    return render_template('index2.html')

#Procesamos la respuesta
@app.route('/opcion', methods=['POST'])
def pregunta():
    respuesta = request.form['opcion']
    if respuesta.lower() == 'si':
        return render_template('opcion_si.html')
    elif respuesta.lower() == 'no':
        return render_template('opcion_no.html')

#Recuperamos los datos de crecimiento y tratamento para resolver el modelo de ecuaciones diferenciales que modelan el crecimiento tumoral 
@app.route('/procesar_si', methods=['POST'])
def procesar_si():
    tipo = int(request.form['tipo'])
    v1 = float(request.form['v1'])
    v2 = float(request.form['v2'])
    paso1 = float(request.form['paso1'])
    rad = float(request.form['rad'])
    dir1 = float(request.form['dir1'])
    dtr = int(request.form['dtr'])
    der = int(request.form['der'])
    ser = int(request.form['ser'])
    med = float(request.form['med'])
    dim = int(request.form['dim'])
    dtm = int(request.form['dtm'])
    dem = int(request.form['dem'])
    sem = int(request.form['sem'])
    s_2 = float(request.form['cel1'])
    s_1 = float(request.form['cel2'])
    dic = int(request.form['dic'])
    dec = int(request.form['dec'])
    sec = int(request.form['sec'])
    
    def Crecimientosi(V0,V1,T1):
        T0 = 0.0
        b = 0.072
        V_max = 15

        CreEx = np.log(V1/V0)/(T1-T0)
        Creci = (V1-V0)/(T1-T0)
        CreciLo = (1/T1)*np.log( 1/(((V_max*V0)/V1 -V0)*(1/(V_max-V0))))
        CVon = (V1**(1/3) - (V1**(1/3))*np.exp(-b*T1/3))/(1/b - np.exp(-b*T1/3)/b )
        Gom = (b*np.log(V1/V0))/(1- np.exp(-b*T1))

        return V0, CreEx,Creci,CreciLo, CVon, Gom
        
    V0, CreEx,Creci,CreciLo, CVon, Gom = Crecimientosi(v1, v2,paso1)


    def Tipod(U):
        if tipo == 1 or tipo ==2:
            dU = CreEx*U
            return dU
            #--------------------------------------------------------------------------------------------------------------
        elif tipo == 3 or tipo ==4:
            V_max = 15
            dU = CreciLo*U*(1 - (U/V_max))
            return dU
            #--------------------------------------------------------------------------------------------------------------
        elif tipo == 5 or tipo ==6 or tipo ==7:
            b = 0.072
            dU = CVon*(U**(2/3)) - b*U
            return dU
            #--------------------------------------------------------------------------------------------------------------
        elif tipo == 8 or tipo ==9:
            dU = Creci
            return dU
            #--------------------------------------------------------------------------------------------------------------
        elif tipo == 10 or tipo ==11 or tipo ==12:
            b = 10
            dU =  Gom*U*np.log(b/U)
        return dU
    
    def RK4(ini,f):
        TEx = {'tinit':0, 'tfinal':365, 'h':0.5}
        t = np.arange(TEx['tinit'],TEx['tfinal'],TEx['h'])
        U = np.zeros(len(t))
        U[0] = ini
        for i in range(1,len(t)):
            k1 = f(U[i-1])*TEx['h']
            k2 = f(U[i-1]+k1/2)*TEx['h']
            k3 = f(U[i-1]+k2/2)*TEx['h']
            k4 = f(U[i-1]+k3)*TEx['h']
            U[i] = U[i-1]+(k1+2*k2+2*k3+k4)/6
        return t, U.transpose()

    t, solucion = RK4(float(V0),Tipod)
    
    #RK4 para los tratamientos
    def RK4tratamientos(ic,f1,f2, f3, tpr1, descanso1,ddt1, sesiones): 
        tpr = 2*tpr1
        descanso = 2*descanso1
        ddt = 2*ddt1
        p = {'tinit':0, 'tfinal':365, 'h':0.5}
        t = np.arange(p['tinit'],p['tfinal'],p['h'])
        U = np.zeros((len(t),len(ic)))
        U[0] = ic
        for i in range(1,len(t)):
            for j in range(0, sesiones):
                if i < tpr:
                    #f1 es la funcion de crecimiento normal del tumor sin tratamiento 
                    k1 = f1(U[i-1])*p['h']
                    k2 = f1(U[i-1]+k1/2)*p['h']
                    k3 = f1(U[i-1]+k2/2)*p['h']
                    k4 = f1(U[i-1]+k3)*p['h']
                    U[i] = U[i-1]+(k1+2*k2+2*k3+k4)/6
                elif  tpr + j*(ddt +descanso) <= i < tpr + ddt*(j+1) + descanso*j:
                    #f2 es la funcion con la primera inyeccion de las celulas dentriticas 
                    #f2 es la funcion con la primera inyeccion de las celulas dentriticas 
                    k1 = f2(U[i-1])*p['h']
                    k2 = f2(U[i-1]+k1/2)*p['h']
                    k3 = f2(U[i-1]+k2/2)*p['h']
                    k4 = f2(U[i-1]+k3)*p['h']
                    U[i] = U[i-1]+(k1+2*k2+2*k3+k4)/6 
                elif  tpr + ddt*(j+1) + descanso*j <= i < tpr + ddt*(j+1) + descanso*(j+1):
                    #f3 es la func
                    k1 = f3(U[i-1])*p['h']
                    k2 = f3(U[i-1]+k1/2)*p['h']
                    k3 = f3(U[i-1]+k2/2)*p['h']
                    k4 = f3(U[i-1]+k3)*p['h']
                    U[i] = U[i-1]+(k1+2*k2+2*k3+k4)/6
        return t, U
    
    #Radioterapia--------------------------------------------------
    gamma = 0.24 #dia-1
    d_r = 0.0432
    Dosis = rad
    
    def Radioterapiasi(M):
        U,e_r,R = M
        dU = Tipod(U) - e_r*R*U
        de_r = d_r*e_r*(1-e_r)
        dR = -gamma*R
        return np.array([dU, de_r, dR])

    def efirad(Dosis):
        e_0 = 0.0912 /(Dosis**(1/2))
        return e_0

    #funcion considerando la dosis de medicamento
    def Radioterapiact(M):
        U,e_r,R = M
        dU = Tipod(U) - e_r*R*U
        de_r = E_01 + d_r*e_r*(1-e_r)
        dR = Dosis -gamma*R
        return np.array([dU, de_r, dR])
    
    E_01 = efirad(rad)
    
    t, radioterapia = RK4tratamientos([V0,0,0],Tipod,Radioterapiact,Radioterapiasi, dir1, der, dtr, ser)
    
    #Quimioterapia -------------------------------------------------------
    beta = 0.24 #dia-1
    d_c = 0.3504
    concentracion = med
    
    def Quimioterapiasi(M):
        U,e_c,Q = M
        dU = Tipod(U) - e_c*Q*U
        de_c = d_c*e_c*(1-e_c)
        dQ = -beta*Q
        return np.array([dU, de_c, dQ])

    def efiini(Concentracion):
        e_0 = 0.1824 /(Concentracion**0.95)
        return e_0

    #funcion considerando la dosis de medicamento
    def Quimioterapiact(M):
        U,e_c,Q = M
        dU = Tipod(U) - e_c*Q*U
        de_c = e_0 + d_c*e_c*(1-e_c)
        dQ = concentracion -beta*Q
        return np.array([dU, de_c, dQ])

    e_0 = efiini(concentracion)
    
    t, quimioterapia = RK4tratamientos([V0,0,0],Tipod,Quimioterapiact,Quimioterapiasi, dim, dem, dtm, sem)
    
    #Inmunoterapia--------------------------------------
    g_1 = (2*10**7)/1000
    g_2 = (1*10**5)/1000
    g_3 = (1*10**3)/1000
    mu_2 = 0.03
    p_1=0.1245
    r_2=0.18
    b=1*10**(-9)
    a=1
    mu_3=10
    c= 8.55*10**(-5)
    p_2 = 5
    s_1 = 0.01
    s_2 = 0.001
    def Inmunoterapiaci(M):
        U,E,I = M
        dU = Tipod(U) -  (a*E*U)/(g_2 + U)
        dE = c*U - mu_2*E + (p_1*E*I)/(g_1 + I) + s_1 
        dI = - mu_3*I + (p_2*E*U)/(g_3 + U) + s_2 
        return np.array([dU, dE, dI])

    def Inmunoterapiasi(M):
        U,E,I = M
        dU = Tipod(U) -  (a*E*U)/(g_2 + U)
        dE = c*U - mu_2*E + (p_1*E*I)/(g_1 + I)  
        dI = - mu_3*I + (p_2*E*U)/(g_3 + U) 
        return np.array([dU, dE, dI])

    t, inmunoterapia= RK4tratamientos([V0,0,0],Tipod,Inmunoterapiaci,Inmunoterapiasi, dic, dec, 1, sec)
    

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes[0, 0].set_title("Sin tratamiento")
    axes[0, 0].set_xlabel("Tiempo (dias)")
    axes[0, 0].set_ylabel("Volumen ($cm^3$)")
    axes[0, 0].plot(t, solucion, 'palevioletred')


    axes[0, 1].set_title("Radioterapia")
    axes[0, 1].set_xlabel("Tiempo (dias)")
    axes[0, 1].set_ylabel("Volumen ($cm^3$)")
    axes[0, 1].plot(t, radioterapia[:,0], 'teal')


    axes[1, 0].set_title("Quimioterapia")
    axes[1, 0].set_xlabel("Tiempo (dias)")
    axes[1, 0].set_ylabel("Volumen ($cm^3$)")
    axes[1, 0].plot(t, quimioterapia[:,0], 'yellowgreen')


    axes[1, 1].set_title("Inmunoterapia")
    axes[1, 1].set_xlabel("Tiempo (dias)")
    axes[1, 1].set_ylabel("Volumen ($cm^3$)")
    axes[1, 1].plot(t, inmunoterapia[:,0], 'skyblue')
    
    fig.tight_layout()
    plt.savefig('static/my_plot.png')
    
    return render_template('result.html', plot_url ="static/my_plot.png")

@app.route('/procesar_no', methods=['POST'])
def procesar_no():
    tipo = int(request.form['tipo'])
    l1 = float(request.form['l1'])
    l2 = float(request.form['l2'])
    a1 = float(request.form['a1'])
    a2 = float(request.form['a2'])
    paso2 = float(request.form['paso2'])
    rad = float(request.form['rad'])
    dir1 = float(request.form['dir1'])
    dtr = int(request.form['dtr'])
    der = int(request.form['der'])
    ser = int(request.form['ser'])
    med = float(request.form['med'])
    dim = int(request.form['dim'])
    dtm = int(request.form['dtm'])
    dem = int(request.form['dem'])
    sem = int(request.form['sem'])
    s_2 = float(request.form['cel1'])
    s_1 = float(request.form['cel2'])
    dic = int(request.form['dic'])
    dec = int(request.form['dec'])
    sec = int(request.form['sec'])

    # Aquí puedes realizar cualquier procesamiento o cálculos necesarios con los datos
    def Crecimientono(L0, L1, A0,A1,T1):
        V0 = (L0*(A0**2))/2
        V1 = (L1*(A1**2))/2
        T0 = 0.0
        b = 0.072
        V_max = 15

        CreEx = np.log(V1/V0)/(T1-T0)
        Creci = (V1-V0)/(T1-T0)
        CreciLo = (1/T1)*np.log( 1/(((V_max*V0)/V1 -V0)*(1/(V_max-V0))))
        CVon = (V1**(1/3) - (V1**(1/3))*np.exp(-b*T1/3))/(1/b - np.exp(-b*T1/3)/b )
        Gom = (b*np.log(V1/V0))/(1- np.exp(-b*T1))

        return V0, CreEx,Creci,CreciLo, CVon, Gom
    
    V0, CreEx,Creci,CreciLo, CVon, Gom = Crecimientono(l1, l2, a1,a2,paso2)
    
    def Tipod(U):
        if tipo == 1 or tipo ==2:
            dU = CreEx*U
            return dU
            #--------------------------------------------------------------------------------------------------------------
        elif tipo == 3 or tipo ==4:
            V_max = 15
            dU = CreciLo*U*(1 - (U/V_max))
            return dU
            #--------------------------------------------------------------------------------------------------------------
        elif tipo == 5 or tipo ==6 or tipo ==7:
            b = 0.072
            dU = CVon*(U**(2/3)) - b*U
            return dU
            #--------------------------------------------------------------------------------------------------------------
        elif tipo == 8 or tipo ==9:
            dU = Creci
            return dU
            #--------------------------------------------------------------------------------------------------------------
        elif tipo == 10 or tipo ==11 or tipo ==12:
            b = 10
            dU =  Gom*U*np.log(b/U)
        return dU
    
    def RK4(ini,f):
        TEx = {'tinit':0, 'tfinal':365, 'h':0.5}
        t = np.arange(TEx['tinit'],TEx['tfinal'],TEx['h'])
        U = np.zeros(len(t))
        U[0] = ini
        for i in range(1,len(t)):
            k1 = f(U[i-1])*TEx['h']
            k2 = f(U[i-1]+k1/2)*TEx['h']
            k3 = f(U[i-1]+k2/2)*TEx['h']
            k4 = f(U[i-1]+k3)*TEx['h']
            U[i] = U[i-1]+(k1+2*k2+2*k3+k4)/6
        return t, U.transpose()
    
    
    t, solucion = RK4(float(V0),Tipod)
    
    #RK4 para los tratamientos
    def RK4tratamientos(ic,f1,f2, f3, tpr1, descanso1,ddt1, sesiones): 
        tpr = 2*tpr1
        descanso = 2*descanso1
        ddt = 2*ddt1
        p = {'tinit':0, 'tfinal':365, 'h':0.5}
        t = np.arange(p['tinit'],p['tfinal'],p['h'])
        U = np.zeros((len(t),len(ic)))
        U[0] = ic
        for i in range(1,len(t)):
            for j in range(0, sesiones):
                if i < tpr:
                    #f1 es la funcion de crecimiento normal del tumor sin tratamiento 
                    k1 = f1(U[i-1])*p['h']
                    k2 = f1(U[i-1]+k1/2)*p['h']
                    k3 = f1(U[i-1]+k2/2)*p['h']
                    k4 = f1(U[i-1]+k3)*p['h']
                    U[i] = U[i-1]+(k1+2*k2+2*k3+k4)/6
                elif  tpr + j*(ddt +descanso) <= i < tpr + ddt*(j+1) + descanso*j:
                    #f2 es la funcion con la primera inyeccion de las celulas dentriticas 
                    #f2 es la funcion con la primera inyeccion de las celulas dentriticas 
                    k1 = f2(U[i-1])*p['h']
                    k2 = f2(U[i-1]+k1/2)*p['h']
                    k3 = f2(U[i-1]+k2/2)*p['h']
                    k4 = f2(U[i-1]+k3)*p['h']
                    U[i] = U[i-1]+(k1+2*k2+2*k3+k4)/6 
                elif  tpr + ddt*(j+1) + descanso*j <= i < tpr + ddt*(j+1) + descanso*(j+1):
                    #f3 es la func
                    k1 = f3(U[i-1])*p['h']
                    k2 = f3(U[i-1]+k1/2)*p['h']
                    k3 = f3(U[i-1]+k2/2)*p['h']
                    k4 = f3(U[i-1]+k3)*p['h']
                    U[i] = U[i-1]+(k1+2*k2+2*k3+k4)/6
        return t, U
    
    #Radioterapia--------------------------------------------------
    gamma = 0.24 #dia-1
    d_r = 0.0432
    Dosis = rad
    
    def Radioterapiasi(M):
        U,e_r,R = M
        dU = Tipod(U) - e_r*R*U
        de_r = d_r*e_r*(1-e_r)
        dR = -gamma*R
        return np.array([dU, de_r, dR])

    def efirad(Dosis):
        e_0 = 0.0912 /(Dosis**(1/2))
        return e_0

    #funcion considerando la dosis de medicamento
    def Radioterapiact(M):
        U,e_r,R = M
        dU = Tipod(U) - e_r*R*U
        de_r = E_01 + d_r*e_r*(1-e_r)
        dR = Dosis -gamma*R
        return np.array([dU, de_r, dR])
    
    E_01 = efirad(rad)
    
    t, radioterapia = RK4tratamientos([V0,0,0],Tipod,Radioterapiact,Radioterapiasi, dir1, der, dtr, ser)
    
    #Quimioterapia -------------------------------------------------------
    beta = 0.24 #dia-1
    d_c = 0.3504
    concentracion = med
    
    def Quimioterapiasi(M):
        U,e_c,Q = M
        dU = Tipod(U) - e_c*Q*U
        de_c = d_c*e_c*(1-e_c)
        dQ = -beta*Q
        return np.array([dU, de_c, dQ])

    def efiini(Concentracion):
        e_0 = 0.1824 /(Concentracion**0.95)
        return e_0

    #funcion considerando la dosis de medicamento
    def Quimioterapiact(M):
        U,e_c,Q = M
        dU = Tipod(U) - e_c*Q*U
        de_c = e_0 + d_c*e_c*(1-e_c)
        dQ = concentracion -beta*Q
        return np.array([dU, de_c, dQ])

    e_0 = efiini(concentracion)
    
    t, quimioterapia = RK4tratamientos([V0,0,0],Tipod,Quimioterapiact,Quimioterapiasi, dim, dem, dtm, sem)
    
    #Inmunoterapia--------------------------------------
    g_1 = (2*10**7)/1000
    g_2 = (1*10**5)/1000
    g_3 = (1*10**3)/1000
    mu_2 = 0.03
    p_1=0.1245
    r_2=0.18
    b=1*10**(-9)
    a=1
    mu_3=10
    c= 8.55*10**(-5)
    p_2 = 5
    s_1 = 0.01
    s_2 = 0.001
    def Inmunoterapiaci(M):
        U,E,I = M
        dU = Tipod(U) -  (a*E*U)/(g_2 + U)
        dE = c*U - mu_2*E + (p_1*E*I)/(g_1 + I) + s_1 
        dI = - mu_3*I + (p_2*E*U)/(g_3 + U) + s_2 
        return np.array([dU, dE, dI])

    def Inmunoterapiasi(M):
        U,E,I = M
        dU = Tipod(U) -  (a*E*U)/(g_2 + U)
        dE = c*U - mu_2*E + (p_1*E*I)/(g_1 + I)  
        dI = - mu_3*I + (p_2*E*U)/(g_3 + U) 
        return np.array([dU, dE, dI])

    t, inmunoterapia= RK4tratamientos([V0,0,0],Tipod,Inmunoterapiaci,Inmunoterapiasi, dic, dec, 1, sec)
    

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes[0, 0].set_title("Sin tratamiento")
    axes[0, 0].set_xlabel("Tiempo (dias)")
    axes[0, 0].set_ylabel("Volumen ($cm^3$)")
    axes[0, 0].plot(t, solucion, 'palevioletred')


    axes[0, 1].set_title("Radioterapia")
    axes[0, 1].set_xlabel("Tiempo (dias)")
    axes[0, 1].set_ylabel("Volumen ($cm^3$)")
    axes[0, 1].plot(t, radioterapia[:,0], 'teal')


    axes[1, 0].set_title("Quimioterapia")
    axes[1, 0].set_xlabel("Tiempo (dias)")
    axes[1, 0].set_ylabel("Volumen ($cm^3$)")
    axes[1, 0].plot(t, quimioterapia[:,0], 'yellowgreen')


    axes[1, 1].set_title("Inmunoterapia")
    axes[1, 1].set_xlabel("Tiempo (dias)")
    axes[1, 1].set_ylabel("Volumen ($cm^3$)")
    axes[1, 1].plot(t, inmunoterapia[:,0], 'skyblue')
    
    fig.tight_layout()
    plt.savefig('static/my_plot.png')
    
    return render_template('result.html', plot_url ="static/my_plot.png")


if __name__ == '__main__':
    app.run('127.0.0.1', 5000, debug = False)


# In[ ]:




