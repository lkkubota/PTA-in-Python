
import numpy as np
import pandas as pd

def ReturnBuildUp(t_ini,dt,tt,pp):

#Extracts build-up segment from time and pressure array
#Code written by L.Kubota
        
#Parameters
#----------
#tt: time array (k,)
#pp: bhp array (k,) 
#t_ini: time instant in "tt" where build-up starts (float)
#dt: build-up duration (float)
        
#Returns
#-------
#t: build-up time array (m,)
#p: build-up pressure array (m,)   

    t=np.linspace(t_ini,t_ini+dt,500)
    p=np.interp(t,tt,pp)
    return t,p


def calc_derivative(time, pressure, factor_L):

#Bourdet logarithmic pressure derivative
#Code written by R.M. Araujo
        
#Parameters
#----------
#time: build-up time array (k,)
#pressure: build-up bhp array (k,) 
#factor_L: smoothing factor (float)
        
#Returns
#-------
#deriv_pressure: bourdet derivative (m,)
           
    n_points = len(time)
    deriv_pressure = []
    deriv_time = []
    i = 0
    while i < n_points:
        t1 = time[i]
        p1 = pressure[i]
        # find ti-1
        j = i
        while j > 0:
            if time[j] < t1 / np.exp(factor_L):
                break
            j -= 1
          # find ti+1
        k = i
        while k < n_points-1:
            if time[k] > t1 * np.exp(factor_L):
                break
            k += 1

        p0, p2 = pressure[j], pressure[k]
        t0, t2 = time[j], time[k]
        log_t0 = np.log(t0) if t0 != 0 else 0
        log_t1 = np.log(t1) if t1 != 0 else 0
        log_t2 = np.log(t2) if t2 != 0 else 0
        w1 = log_t1 - log_t0
        w2 = log_t2 - log_t1
        m1 = (p1-p0)/w1 if w1 > 0 else 0
        m2 = (p2-p1)/w2 if w2 > 0 else 0
        tdpdt = m1*w2/(w1+w2) + m2*w1/(w1+w2)
        deriv_pressure.append(tdpdt)
        i+=1
    return np.array(deriv_pressure)

def skin_factor_pdg(tr,dt,dp,der,dp_radial_ref,der_radial_ref):
#Code written by L.Kubota
        
#Parameters
#----------
#tr: any time in IARF plateau (float)
#dt: build-up elapsed time (t-t[0]) array from log-log plot(k,)
#dp: delta pressure (p-p[0]) array from log-log plot (k,)
#der: bourdet pressure derivative 'dpdlnt' array from log-log plot (k,) 
#dp_radial_ref: 'dp' extracted at 't=tr' from the log-log plot belonging
#    to the reference build-up (float)
#der_radial_ref: 'dpdlnt' extracted at 't=tr' from the log-log plot belonging
#    to the reference build-up (float)
         
#Returns
#-------
#s: skin factor (float)
           
   
    dp_radial=np.interp(tr,dt,dp)
    der_radial=np.interp(tr,dt,der)
    dp_skin=dp_radial-dp_radial_ref
    s=0.5*dp_skin/der_radial_ref
    return s

def teq_agarwal(t,q,delta_t):

#Agarwal equivalent time for multiple rates 
#Code written by L.Kubota
#
#Parameters
#----------
#t: time array before the start of build-up (k,)
#q: flow-rate array before the start of build-up (k,)
#delta_t: elapsed-time array during build-up (t-t[0]) (m,)

#Returns
#-------
#teq: Agarwal equivalent time(m,)
           
      
    t=np.array(t)
    q=np.array(q)
    t=np.append(np.zeros(1),t)
    q=np.append(np.zeros(1),q)

    dt=t[-1]-t[:-1]
    dq=np.diff(q)
    dq=dq/q[-1]

    teq=[]
    i=0
    
    for time in delta_t:
        prod0=1
        for i in range(len(dq)):
            prod=((dt[i]/(time+dt[i]))**dq[i])*prod0
            prod0=prod
        teq.append(prod*time)

    return np.array(teq)

def build_up_cut_out(t,p,t_start,t_end):
      
#Extracts build-up segments from t and p 
#Code written by L.Kubota
        
#Parameters
#----------
#t: time array of entire well-test (k,)
#p: pressure array of entire well-test (k,)
#t_start: time associated with the start of build-up (float)
#t_end: time associated with the end of build-up (float)
        
#Returns
#-------
#t_bu: time array for selected build-up (m,)
#pw_bu: pressure array for selected build-up (m,)
      
    t_bu=t[(np.where((t >= t_start) & (p <= t_end)))]
    pw_bu=p[(np.where((t >= t_start) & (p <= t_end)))]    
    return t_bu,pw_bu


def k_and_s(t_radial,delta_tbu,dp_bu,der,B,visc,h,por,ct,rw):
#Computes permeability and skin for build-ups
#Code written by L.Kubota    
    der_radial=np.interp(t_radial,delta_tbu,der)
    dp_radial=np.interp(t_radial,delta_tbu,dp_bu)
    k=0.5*19.03*(B*visc/h)*(1/der_radial)
    dp_radial_ref=1.151*2*der_radial*(np.log10(t_radial)+np.log10((0.0003484*k)/(por*visc*ct*rw*rw))+0.3514)
    der_radial_ref=der_radial
    dpq_radial_ref=dp_radial
    dp_skin=dp_radial-dp_radial_ref
    s=0.5*dp_skin/der_radial_ref
    return k,s 

def main():
	
	
 if __name__ == "__main__":
	 main()