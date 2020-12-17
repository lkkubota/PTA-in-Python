import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def ReturnBuildUp(t_ini,dt,tt,pp):
    t=np.linspace(t_ini,t_ini+dt,500)
    p=np.interp(t,tt,pp)
    return t,p

def read_file(filename):
  from google.colab import files
  uploaded = files.upload()
  import io
  df = pd.read_csv(io.BytesIO(uploaded[filename]), header = 0, sep = ';')
  t=df.time
  p=df.pressure
  q=df.qo
  w=df.bsw
  return t,q,p,w


def calc_derivative(time, pressure, factor_L):
    # inicializacao do fator_L e variaveis secundarias
    
    n_points = len(time)
    deriv_pressure = []
    deriv_time = []
    i = 0
    while i < n_points:
        t1 = time[i]
        p1 = pressure[i]
        # encontrar o ti-1
        j = i
        while j > 0:
            if time[j] < t1 / np.exp(factor_L):
                break
            j -= 1
        # encontrar o ti+1
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
        # calculo da derivada
        deriv_pressure.append(tdpdt)
        i+=1
        # retorna ambos os arrays (t e dpdt) como arrays numpy
    return np.array(deriv_pressure)

def plot_analysis(t,p,dt,dp,der,s,shut,date,xmin,xmax,ymin,ymax):
    f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,4))
    #f, (ax1,ax2) = plt.subplots(1,2)
    
    f.subplots_adjust(wspace=0.7)
    ax1.scatter(t,p,color='blue', s=100,edgecolors=(0, 0, 0))
    ax1.set_title(f"Shut-in {date}: {shut} hours")
    ax1.set_xlabel('t [h]')
    ax1.set_ylabel('p [kgf/cm²]')
    ax2.scatter(dt,dp,color='red', s=100,edgecolors=(0, 0, 0))
    ax2.scatter(dt,der,color='blue', s=100,edgecolors=(0, 0, 0))
    ax2.set_title(f"Skin: {round(s,1)}")
    ax2.set_xlabel('t [h]')
    ax2.set_ylabel('dp/q')
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(xmin,xmax)
    ax2.set_ylim(ymin,ymax)

    plt.tight_layout()
    plt.grid()
    plt.show()  
    
def ShowAllLogLog(df_t,df_dp,df_der,xmin,xmax,ymin,ymax,well):
    plt.plot(figsize=(6,6))
    plt.loglog(df_t.iloc[:,:-1],df_der.iloc[:,:-1],'ko',df_t.iloc[:,:-1],df_dp.iloc[:,:-1],'ko',alpha=0.04)
    plt.loglog(df_t.iloc[:,-1],df_der.iloc[:,-1],'ro',alpha=0.5,label='last buil-up')    
    plt.loglog(df_t.iloc[:,-1],df_dp.iloc[:,-1],'ro',alpha=0.5)    
    
    plt.legend()
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.xlabel('time [h]')
    plt.ylabel('dp/q [kgf/cm²/m³/d]')
    plt.grid()
    plt.title(f"{well}: Diagnostic Plots")
    plt.savefig(f"{well} Diagnostic Plots",dpi=200)
    plt.tight_layout()
    plt.show()


def PlotBHPandRateandPavg(tt,pp,qq,pavg,t2,qmin,qmax):
    fig, ax1 = plt.subplots(figsize=(10,5))

    ax1.scatter(tt, pp, c='b',s=10)
    ax1.scatter(t2,pavg,c='r',s=40,label='pavg [kgf/cm²]')
    ax1.set_xlabel('time (hrs)')

    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('pressure [kgf/cm²]', color='k')
    #ax1.set_ylim(400.0,600.0)
    ax1.tick_params('y', colors='b')

    ax2 = ax1.twinx()
    ax3 = ax1.figure.add_axes(ax1.get_position(True), sharex=ax1, sharey=ax1,
            frameon=False)
    ax3 = ax2.figure.add_axes(ax2.get_position(True), sharex=ax2, sharey=ax2,
            frameon=False)

    ax3.xaxis.set_visible(False)
    ax3.yaxis.set_visible(False)

    #ax2.plot(tt, qq, 'k-')
    ax2.scatter(tt, qq, c='g',s=8)
    ax2.set_ylabel('rate [m³/d]', color='k')
    ax2.set_ylim(qmin,qmax)
    ax2.tick_params('y', colors='k')
    ax1.legend(loc='lower right')
    

    for i in range(len(pavg)):
        ax1.annotate(str(round(pavg[i],1)), xy=(t2[i],pavg[i]),fontsize=14)
    plt.savefig("pavg",dpi=300)
    plt.tight_layout()
    plt.show()
    
def PlotBHPandRateandSkin(tt,pp,skin,t2,smin,smax):
    fig, ax1 = plt.subplots(figsize=(10,5))

    ax1.scatter(tt, pp, c='b',s=10)
    ax1.set_xlabel('time (hrs)')

    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('pressure [kgf/cm²]', color='k')
    ax1.tick_params('y', colors='b')

    ax2 = ax1.twinx()
    ax3 = ax1.figure.add_axes(ax1.get_position(True), sharex=ax1, sharey=ax1,
            frameon=False)
    ax3 = ax2.figure.add_axes(ax2.get_position(True), sharex=ax2, sharey=ax2,
            frameon=False)

    ax3.xaxis.set_visible(False)
    ax3.yaxis.set_visible(False)

    
    ax2.scatter(t2, skin, c='g',s=40,label='skin')
    ax2.set_ylabel('skin', color='k')
    ax2.set_ylim(smin,smax)
    ax2.tick_params('y', colors='k')
    ax2.legend(loc='lower right')

    for i in range(len(skin)):
        ax2.annotate(str(round(skin[i],1)), xy=(t2[i],skin[i]),fontsize=14)
    plt.savefig("skin",dpi=300)
    plt.tight_layout()
    plt.show()    

def PlotBHPandRateandIP(tt,pp,IP,t2,IPmin,IPmax):
    fig, ax1 = plt.subplots(figsize=(10,5))

    ax1.scatter(tt, pp, c='b',s=10)
    ax1.set_xlabel('time (hrs)')

    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('pressure [kgf/cm²]', color='k')
    ax1.tick_params('y', colors='b')

    ax2 = ax1.twinx()
    ax3 = ax1.figure.add_axes(ax1.get_position(True), sharex=ax1, sharey=ax1,
            frameon=False)
    ax3 = ax2.figure.add_axes(ax2.get_position(True), sharex=ax2, sharey=ax2,
            frameon=False)

    ax3.xaxis.set_visible(False)
    ax3.yaxis.set_visible(False)

    
    ax2.scatter(t2, IP, c='g',s=40,label='IP ')
    ax2.set_ylabel('IP [m³/d/kgf/cm²]', color='k')
    ax2.set_ylim(IPmin,IPmax)
    ax2.tick_params('y', colors='k')
    ax2.legend(loc='lower right')

    for i in range(len(IP)):
        ax2.annotate(str(round(IP[i],1)), xy=(t2[i],IP[i]),fontsize=14)
    plt.savefig("ip",dpi=300)
    plt.tight_layout()
    plt.show()

def SkinBSW(skin,water,well):
    water=np.array(water)
    water=water*100
    plt.plot(figsize=(6,6))
    plt.plot(skin,water,'ro',markersize=10,markeredgecolor=(0,0,0))
    plt.title(f"{well} Skin vs BSW")
    plt.ylabel("BSW [%]",fontsize=14)
    plt.xlabel("skin",fontsize=14)
    plt.grid()
    plt.savefig(f"Skin_vs_BSW",dpi=300)
    plt.tight_layout()
    plt.show()

def SkinIP(skin,IP,well,ymin,ymax):
    
    plt.plot(figsize=(6,6))
    plt.plot(skin,IP,'ro',markersize=10,markeredgecolor=(0,0,0))
    plt.title(f"{well} Skin vs IP")
    plt.ylabel("IP [m³/d/kgf/cm²]",fontsize=14)
    plt.ylim(ymin,ymax)
    plt.xlabel("skin",fontsize=14)
    plt.grid()
    plt.savefig(f"Skin_vs_IP",dpi=300)
    plt.tight_layout()
    plt.show()





def main():
	
	
 if __name__ == "__main__":
	 main()