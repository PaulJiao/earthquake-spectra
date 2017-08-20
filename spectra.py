import numpy as np
import matplotlib.pyplot as plt
class spectra:
    accpath=''
    accwave=np.zeros((1,1))
    damp=np.array([[0,0.05,0.1]])
    dt=0.01
    def __init__(self,accpath,dt=0.01):
        self.accpath=accpath
        self.accwave=np.loadtxt(self.accpath)
        self.accwave=self.accwave.reshape((self.accwave.size,1))
        self.dt=dt
    def figure(self):
        data=self.spectradata()
        period=data[:,0]
        data1=data[:,1:3]
        data2=data[:,5:8]
        data3=data[:,10:]
        plt.subplot(1,3,1)
        plt.plot(period,data1[:,1],'-')
        plt.subplot(1,3,2)
        plt.plot(period,data2[:,1],'-')
        plt.subplot(1,3,3)
        plt.plot(period,data3[:,1],'-')
        plt.show()

    

    def spectradata(self):
        matresult=np.zeros((50,1))
        tmpt1=matresult
        for i in self.damp.flat:
            tmpt=self.__getdata__(i)
            matresult=np.append(matresult,tmpt,1)
        return np.delete(matresult,tmpt1,1)

    def __getdata__(self,damp):
        period=np.linspace(0.0001,10,50).reshape(50,1)
        omega=1./period*2*np.pi
        results=np.empty((1,3))
        temy=results
        for oneomega in omega.flat:
            results=np.append(results,self.__getMAx__(oneomega,damp),axis=0)
        
        results=np.delete(results,temy,axis=0)       
        results=np.append(period,results,1)       
        return results

    def __getMAx__(self,omega,damp):
        dt=self.dt
        s=np.sin(omega*np.sqrt(1-damp**2)*dt)
        c=np.cos(omega*np.sqrt(1-damp**2)*dt)
        et=np.exp(-damp*omega*dt)
        damp1=damp/np.sqrt(1-damp**2)
        dampo2=(2*damp**2-1)/omega**2/dt
        dampo3=2*damp/omega**3/dt
        a11=et*(damp1*s+c)
        a12=et*s/omega/np.sqrt(1-damp**2)
        a21=-omega/np.sqrt(1-damp**2)*et*s
        a22=et*(c-damp1*s)
        b11=et*((dampo2+damp/omega)*s/omega/np.sqrt(1-damp**2)+(dampo3+1/omega**2)*c)-dampo3
        b12=-et*(dampo2*s/omega/np.sqrt(1-damp**2)+dampo3*c)-1/omega**2+dampo3
        b21=et*((dampo2+damp/omega)*(c-damp1*s)-(dampo3+1/omega**2)*(omega*np.sqrt(1-damp**2)*s+damp*omega*c))+1/omega**2/dt
        b22=-et*(dampo2*(c-damp1*s)-dampo3*(omega*np.sqrt(1-damp**2)*s+omega*damp*c))-1/omega**2/dt
        A=np.array([[a11,a12],[a21,a22]])
        B=np.array([[b11,b12],[b21,b22]])
        vol=np.zeros((self.accwave.size,1))
        disp=np.zeros((self.accwave.size,1))
        acc2=np.zeros((self.accwave.size,1))
        for i in range(self.accwave.size-1):
            acc2[i,0]=-(2*damp*omega*vol[i,0]+(omega**2)*disp[i,0])
            tmpt=self.__onestep__(self.accwave[i,0],self.accwave[i+1,0],vol[i,0],disp[i,0],A,B)
            disp[i+1,0]=tmpt[0,0]
            vol[i+1,0]=tmpt[1,0]
            acc2[i+1,0]=-(2*damp*omega*vol[i+1,0]+omega**2*disp[i+1,0])
        print(tmpt)
        SD=np.amax(np.absolute(disp))
        SV=np.amax(np.absolute(vol))
        SA=np.amax(np.absolute(acc2))
        result=np.array([[SD,SV,SA]])
        return result
    def __onestep__(self,acc1,acc2,vol1,disp1,A,B):
        x1=np.array([[disp1],[vol1]])
        a1=np.array([[acc1],[acc2]])
        x2=np.dot(A,x1)+np.dot(B,a1)
        return x2
        
a=spectra('EQ.txt',0.01) 
a.figure()
#print(a.spectradata())
#print(a.spectradata())
#k=a.spectradata()       
#print(k.shape()) 
        
        
        

            







    



