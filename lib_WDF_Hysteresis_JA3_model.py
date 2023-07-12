import numpy as np
import warnings
warnings.simplefilter('error')
from scipy.optimize import newton


#
# lib_WDF_Hysteresis_JA_model
#

#WDF_Hysteresis_JA_model_version = "3.0.0.2021.10.22"
WDF_Hysteresis_JA_model_version = "3.1.0.2023.7.10 from 3.0.0.2021.10.22"
USE_JA_VER_3_1 = True

def get_WDF_Hystereys_JA_model_version():
    return WDF_Hysteresis_JA_model_version


class Hysteresis_JA_model:

    def __init__(self, Ms, hys_a, hys_kp, alpha, cr, lambda1, Cal_Gain1, Cal_Gain2, Hs, Ns, Ld, Sd):
 
        self.Ms       = Ms
        self.hys_a    = hys_a
        self.kp       = hys_kp
        self.kp_b     = 1/hys_kp
        self.alpha    = alpha
        self.cr       = cr
        self.lambda1  = lambda1

        self.u0       = 1.256637062e-6

        self.Cal_Gain1 = Cal_Gain1
        self.Cal_Gain2 = Cal_Gain2

        self.Hs = Hs
        self.Ns = Ns
        self.Ld = Ld
        self.Sd = Sd

        self.H   = 0
        self.M   = 0
        self.B   = 0
        self.Phi = 0

        self.dB   = 0
        self.dPhi = 0
       
        self.Old_H   = 0
        self.Old_M   = 0
        self.Old_B   = 0
        self.Old_Phi = 0


        self.Hd   = 0
        self.He   = 0
        self.Man  = 0

        self.Old_He = 0
        self.Old_Man = 0

        if USE_JA_VER_3_1:
            self.limit_HeAbs = 1e-4 * self.hys_a    # if He_Abs < self.limit_HeAbs: Man = 0 (Man=3.33e-5)


    def getHysteresis(self, Current,Fs):

        NI =  Current*self.Ns*self.Cal_Gain1 # scaling

        self.H =  NI/self.Ld
#        print("H=", self.H)

        self.Hd =self.H - self.Old_H

        self.He = self.H+self.alpha*self.Old_M # https://doc.comsol.com/5.5/doc/com.comsol.help.acdc/acdc_ug_theory.05.14.html
        self.He = self.He + 3.0*self.lambda1/(self.u0*self.Ms**2)*self.Old_M
        dHe = self.He - self.Old_He


        He_abs = np.abs(self.He)

        if USE_JA_VER_3_1:
            if He_abs < self.limit_HeAbs:   # この条件の時に　Man=3.33e-5  ==> 0 にする
                self.Man = 0
            else:
                tmp = He_abs/self.hys_a
                tmp2 = 1/np.tanh(tmp)-1/tmp
                self.Man = self.Ms * tmp2
                if self.He < 0:
                    self.Man = -self.Man
        else:
            if self.He ==0:
                self.Man = 0
            else:
                self.Man = self.Ms * (1/np.tanh(He_abs/self.hys_a)-self.hys_a/He_abs)*self.He/He_abs

        dMan = self.Man - self.Old_Man

        X = self.kp_b*(self.Man -self.Old_M)
        if USE_JA_VER_3_1:
            if X ==0:
                dM = 0.0
            else:
                try:
                    XdHe = X*dHe
                except:
                    print('[!!! Overflow !!! in X*dHe] X=',X,'dHe=',dHe)
                    return False, 0,0,0,0
                
                tmp = max(XdHe, 0)
                if X < 0:
                    tmp = -tmp
                dM = tmp +self.cr*dMan
        else:
            if X ==0:
                dM = 0.0
            else:
                dM = max(X*dHe, 0)*X/abs(X)+self.cr*dMan

        self.M = self.Old_M + dM

        self.B = self.u0*(self.H + self.M)

        self.dB = (self.B - self.Old_B)*Fs

 
        self.Phi = self.B*self.Sd
  
        self.Old_H   = self.H
        self.Old_M   = self.M
        self.Old_B   = self.B
        self.Old_He  = self.He
        self.Old_Man = self.Man

        self.Old_Phi = self.Phi

        return self.Phi


    def getHysteresis2(self, Current,Fs):

        NI =  Current*self.Ns*self.Cal_Gain1 # scaling

        self.H =  NI/self.Ld
        self.B = self.u0*(self.H + self.M)
        self.Phi = self.B*self.Sd
  
        return self.Phi

