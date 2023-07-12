import numpy as np
from functools import lru_cache
from scipy.optimize import newton
import warnings
warnings.simplefilter('error')

#
# lib_WDF_Triode
#

#WDF_Triode_version = "1.0.5.2021.08.24"
WDF_Triode_version = "2.1.0.2023.4.6 modified from 1.0.5.2021.08.24"
WDF_Triode_version = "2.1.1.2023.4.16 modified from 1.0.5.2021.08.24"
# add 'warnings L4-5'
WDF_Triode_version = "2.1.2.2023.4.18 modified from 1.0.5.2021.08.24"   # add [Debug Input]
WDF_Triode_version = "2.1.3.2023.5.18 modified from 1.0.5.2021.08.24"   # tol: 1-e4 ==> 1-e5
WDF_Triode_version = "2.2.0.2023.6.29 modified from 1.0.5.2021.08.24"   # add triodeNL_core

def get_WDF_Triode_version():
    return WDF_Triode_version


class Triode_model:

    def __init__(self, G, muc, alpha, Ego, Cgp =0, Cgk=0, Cpk = 0):
 
        self.G = G
        self.muc = muc
        self.alpha = alpha
        self.Ego = Ego

        self.tr_a = 1.0/(1.0-alpha)
        self.tr_b = 1.5 - self.tr_a
        self.tr_c = 3.0*alpha -1.0
        self.G_p  = G*(self.tr_c*self.tr_a/3)**self.tr_b
        self.mum = self.tr_a/1.5*muc

        self.Ig_ratio =0.5/(1+1/self.mum)**1.5
        self.G_lim = self.G_p*(1+1/self.mum)**1.5

        self.gcf = self.Ig_ratio * self.G_lim
        self.igcf = (1-self.Ig_ratio) * self.G_lim

        self.Cgp = Cgp
        self.Cgk = Cgk
        self.Cpk = Cpk

        self.dVgk = 0
        self.dVpk = 0

        self.plusIg_cnt = 0 # debug


    def getIp(self, Vgk, Vpk):

        Vgg = Vgk+self.Ego
        
        if Vpk < 0:
            Vpk1 = 0
        else:
            Vpk1 = Vpk

        if Vgg+Vpk1/self.muc < 0:
            VggVpk1 = 0
        else:
            VggVpk1 = Vgg+Vpk1/self.muc

        M1  = ((self.tr_c/2/self.muc)*Vpk1+1e-10)**self.tr_b
        M2  = (1.5/self.tr_a*VggVpk1+1e-10)**self.tr_a

        if Vgg+Vpk1/self.mum < 0:
            VggVpk2 = 0
        else:
            VggVpk2 = Vgg+Vpk1/self.mum


        if Vgg > 0:
            Ik = self.G_p*(VggVpk2+1e-10)**1.5
        else:
            Ik = self.G*M1*M2


        gcf = self.gcf

        if Vgk < 0:
            Vgk = 0

        if Vgk > 0:
           Ig = gcf*Vgk**1.5*((Vgk/ (Vpk1+Vgk) )* 1.2 + 0.4)
        else:
           Ig = 0


        Ip_tmp1 = Ik -Ig -(self.igcf*Vpk1)**1.5
        if Ip_tmp1 < 0:
            Ip_tmp1 = 0
        
        Ip_tmp2 = Ik-Ig-Ip_tmp1

        if Ip_tmp2 < 0:
            Ip_tmp2 = 0

        Ip = Ip_tmp2+1e-10*Vpk


        return Ip



    def getIg(self, Vgk, Vpk): 

        gcf = self.gcf

        if Vgk < 0:
            Vgk = 0
        else:
            self.plusIg_cnt = self.plusIg_cnt +1 # debug
            #if self.plusIg_cnt == 1:
            #    print("Vgk=",Vgk)

        if Vpk < 0:
            Vpk = 0

        if Vgk > 0:
           Ig = gcf*Vgk**1.5*((Vgk/ (Vpk+Vgk) )* 1.2 + 0.4)
        else:
           Ig = 0

                 
        return Ig


    def getMu(self, Vgk, Vpk): 

        Vgg = Vgk+self.Ego

        Ig = self.getIg(Vgk, Vpk)

        Ip = self.getIp(Vgk, Vpk)

        if Ip <= 0:
            mu = 0
        else:
            if Vgg <= 0:
                mu = 1/((3-3*self.alpha)/2/self.muc + (1-3*self.alpha)/2 * Vgg/Vpk)
            else:
                mu = self.mum

        return mu

    def triodeNL(self, ap, R0, ag, Rg, ak, Rk, Vgk, Vpk, Fs):
        result, bg, bp, bk, Vgk, Vpk, plusIg_cnt = self.triodeNL_core(ap, R0, ag, Rg, ak, Rk, Vgk, Vpk, Fs)
        return bg, bp, bk, Vgk, Vpk

    def triodeNL_core(self, ap, R0, ag, Rg, ak, Rk, Vgk, Vpk, Fs):

        DISPLAY_PARAM = False   # for Debug 2023.4.4
        debug_input = [ap, R0, ag, Rg, ak, Rk, Vgk, Vpk, Fs]

        T = 1/Fs


        Old_Vgk = Vgk
        Old_Vpk = Vpk

        Icgk = self.Cgk*self.dVgk
        Icpk = self.Cpk*self.dVpk

        Ipg = self.Cgp*(self.dVpk - self.dVgk)
        Igp = -Ipg

        Icgk = Icgk - Ipg
        Icpk = Icpk - Igp

        agk = ag - ak
        apk = ap - ak

        self.plusIg_cnt = 0 # debug
        def fg(Vgkx, Vpkx):
            return Vgkx + Rg * (self.getIg(Vgkx, Vpkx)+Icgk)-agk

        fg_tol = 1e-5   #1e-4
        fg0 = fg(0,Vpk)
        if fg0 > 0:
            fg_start = -1
        else:
            fg_start = 1
        #fg_start = Vgk
        #print('Vgk,Vpk,Rg,Icgk,agk,gcf,',Vgk,',',Vpk,',',Rg,',',Icgk,',',agk,self.gcf)
        debug_Vgk = [Vgk,Vpk,Rg,Icgk,agk,self.gcf,fg_start]
        try:
            Vgk = newton(fg, fg_start,  args=(Vpk,), tol=fg_tol, maxiter=100)
            #Vgk = newton(fg, Vgk,  args=(Vpk,), tol=1e-4, maxiter=100)
            y = fg(Vgk,Vpk)
            if y > 1e-4 or y < -1e-4:
                print('[Wrong result] Vgk=',Vgk,'y=',y,'fg0=',fg0,'fg_start=',fg_start)
                return False, 0.0, 0.0, 0.0, 0.0, 0.0, 0
        except:
            print('[Debug Input] ap, R0, ag, Rg, ak, Rk, Vgk, Vpk, Fs',debug_input)
            debug_Vgk.append(self.plusIg_cnt)
            print("[Debug Vgk] VgkIN,Vpk,Rg,Icgk,agk,gcf,fg_start,plusIg_cnt:",debug_Vgk)
            print("Error: Vgk = newton(fg, Vgk,  args=(Vpk,), tol=1e-4, maxiter=100)")
            return False, 0.0, 0.0, 0.0, 0.0, 0.0, 0


        #if self.plusIg_cnt >0: # debug
        #    print("fg: plusIg_cnt=",self.plusIg_cnt)

        def fp(Vpkx, Vgkx): 
           return Vpkx + R0 * (self.getIp(Vgkx, Vpkx) + Icpk) - apk

        fp_tol = 1e-5   #1e-4
        fp0 = fp(0,Vgk)
        if fp0 > 0:
            fp_start = -1
        else:
            fp_start = 1
        debug_Vpk = [Vgk,Vpk,R0,Icpk,apk,self.gcf,fg_start,self.Ego,self.muc,self.tr_a,self.tr_b,self.tr_c,self.mum,self.G_p,self.G,self.gcf,self.igcf]
        if DISPLAY_PARAM:
            print('Vgk,Vpk,R0,Icpk,apk,',Vgk,',',Vpk,',',R0,',',Icpk,',',apk)
            print('Egp,muc,tr_a,tr_b,tr_c,mum,G_p,G,gcf,igcf,',self.Ego,self.muc,self.tr_a,self.tr_b,self.tr_c,self.mum,self.G_p,self.G,self.gcf,self.igcf)
        try:
            Vpk =  newton(fp, fp_start,  args=(Vgk,), tol=fp_tol, maxiter=100)
            #Vpk =  newton(fp, Vpk,  args=(Vgk,), tol=1e-4, maxiter=100)
            y = fp(Vpk,Vgk)
            if y > 1e-4 or y < -1e-4:
                print('[Wrong result] Vpk=',Vpk,'y=',y,'fp0=',fp0,'fp_start=',fp_start)
                return False, 0.0, 0.0, 0.0, 0.0, 0.0, 0
        except:
            print('[Debug Input] ap, R0, ag, Rg, ak, Rk, Vgk, Vpk, Fs',debug_input)
            debug_Vgk.append(self.plusIg_cnt)
            print("[Debug Vgk] VgkIN,Vpk,Rg,Icgk,agk,gcf,fg_start,plusIg_cnt:",debug_Vgk)
            print('[Debug Vpk] Vgk,VpkIN,R0,Icpk,apk,gcf,fg_start,Ego,muc,tr_a,tr_b,tr_c,mum,G_p,G,gcf,igcf:',debug_Vpk)
            print("Error: Vpk = newton(fp, Vpk,  args=(Vgk,), tol=1e-4, maxiter=100)")
            return False, 0.0, 0.0, 0.0, 0.0, 0.0, 0

        Ig = self.getIg(Vgk, Vpk) + Icgk
        Ip = self.getIp(Vgk, Vpk) + Icpk


        Ik = Ig + Ip

        bk =  ak + 2*Rk*Ik

        bp = Vpk - R0* Ip + bk
        bg = Vgk - Rg* Ig + bk


        self.dVgk = (Vgk - Old_Vgk)/T
        self.dVpk = (Vpk - Old_Vpk)/T


        return True, bg, bp, bk, Vgk, Vpk, self.plusIg_cnt

    def triodeNL0Ig(self, ap, R0, ag, Rg, ak, Rk, Vgk, Vpk, Fs):

        DISPLAY_PARAM = False   # for Debug 2023.4.4
        debug_input = [ap, R0, ag, Rg, ak, Rk, Vgk, Vpk, Fs]

        T = 1/Fs


        Old_Vgk = Vgk
        Old_Vpk = Vpk

        Icgk = self.Cgk*self.dVgk
        Icpk = self.Cpk*self.dVpk

        Ipg = self.Cgp*(self.dVpk - self.dVgk)
        Igp = -Ipg

        Icgk = Icgk - Ipg
        Icpk = Icpk - Igp

        agk = ag - ak
        apk = ap - ak

        self.plusIg_cnt = 0 # debug
        def fg(Vgkx, Vpkx):
            return Vgkx + Rg * (Icgk)-agk
            #return Vgkx + Rg * (self.getIg(Vgkx, Vpkx)+Icgk)-agk

        fg0 = fg(0,Vpk)
        if fg0 > 0:
            fg_start = -1
        else:
            fg_start = 1
        #fg_start = Vgk
        #print('Vgk,Vpk,Rg,Icgk,agk,gcf,',Vgk,',',Vpk,',',Rg,',',Icgk,',',agk,self.gcf)
        debug_Vgk = [Vgk,Vpk,Rg,Icgk,agk,self.gcf,fg_start]
        try:
            Vgk = newton(fg, fg_start,  args=(Vpk,), tol=1e-4, maxiter=100)
            #Vgk = newton(fg, Vgk,  args=(Vpk,), tol=1e-4, maxiter=100)
            y = fg(Vgk,Vpk)
            if y > 1e-4 or y < -1e-4:
                print('[Wrong result on Vgk] Vgk=',Vgk,'y=',y)
                print('[Debug Input] ap, R0, ag, Rg, ak, Rk, Vgk, Vpk, Fs',debug_input)
                debug_Vgk.append(self.plusIg_cnt)
                print("[Debug Vgk] VgkIN,Vpk,Rg,Icgk,agk,gcf,fg_start,plusIg_cnt:",debug_Vgk)
                return False, 0.0, 0.0, 0.0, 0.0, 0.0, 0
        except:
            print('[Debug Input] ap, R0, ag, Rg, ak, Rk, Vgk, Vpk, Fs',debug_input)
            debug_Vgk.append(self.plusIg_cnt)
            print("[Debug Vgk] VgkIN,Vpk,Rg,Icgk,agk,gcf,fg_start,plusIg_cnt:",debug_Vgk)
            print("Error: Vgk = newton(fg, Vgk,  args=(Vpk,), tol=1e-4, maxiter=100)")
            return False, 0.0, 0.0, 0.0, 0.0, 0.0, 0


        #if self.plusIg_cnt >0: # debug
        #    print("fg: plusIg_cnt=",self.plusIg_cnt)

        def fp(Vpkx, Vgkx): 
           return Vpkx + R0 * (self.getIp(Vgkx, Vpkx) + Icpk) - apk

        fp0 = fp(0,Vgk)
        if fp0 > 0:
            fp_start = -1
        else:
            fp_start = 1
        debug_Vpk = [Vgk,Vpk,R0,Icpk,apk,self.gcf,fg_start,self.Ego,self.muc,self.tr_a,self.tr_b,self.tr_c,self.mum,self.G_p,self.G,self.gcf,self.igcf]
        if DISPLAY_PARAM:
            print('Vgk,Vpk,R0,Icpk,apk,',Vgk,',',Vpk,',',R0,',',Icpk,',',apk)
            print('Egp,muc,tr_a,tr_b,tr_c,mum,G_p,G,gcf,igcf,',self.Ego,self.muc,self.tr_a,self.tr_b,self.tr_c,self.mum,self.G_p,self.G,self.gcf,self.igcf)
        try:
            Vpk =  newton(fp, fp_start,  args=(Vgk,), tol=1e-4, maxiter=100)
            #Vpk =  newton(fp, Vpk,  args=(Vgk,), tol=1e-4, maxiter=100)
            y = fp(Vpk,Vgk)
            if y > 1e-4 or y < -1e-4:
                print('[Wrong result on Vpk] Vpk=',Vpk,'y=',y)
                print('[Debug Input] ap, R0, ag, Rg, ak, Rk, Vgk, Vpk, Fs',debug_input)
                debug_Vgk.append(self.plusIg_cnt)
                print("[Debug Vgk] VgkIN,Vpk,Rg,Icgk,agk,gcf,fg_start,plusIg_cnt:",debug_Vgk)
                print('[Debug Vpk] Vgk,VpkIN,R0,Icpk,apk,gcf,fg_start,Ego,muc,tr_a,tr_b,tr_c,mum,G_p,G,gcf,igcf:',debug_Vpk)
                return False, 0.0, 0.0, 0.0, 0.0, 0.0, 0
        except:
            print('[Debug Input] ap, R0, ag, Rg, ak, Rk, Vgk, Vpk, Fs',debug_input)
            debug_Vgk.append(self.plusIg_cnt)
            print("[Debug Vgk] VgkIN,Vpk,Rg,Icgk,agk,gcf,fg_start,plusIg_cnt:",debug_Vgk)
            print('[Debug Vpk] Vgk,VpkIN,R0,Icpk,apk,gcf,fg_start,Ego,muc,tr_a,tr_b,tr_c,mum,G_p,G,gcf,igcf:',debug_Vpk)
            print("Error: Vpk = newton(fp, Vpk,  args=(Vgk,), tol=1e-4, maxiter=100)")
            return False, 0.0, 0.0, 0.0, 0.0, 0.0, 0

        Ig = self.getIg(Vgk, Vpk) + Icgk
        Ip = self.getIp(Vgk, Vpk) + Icpk


        Ik = Ig + Ip

        bk =  ak + 2*Rk*Ik

        bp = Vpk - R0* Ip + bk
        bg = Vgk - Rg* Ig + bk


        self.dVgk = (Vgk - Old_Vgk)/T
        self.dVpk = (Vpk - Old_Vpk)/T


        return True, bg, bp, bk, Vgk, Vpk, self.plusIg_cnt

    def triodeNL0IgOrg(self, ap, R0, ag, Rg, ak, Rk, Vgk, Vpk, Fs):

        T = 1/Fs

        Old_Vgk = Vgk
        Old_Vpk = Vpk

        Icgk = self.Cgk*self.dVgk
        Icpk = self.Cpk*self.dVpk

        Ipg = self.Cgp*(self.dVpk - self.dVgk)
        Igp = -Ipg

        Icgk = Icgk - Ipg
        Icpk = Icpk - Igp

        agk = ag - ak
        apk = ap - ak

        self.plusIg_cnt = 0 # debug
        def fg(Vgkx, Vpkx):
            return Vgkx + Rg * (Icgk)-agk
            #return Vgkx + Rg * (self.getIg(Vgkx, Vpkx)+Icgk)-agk

        try:
            Vgk = newton(fg, Vgk,  args=(Vpk,), tol=1e-4, maxiter=100)
        except:
            print("Error: Vgk = newton(fg, Vgk,  args=(Vpk,), tol=1e-4, maxiter=100)")
            return False, 0.0, 0.0, 0.0, 0.0, 0.0, 0


        #if self.plusIg_cnt >0: # debug
        #    print("fg: plusIg_cnt=",self.plusIg_cnt)

        def fp(Vpkx, Vgkx): 
           return Vpkx + R0 * (self.getIp(Vgkx, Vpkx) + Icpk) - apk

        Vpk =  newton(fp, Vpk,  args=(Vgk,), tol=1e-4, maxiter=100)

        Ig = self.getIg(Vgk, Vpk) + Icgk
        Ip = self.getIp(Vgk, Vpk) + Icpk


        Ik = Ig + Ip

        bk =  ak + 2*Rk*Ik

        bp = Vpk - R0* Ip + bk
        bg = Vgk - Rg* Ig + bk


        self.dVgk = (Vgk - Old_Vgk)/T
        self.dVpk = (Vpk - Old_Vpk)/T

        #if self.plusIg_cnt > 0:
        #    print('Called')

        return True, bg, bp, bk, Vgk, Vpk, self.plusIg_cnt


