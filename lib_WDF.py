#
# lib_WDF
# WUが末端からルート, WDがルートから末端
# WD : 入射波(素子に到来する波成分), WU : 反射波(素子から離れる波成分)

WDF_version = "1.2.0.2021.08.17"

def get_WDF_version():
    return WDF_version

class WDF(object):

    def __init__(self, Rp=0):
        self.Rp = Rp
        self.WU = 0
        self.WD = 0

    def Voltage(self):
        return (self.WU + self.WD)/2 # WD : 入射波(素子に到来する波成分), WU : 反射波(素子から離れる波成分)

    def Current(self):
        return (self.WD - self.WU)/( 2* self.Rp )


class OnePort(WDF):

    def __init__(self, Rp):
        self.Rp = Rp
        self.WU = 0
        self.WD = 0


    def WaveDown(self,val):
        self.WD = val;


class Adaptor(WDF):

   def __init__(self):
        self.PortLeft = 0
        self.PortRight = 0

class Resistor(OnePort):

    def __init__(self, Rp):
        self.Rp = Rp
        self.WU = 0
        self.WD = 0

    def WaveUp(self):
        self.WU = 0 # 反射波なし, b[n] = 0
        return self.WU


class Capacitor(OnePort):

    def __init__(self, Rp, Fs):
        T = 1/Fs
        self.Rp = T/(2*Rp)
        self.State = 0
        self.WU = 0
        self.WD = 0

    def WaveUp(self):
        self.WU = self.State # b[n] = a[n-1]
        return self.WU

    def WaveDown(self,val):
        self.WD = val
        self.State = val

    def setWD(self,val):
        self.WD = val
        self.State = val

    def getState(self):
        return self.State


class Inductor(OnePort):

    def __init__(self, Rp, Fs):
        T = 1/Fs
        self.Rp = 2*Rp/T
        self.State = 0
        self.WU = 0
        self.WD = 0

    def WaveUp(self): # get the up-going wave
        self.WU = -self.State # b[n] = -a[n-1]
        return self.WU

    def WaveDown(self,val):
        self.WD = val
        self.State = val

    def setWD(self,val):
        self.WD = val
        self.State = val

    def getState(self):
        return self.State


class Diode(OnePort):

    def __init__(self, Rp, Vp=0):
        self.Rp = Rp
        self.Rpp = Rp
        self.WU = 0
        self.WD = 0

        self.Vrd = 0
        self.Vp = 0

    def WaveUp(self):
        kd = 0.005
        ud = 2.1

#        if self.Vrd > 10**(-60):
#            Rd = 1/kd * self.Vrd**(1-ud)
#        else:
#            Rd = 1/kd *10**(-60*(1-ud))

        if self.Vrd >= self.Vp:
            Rd = self.Rpp
        else:
            Rd = 100e9

        self.Rp = Rd
        self.WU = 0

        self.Vrd = (self.WD + self.WU)/2

        return 0


class VoltageSource(OnePort): # 11/25 米田さんより挙動がおかしいかも要確認, 内部抵抗無し電圧源

    def __init__(self, E, Rp):
        self.E = E
        self.Rp = Rp
        self.WU = 0
        self.WD = 0

    def WaveUp(self):
        self.WU = 2*self.E-self.WD
        return self.WU


class CurrentSource(OnePort): # 11/25 米田さんより挙動がおかしいかも要確認, 内部抵抗無し電圧源

    def __init__(self, Is, Rp):
        self.Is = Is
        self.Rp = Rp
        self.WU = 0
        self.WD = 0

    def WaveUp(self):
        self.WU = self.WD + 2*self.Is*self.Rp
        return self.WU


class TerminatedVs(OnePort): # 内部抵抗付き電圧源

    def __init__(self, E, Rp):
        self.E = E
        self.Rp = Rp
        self.WU = 0
        self.WD = 0

    def WaveUp(self):
        self.WU = self.E # b[n] = Vs
        return self.WU


class TerminatedIs(OnePort): # 内部抵抗付き電流源

    def __init__(self, Is, Rp):
        self.Is = Is
        self.Rp = Rp
        self.WU = 0
        self.WD = 0

    def WaveUp(self):
        self.WU = self.Is*self.Rp
        return self.WU


class OpenCircuit(OnePort):

    def __init__(self, Rp):
        self.Rp = Rp
        self.WU = 0
        self.WD = 0

    def WaveUp(self):
        self.WU = self.WD
        return self.WU


class ShortCircuit(OnePort):

    def __init__(self, Rp):
        self.Rp = Rp
        self.WU = 0
        self.WD = 0

    def WaveUp(self):
        self.WU = -self.WD
        return self.WU

# 直列接続
class Series(Adaptor):

    def __init__(self, PortLeft,PortRight):

        self.WD = 0
        self.WU = 0

        self.PortLeft = PortLeft
        self.PortRight = PortRight
        self.Rp = PortLeft.Rp+PortRight.Rp # 直列抵抗


    def WaveUp(self):
        self.WU = -(self.PortLeft.WaveUp()+self.PortRight.WaveUp()) # 反射波は-になる, b3 = -(a1 + a2)
        return self.WU


    def WaveDown(self,WaveFromParent):
        self.WD = WaveFromParent
  
#        lrW = (WaveFromParent+self.PortLeft.WU+self.PortRight.WU)
        lrW = (WaveFromParent - self.WU) # a3 - b3
        left = self.PortLeft.WU - ((self.PortLeft.Rp/self.Rp) * lrW)    # b1 = a1 - γ1s * (a3 - b3)
        right = self.PortRight.WU - ((self.PortRight.Rp/self.Rp) * lrW) # b2 = a2 - γ2s * (a3 - b3)
             
        self.PortLeft.WaveDown(left)
        self.PortRight.WaveDown(right)

    def UpdateRp(self):
        if issubclass(type(self.PortLeft), Adaptor) == True:
            self.PortLeft.UpdateRp()

        if issubclass(type(self.PortRight), Adaptor) == True:
            self.PortRight.UpdateRp()

        self.Rp = self.PortLeft.Rp+self.PortRight.Rp


# 並列接続
class Parallel(Adaptor):

    def __init__(self, PortLeft,PortRight):
        self.WD = 0
        self.WU = 0
        self.a1 = 0
        self.a2 = 0
 
        self.PortLeft = PortLeft
        self.PortRight = PortRight
            
        R1 = PortLeft.Rp
        R2 = PortRight.Rp
        R = (R1 * R2)/(R1 + R2) # 並列接続
        self.Rp = R

        Gu = 1.0 / R
        Gl = 1.0 / R1
        Gr = 1.0 / R2
        self.dl = 2.0 * Gl / ( Gu + Gl + Gr )
        self.dr = 1.0 - self.dl
                       
    def WaveUp(self):
        self.a1 = self.PortLeft.WaveUp()
        self.a2 = self.PortRight.WaveUp()

        self.WU = self.dl * self.a1 + self.dr * self.a2 # b3 = γ1p*a1 + γ2p*a2
        return self.WU

    def WaveDown(self,WaveFromParent):
        self.WD = WaveFromParent
            
        a3 = WaveFromParent

#        left =  ( self.dl - 1 ) * self.a1 + self.dr * self.a2 + a3
#        right = self.dl * self.a1 + ( self.dr - 1 ) * self.a2 + a3

        left = self.WU + a3 - self.a1  # b1 = b3 + a3 - a1
        right = self.WU + a3 - self.a2 # b2 = b3 + a3 - a2


        self.PortLeft.WaveDown(left)
        self.PortRight.WaveDown(right)


    def UpdateRp(self):
        if issubclass(type(self.PortLeft), Adaptor) == True:
            self.PortLeft.UpdateRp()

        if issubclass(type(self.PortRight), Adaptor) == True:
            self.PortRight.UpdateRp()

        self.Rp = self.PortLeft.Rp+self.PortRight.Rp

        R1 = self.PortLeft.Rp
        R2 = self.PortRight.Rp
        R = (R1 * R2)/(R1 + R2)
        self.Rp = R

        Gu = 1.0 / R
        Gl = 1.0 / R1
        Gr = 1.0 / R2
        self.dl = 2.0 * Gl / ( Gu + Gl + Gr )
        self.dr = 1.0 - self.dl


class IdealTransformer(Adaptor):

    def __init__(self, PortRight, N):
        self.WD = 0
        self.WU = 0
        self.A2 = 0

        self.PortRight = PortRight
        self.N = N
        self.Rp = N**2*PortRight.Rp

    def WaveUp(self):
        A1 = self.PortRight.WaveUp()
        self.WU = self.N*A1
        return self.WU

    def WaveDown(self, WaveFromParent):
        self.A2 = WaveFromParent
        self.WD = (1/self.N)*self.A2
        self.PortRight.WaveDown(self.WD)

    def UpdateRp(self):
        if issubclass(type(self.PortRight), Adaptor) == True:
            self.PortRight.UpdateRp()

        self.Rp = self.N**2*self.PortRight.Rp

    def Voltage(self):
        return (self.WU + self.A2)/2

    def Current(self):
        return (self.A2 - self.WU)/( 2* self.Rp )


class PolarityInverter(Adaptor):

    def __init__(self, PortRight):
        self.WD = 0
        self.WU = 0
        self.A2 = 0

        self.PortRight = PortRight
        self.Rp = PortRight.Rp

    def WaveUp(self):
        A1 = self.PortRight.WaveUp()
        self.WU = -A1
        return self.WU

    def WaveDown(self, WaveFromParent):
        self.A2 = WaveFromParent
        self.WD = -self.A2
        self.PortRight.WaveDown(self.WD)

    def UpdateRp(self):
        if issubclass(type(self.PortRight), Adaptor) == True:
            self.PortRight.UpdateRp()

        self.Rp = self.PortRight.Rp

    def Voltage(self):
        return (self.WU + self.A2)/2

    def Current(self):
        return (self.A2 - self.WU)/( 2* self.Rp )

