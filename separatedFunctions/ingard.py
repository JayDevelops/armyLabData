import pandas as pd  # Install pandas in python or use Anaconda environment
import math
from DetectionScripts import data_dict
from DetectionScripts import ADM_Functions

def Ingard(freq):
    Hs = data_dict.det_cons['source_height_det'] 
    Hr = data_dict.det_cons['listener_height_det'] 
    R = data_dict.given_cons['m_measure_distance']
    celsius_degrees_T1 = data_dict.det_cons['celsius_degrees_det']
    kelvin_Too = 273.15
    Co = 331.32
    Rhoa = 0.747
    Rhod = 0.747
    L = 1.1
    Pi = 4 * math.atan(1)
    Em2 = data_dict.det_cons['em2_det']
    Sigma = data_dict.meas_cons['sigma_meas']
    Log10Div10 = 0.230258509
    TenDivLog10 = 1 / Log10Div10
    ground_effect = []
  
    HspHr = Hs + Hr
    HsmHr = Hs - Hr
    R1 = math.sqrt(R * R + HsmHr * HsmHr)
    R2 = math.sqrt(R * R + HspHr * HspHr)
    R12 = R1 / R2
    Dr = R2 - R1
    SinTheta1 = HspHr / R2
    CosTheta1 = R / R2
    TanTheta1 = HspHr / R
    Tk = kelvin_Too + celsius_degrees_T1
    Cs = Co * math.sqrt(Tk / kelvin_Too)
    Bef = math.exp(math.log(2) / 6)
    Mu = (Bef - 1 / Bef) / 2
    Eta = (Bef + 1 / Bef) / 2
    I1c = math.sqrt(Pi) * Em2 * R1 * L

    for i in range(24):
        K1 = 2 * Pi * freq[i] / Cs
        K1dr = K1 * Dr
        Del = R1 / (K1 * L * L)
        Omega = math.sqrt(1 + 1 / (Del * Del)) - 1
        Dom1 = Del * Omega
        Dom2 = Del * math.sqrt(2 * Omega)
        At = math.atan(Dom1 / (1 - Dom2)) - math.atan(Dom1 / (1 + Dom2))
        L1 = 0.5 * Dom1 * math.log((1 + Dom2) / abs(1 - Dom2)) + At
        I1 = I1c * K1 * K1
        I2 = 0.5 * I1 * L1 / ((Dom1 + Del) * Dom2)
        Sd2 = 0.5 * (I1 + I2)
        X = 0.5 * (I1 - I2)
        if(X > 1):
            Ea2 = 0.27 * math.exp(math.log(X) * 0.33)
        else:
            Ea2 = X / (1 + 11 * X / 4)

        #Functions we need for this portion: ComplexDiv, ComplexMul, Polar
        Lnfos = math.log(freq[i] / Sigma)
        Rrc = 1 + 9.08 * math.exp(-0.75 * Lnfos)
        Xrc = 11.9 * math.exp(-0.73 * Lnfos)
        Rzr, Rzi = ComplexDiv(1, 0, Rrc, Xrc)
        Rpr, Rpi = ComplexDiv(SinTheta1 - Rzr, -Rzi, SinTheta1 + Rzr, Rzi)
        Srzr, Srzi = ComplexMul(SinTheta1 + Rzr, Rzi, SinTheta1 + Rzr, Rzi)
        Czr, Czi = ComplexDiv(Srzr, Srzi, 1 + SinTheta1 * Rzr, SinTheta1 * Rzi)
        Wr, Wi = ComplexMul(0, 0.5 * K1 * R2, Czr, Czi)
        Wm, Wp = Polar(Wr, Wi) #Means Wm = Sqr((Wr * Wr) + (Wi * Wi)) // and Wp = math.atan2(Wr, Wi)
    
        if Wm < 6:
            Intg = 3
            Fact = 1
            Wsr = Wr
            Wsi = Wi
            W1r = Wr
            W1i = Wi

            for J in range(1, 12):
                Fact = Fact * J
                Cons = Fact * Intg
                W2r, W2i = ComplexMul(Wr, Wi, W1r, W1i)
                Wsr = Wsr + W2r / Cons
                Wsi = Wsi + W2i / Cons
                W1r = W2r
                W1i = W2i
                Intg = Intg + 2

            Ewr = math.exp(-Wr)
            Ewrr = Ewr * math.cos(Wi)
            Ewri = -Ewr * math.sin(Wi)
            Whm = math.sqrt(Pi * Wm)
            Whp = Wp / 2
            Whr = Whm * math.cos(Whp + Pi / 2)
            Whi = Whm * math.sin(Whp + Pi / 2)
            Ws1r, Ws1i = ComplexMul(Ewrr, Ewri, Whr - 2 * Wsr, Whi - 2 * Wsi)
            Fr = 1 + Ws1r
            Fi = Ws1i
        else:
            W1r, W1i = ComplexDiv(1, 0, 2 * Wr, 2 * Wi)
            W2r, W2i = ComplexMul(W1r, W1i, W1r, W1i)
            Wsr = W1r + 3 * W2r
            Wsi = W1i + 3 * W2i
            W3r, W3i = ComplexMul(W1r, W1i, W2r, W2i)
            Wsr = Wsr + 15 * W3r
            Wsi = Wsi + 15 * W3i
            Fr = -Wsr
            Fi = -Wsi
        
        Q1r, Q1i = ComplexMul(Fr, Fi, 1 - Rpr, -Rpi)
        Fm, Fp = Polar(Fr, Fi)
        Qr = Rpr + Q1r
        Qi = Rpi + Q1i
        Qm, Qp = Polar(Qr, Qi)
        Mkdr = Mu * K1dr
        Sinc = math.sin(Mkdr) / Mkdr
        Cost = math.cos(Eta * K1dr + Qp)
        Qmr = Qm * R12
        I3 = (1 + Ea2) * (1 + Qmr * Qmr) + 2 * Qmr * (1 + Ea2 * Rhoa) * Cost * Sinc * math.exp(-Sd2 * (1 - Rhod))
        ground_effect.append(TenDivLog10 * math.log(I3))
    
    return ground_effect
