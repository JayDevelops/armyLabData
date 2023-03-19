' Option Explicit
Option Base 0
Public Pi, Log10, Log10Div10, TenDivLog10
Public I As Integer, J As Integer, K As Integer
Public Done As Integer, LastCol As Integer, ShapeNo As Integer
Public celsius_degrees_det, relative_humid_percent, R, Hs, Hr, zerDegCelsius_inKelvin, Em2, Em2Det, Em2Ref, Sigma, SigmaRef, SigmaDet, Co, Cs
Public source_height_meas_cond, mic_height_meas_cond, observer_efficiency, calculate_d, hit_prob, false_alarm_rate, source_height_det, listener_height_det, Hba, Iwthr1, M2, D5, D6, precision_fraction, B1, WindSpeed
Public X1, X2, X3, distance_from_source, b8, Ka, distance_frmSource_to_fol, extent_of_foliage, avg_leaf_width, leaf_area_per_unit, Trg, Bkg, Hth, Ai, PCor, Pc1, Pc2, Pc3, Pc4
Public Ea, dBa


Public m_meas_distance_D3 As Single     ' m Meas Dist
Public detection_dist As Single     ' Detection Dist
Public barrier_number As Integer     ' Barrier number
Public foliated_zone_nums As Integer     ' Foliated-zone numbers
Public target_spectrum(26) As Single  ' Target spectrum
Public max_modBackground_thresh(26) As Single  ' MAX(mod_background_noise,hearing_threshold)
Public mod_background_noise(26) As Single ' Modified Background Noise
Public mod_background_noise_2(26) As Single ' Modified Background Noise
Public prop_loss_cum(26) As Single  ' Propagation Loss cumulative
Public freq(26) As Integer ' Freq
Public hearing_threshold(26) As Single  ' Hearing Threshold
Public background_noise_spec(26) As Single  ' Background noise spectrum
Public listener_spectrum(26) As Single  ' Listener spectrum
Public prop_loss_indiv(26) As Single  ' Propagation Loss individual
Public prop_loss_cum_4graph(26) As Single  ' Propagation Levels cumulative
Public A_weight_levels(26) As Single  ' A weight levels
Public signal_noise_ratio(26) As Single  ' signal_noise_ratio in AI calculation
Public AiWt(24) As Single ' AI weights
Public detection_band(24) As Single        'Detection band
Public atmos_absorption(24) As Single    'Atmosp Absorption
Public ground_effect(26) As Single        'Ground Effect
Public AtmAbsRef(26) As Single 'Reference Atmosp Absorption during measurement
Public ground_effect_reference(26) As Single     'Reference Ground Effect during measurement
Public barrier_attenuation(26) As Single        'Barrier atten
Public foliage_attenuation(26) As Single        'Foliage tten
Public wind_refrac_loss(26) As Single      'wind_refrac_loss refraction loss (from early version)
Public B3(23) As Single
Public one_third_oct_band_weight(23, 4) As Single     ' 1/3 OB Weightings for masking bands
Public Threshold As Variant
Sub InverseDistance()
  Pinv = 2 * TenDivLog10 * Log(m_meas_distance / detection_dist) ' Inverse distance loss from m_meas_distan to detection_dist
  For I = 0 To 23
    prop_loss_indiv(I) = Pinv
    If prop_loss_indiv(I) = 0 Then prop_loss_indiv(I) = -0.001
    prop_loss_cum(I) = prop_loss_indiv(I)
  Next I
End Sub
Sub AnsiHumidity()
  T = zerDegCelsius_inKelvin + celsius_degrees_det
  Pstar = 1
  Acr1 = 1.6E-08 * (1.377 * T / (T + 110.4)) / Pstar
  Thto = 2239.1 / T
  Exthto = Exp(-Thto)
  Rto = Thto / (1 - Exthto)
  Mumaxo = 0.20948 * 4 * Pi * Rto * Rto * Exthto / 35
  Thtn = 3352! / T
  Exthtn = Exp(-Thtn)
  Rtn = Thtn / (1 - Exthtn)
  Mumaxn = 0.78084 * 4 * Pi * Rtn * Rtn * Exthtn / 35
  Cs = Co * Sqr(T / zerDegCelsius_inKelvin)
  H = relative_humid_percent * Exp(Log10 * (20.5318 - 2939 / T - 2.13759744 * Log(T))) / Pstar
  Fox = (24 + 44100 * H * (0.05 + H) / (0.391 + H)) * Sqr(293 / T)    ' Is 293 used as 273 + 20degC?
  Fni = (9 + 350 * H) * (293 / T)
  calculate_d = R * 0.01
  For I = 0 To 23
    Freq2 = freq(I)
    Freq2 = Freq2 * Freq2
    Acr = Acr1 * Freq2
    X2o = Freq2 / (Fox * Fox)
    Avibo = 868.5 * Mumaxo * Fox / Cs * X2o / (1 + X2o)
    X2n = Freq2 / (Fni * Fni)
    Avibn = 868.5 * Mumaxn * Fni / Cs * X2n / (1 + X2n)
    atmos_absorption(I) = -calculate_d * (Acr + Avibo + Avibn)
  Next I
End Sub
Sub Atmosphere()
  AnsiHumidity
  For I = 0 To 23
    prop_loss_indiv(I) = atmos_absorption(I) - AtmAbsRef(I)
    If prop_loss_indiv(I) = 0 Then prop_loss_indiv(I) = -0.001
    prop_loss_cum(I) = prop_loss_cum(I) + prop_loss_indiv(I)
   Next I
End Sub
Sub Ingard()
  HspHr = Hs + Hr
  HsmHr = Hs - Hr
  R1 = Sqr(R * R + HsmHr * HsmHr)
  R2 = Sqr(R * R + HspHr * HspHr)
  R12 = R1 / R2
  Dr = R2 - R1
  SinTheta1 = HspHr / R2
  CosTheta1 = R / R2
  TanTheta1 = HspHr / R
  Tk = zerDegCelsius_inKelvin + celsius_degrees_det
  Cs = Co * Sqr(Tk / zerDegCelsius_inKelvin)
  Bef = Exp(Log(2) / 6)
  Mu = (Bef - 1 / Bef) / 2
  Eta = (Bef + 1 / Bef) / 2
  Rhoa = 0.747
  Rhod = 0.747
  L = 1.1
  I1c = Sqr(Pi) * Em2 * R1 * L
  For I = 0 To 23
    freq = freq(I)
    K1 = 2 * Pi * freq / Cs
    K1dr = K1 * Dr
    Del = R1 / (K1 * L * L)
    Omega = Sqr(1 + 1 / (Del * Del)) - 1
    Dom1 = Del * Omega
    Dom2 = Del * Sqr(2 * Omega)
    At = Atn(Dom1 / (1 - Dom2)) - Atn(Dom1 / (1 + Dom2))
    L1 = 0.5 * Dom1 * Log((1 + Dom2) / Abs(1 - Dom2)) + At
    I1 = I1c * K1 * K1
    I2 = 0.5 * I1 * L1 / ((Dom1 + Del) * Dom2)
    Sd2 = 0.5 * (I1 + I2)
    X = 0.5 * (I1 - I2)
    If X > 1 Then
      Ea2 = 0.27 * Exp(Log(X) * 0.33)
    Else
      Ea2 = X / (1 + 11 * X / 4)
    End If
    Lnfos = Log(freq / Sigma)
    Rrc = 1 + 9.08 * Exp(-0.75 * Lnfos)
    Xrc = 11.9 * Exp(-0.73 * Lnfos)
    Done = ComplexDiv(1, 0, Rrc, Xrc, Rzr, Rzi)
    Done = ComplexDiv(SinTheta1 - Rzr, -Rzi, SinTheta1 + Rzr, Rzi, Rpr, Rpi)
    Done = ComplexMul(SinTheta1 + Rzr, Rzi, SinTheta1 + Rzr, Rzi, Srzr, Srzi)
    Done = ComplexDiv(Srzr, Srzi, 1 + SinTheta1 * Rzr, SinTheta1 * Rzi, Czr, Czi)
    Done = ComplexMul(0, 0.5 * K1 * R2, Czr, Czi, Wr, Wi)
    Done = Polar(Wr, Wi, Wm, Wp)
    If Wm < 6 Then
      Intg = 3
      Fact = 1
      Wsr = Wr
      Wsi = Wi
      W1r = Wr
      W1i = Wi
      For J = 1 To 11
        Fact = Fact * J
        Cons = Fact * Intg
        Done = ComplexMul(Wr, Wi, W1r, W1i, W2r, W2i)
        Wsr = Wsr + W2r / Cons
        Wsi = Wsi + W2i / Cons
        W1r = W2r
        W1i = W2i
        Intg = Intg + 2
      Next J
      Ewr = Exp(-Wr)
      Ewrr = Ewr * Cos(Wi)
      Ewri = -Ewr * Sin(Wi)
      Whm = Sqr(Pi * Wm)
      Whp = Wp / 2
      Whr = Whm * Cos(Whp + Pi / 2)
      Whi = Whm * Sin(Whp + Pi / 2)
      Done = ComplexMul(Ewrr, Ewri, Whr - 2 * Wsr, Whi - 2 * Wsi, Ws1r, Ws1i)
      Fr = 1 + Ws1r
      Fi = Ws1i
    Else
      Done = ComplexDiv(1, 0, 2 * Wr, 2 * Wi, W1r, W1i)
      Done = ComplexMul(W1r, W1i, W1r, W1i, W2r, W2i)
      Wsr = W1r + 3 * W2r
      Wsi = W1i + 3 * W2i
      Done = ComplexMul(W1r, W1i, W2r, W2i, W3r, W3i)
      Wsr = Wsr + 15 * W3r
      Wsi = Wsi + 15 * W3i
      Fr = -Wsr
      Fi = -Wsi
    End If
    Done = ComplexMul(Fr, Fi, 1 - Rpr, -Rpi, Q1r, Q1i)
    'Done = Polar(Fr, Fi, Fm, Fp)
    Qr = Rpr + Q1r
    Qi = Rpi + Q1i
    Done = Polar(Qr, Qi, Qm, Qp)
    Mkdr = Mu * K1dr
    Sinc = Sin(Mkdr) / Mkdr
    Cost = Cos(Eta * K1dr + Qp)
    Qmr = Qm * R12
    I3 = (1 + Ea2) * (1 + Qmr * Qmr) + 2 * Qmr * (1 + Ea2 * Rhoa) * Cost * Sinc * Exp(-Sd2 * (1 - Rhod))
    ground_effect(I) = TenDivLog10 * Log(I3)
  Next I
End Sub
Sub GroundEffect()
  Hs = source_height_det
  Hr = listener_height_det
  R = detection_dist
  Sigma = SigmaDet
  Em2 = Em2Det
  TrgHgt = Hs
  DetHgt = Hr
  If Iwthr1 = 0 Then
    Ingard    ' ground_effect
    If WindSpeed >= 0 Then ' propagating into wind_refrac_loss, keep calculated ground effect
      For I = 0 To 23
        prop_loss_indiv(I) = ground_effect(I) - ground_effect_reference(I)
        If prop_loss_indiv(I) = 0 Then prop_loss_indiv(I) = -0.001
        prop_loss_cum(I) = prop_loss_cum(I) + prop_loss_indiv(I)
      Next I
    Else     ' propagating with wind_refrac_loss, ground effect is in wind_refrac_loss attenuation already
      For I = 0 To 23
        If prop_loss_indiv(I) = 0 Then prop_loss_indiv(I) = -0.001
        prop_loss_indiv(I) = ground_effect(I) - ground_effect_reference(I)
        prop_loss_cum(I) = prop_loss_cum(I) + prop_loss_indiv(I)
      Next I
    End If
  Else
    If Iwthr1 < 13 Then
      'Upward(R,Attn)
      'For I = 0 to 23 ground_effect(I)=Attn^(I)
      'For I=0 to 23 do prop_loss_cum(I)=prop_loss_cum(I)+ground_effect(I)-ground_effect_reference(I)
      'For I=0 to 23 do prop_loss_cum(I)=prop_loss_cum(I)+ground_effect(I)-ground_effect_reference(I)*1
    Else ' Iwthr1>13
      'Downwd(R,Attn)
      'For I = 0 to 23 ground_effect(I)=Attn^(I)
      'For I=0 to 23 do prop_loss_cum(I)=prop_loss_cum(I)+ground_effect(I)-ground_effect_reference(I)
      'For I=0 to 23 do prop_loss_cum(I)=prop_loss_cum(I)+ground_effect(I)-ground_effect_reference(I)*1
    End If
  End If
End Sub
Function Barrier()
  For I = 0 To 23
    barrier_attenuation(I) = -0.001
  Next I
  If barrier_number > 0 Then
    If detection_dist >= distance_from_source Then
      Hba = barrier_height_det - (source_height_det + distance_from_source * (listener_height_det - source_height_det) / detection_dist)
      X1 = Sqr(distance_from_source * distance_from_source + (barrier_height_det - source_height_det) * (barrier_height_det - source_height_det))
      X1 = X1 + Sqr((detection_dist - distance_from_source) * (detection_dist - distance_from_source) + (barrier_height_det - listener_height_det) * (barrier_height_det - listener_height_det))
      X1 = X1 - Sqr(detection_dist * detection_dist + (listener_height_det - source_height_det) * (listener_height_det - source_height_det))
      Cs = 331.4 * Sqr(1 + celsius_degrees_det / 273.15)
      For I = 0 To 23
        If Hba < 0 Then
          X2 = 2 * X1 * freq(2) / Cs ' Fresnel number N
          X3 = Sqr(2 * Pi * Abs(X2))
          X2 = -X2
          If X2 <= -0.1916 Then
            barrier_attenuation(I) = -0.01
          Else
            barrier_attenuation(I) = -(5 + 2 * TenDivLog10 * Log(X3 / Tan(X3))) - 0.01
          End If
        Else ' Hba>=0
          X2 = 2 * X1 * freq(I) / Cs ' Fresnel number N
          X3 = Sqr(2 * Pi * Abs(X2))
          X2 = Exp(2 * X3)
          If X2 = 1 Then
            barrier_attenuation(I) = -5 - 0.01
          Else
            barrier_attenuation(I) = -(5 + 2 * TenDivLog10 * Log(X3 * (X2 + 1) / (X2 - 1))) - 0.01
            If barrier_attenuation(I) < -20 Then barrier_attenuation(I) = -20
          End If
        End If
      Next I
    End If
  End If
  For I = 0 To 23
    prop_loss_indiv(I) = barrier_attenuation(I)
    prop_loss_cum(I) = prop_loss_cum(I) + prop_loss_indiv(I)
  Next I
End Function
Function Foliage()
  For I = 0 To 23
    foliage_attenuation(I) = -0.001
  Next I
  If foliated_zone_nums > 0 Then
    If detection_dist > distance_frmSource_to_fol Then
      X2 = detection_dist - distance_frmSource_to_fol
      If X2 > extent_of_foliage Then X2 = extent_of_foliage
      X2 = Sqr(X2)
      Cons = 2.647 / Log(10)
      For I = 10 To 23
        Ka = (2 * Pi * freq(I) / Cs) * avg_leaf_width / 100
        If Ka < 0.401 Then
          foliage_attenuation(I) = -0.01
        ElseIf Ka < 5 Then
          foliage_attenuation(I) = -X2 * Sqr(leaf_area_per_unit) * (Cons * Log(Ka) + 1.05)
        Else
          foliage_attenuation(I) = -X2 * Sqr(leaf_area_per_unit) * 2.9
        End If
      Next I
    End If
  End If
  For I = 0 To 23
    prop_loss_indiv(I) = foliage_attenuation(I)
    prop_loss_cum(I) = prop_loss_cum(I) + prop_loss_indiv(I)
  Next I
End Function

Function NormalDeviate(P)
  T = Sqr(-2 * Log(P))
  NormalDeviate = T - (2.30753 + T * 0.27061) / (1 + T * (0.99229 + T * 0.04481))
End Function
Function Dprime()
  If hit_prob > 0.5 Then
    Z3 = NormalDeviate(1 - hit_prob)
  Else
    Z3 = NormalDeviate(hit_prob)
  End If
  If false_alarm_rate > 0.5 Then
    Z = NormalDeviate(1 - false_alarm_rate)
  Else
    Z = NormalDeviate(false_alarm_rate)
  End If
  calculate_d = Z3 - Z
  If hit_prob > 0.5 Then
    If false_alarm_rate > 0.5 Then
      calculate_d = Z3 - Z
    Else
      calculate_d = Z3 + Z
    End If
  Else
    calculate_d = Z - Z3
  End If
  Dprime = calculate_d
  For I = 0 To 23
    If I > 10 Then
      mod_background_noise(I) = background_noise_spec(I) + TenDivLog10 * Log(calculate_d / observer_efficiency / B3(I))
    Else
      E4 = 0
      For J = 0 To 4
        I0 = I + J - 2
        W4 = one_third_oct_band_weight(I, J)
        If (W4 > 0) And (I0 >= 0) Then E4 = E4 + W4 * Exp(Log10Div10 * background_noise_spec(I0))
      Next J
      mod_background_noise(I) = TenDivLog10 * Log(E4 * calculate_d / observer_efficiency / B3(I))
    End If
    If hearing_threshold(I) > mod_background_noise(I) Then
      max_modBackground_thresh(I) = hearing_threshold(I)
    Else
      max_modBackground_thresh(I) = mod_background_noise(I)
    End If
    mod_background_noise_2(I) = max_modBackground_thresh(I)
  Next I
End Function
Function ComplexDiv(A, B, C, D, E, F) As Integer
  G = 1 / (C * C + D * D)
  E = (A * C + B * D) * G
  F = (B * C - A * D) * G
  ComplexDiv = 1
End Function
Function ComplexMul(A, B, C, D, E, F) As Integer
  E = A * C - B * D
  F = A * D + B * C
  ComplexMul = 1
End Function
Function Polar(X, Y, R, Th)
  R = Sqr(X * X + Y * Y)
  Th = Application.Atan2(X, Y)
  Polar = 1
End Function
Sub SignalNoise()
  For I = 0 To 23
    If I > 10 Then
      detection_band(I) = target_spectrum(I) + prop_loss_cum(I) - max_modBackground_thresh(I)
    Else
      E3 = 0
      For J = 0 To 4
        I0 = I + J - 2
        W4 = one_third_oct_band_weight(I, J)
        If (W4 > 0) And (I0 >= 0) Then
          E3 = E3 + W4 * Exp(Log10Div10 * (target_spectrum(I0) + prop_loss_cum(I0)))
        End If
      Next J
      detection_band(I) = TenDivLog10 * Log(E3) - max_modBackground_thresh(I)
    End If
    listener_spectrum(I) = detection_band(I) + max_modBackground_thresh(I)
  Next I
  M2 = detection_band(0)
  B1 = 0
  For I = 1 To 23
    If detection_band(I) >= M2 Then
      M2 = detection_band(I)
      B1 = I
    End If
  Next I
End Sub
Sub Propagate()
  InverseDistance
  GroundEffect
  Barrier
  Foliage
  'If Iwthr1 = 0 Then Winds
  Atmosphere
  SignalNoise
End Sub
Sub BinarySearch()
  Z9 = -1
  detection_dist = m_meas_distance * 25
  Do
    Z9 = Z9 + 1
    D5 = detection_dist
    detection_dist = 2 * detection_dist
    D6 = detection_dist
    Propagate
  Loop Until M2 < 0
  If Z9 = 0 Then
    Do
      D6 = detection_dist
      detection_dist = detection_dist / 2
      D5 = detection_dist
      Propagate
    Loop Until M2 > 0
  End If
  Do
    detection_dist = (D5 + D6) / 2
    Propagate
    If M2 > 0 Then
      D5 = detection_dist
    Else
      D6 = detection_dist
    End If
  Loop Until Abs(D6 - D5) < precision_fraction * detection_dist
  'Target detectable at detection_dist:3:2 meters +/- (precision_fraction*detection_dist):3:2 meters in freq[B1] Hz Band
  'at (target_spectrum(B1)+prop_loss_cum(B1)):2:1 dB.
  If hearing_threshold(B1) < mod_background_noise(B1) Then
  'Detection limited by background noise.
  Else
  'Detection limited by human threshold of hearing.
  End If
End Sub

Sub BinarySearchA()
  Z9 = -1
  detection_dist = m_meas_distance * 25
  Do
    Z9 = Z9 + 1
    D5 = detection_dist
    detection_dist = 2 * detection_dist
    D6 = detection_dist
    Propagate
    Ea = 0: For I = 0 To 23: Ea = Ea + Exp(Log10Div10 * (target_spectrum(I) + prop_loss_cum(I) + A_weight_levels(I))): Next I
    prop_loss_cum(24) = TenDivLog10 * Log(Ea)
  Loop Until prop_loss_cum(24) < dBa
  If Z9 = 0 Then
    Do
      D6 = detection_dist
      detection_dist = detection_dist / 2
      D5 = detection_dist
      Propagate
      Ea = 0: For I = 0 To 23: Ea = Ea + Exp(Log10Div10 * (target_spectrum(I) + prop_loss_cum(I) + A_weight_levels(I))): Next I
      prop_loss_cum(24) = TenDivLog10 * Log(Ea)
    Loop Until prop_loss_cum(24) > dBa
  End If
  Do
    detection_dist = (D5 + D6) / 2
    Propagate
    Ea = 0: For I = 0 To 23: Ea = Ea + Exp(Log10Div10 * (target_spectrum(I) + prop_loss_cum(I) + A_weight_levels(I))): Next I
    prop_loss_cum(24) = TenDivLog10 * Log(Ea)
    If prop_loss_cum(24) > dBa Then
      D5 = detection_dist
    Else
      D6 = detection_dist
    End If
  Loop Until Abs(D6 - D5) < precision_fraction * detection_dist
End Sub

Sub Reverse()
  Propagate
End Sub
Sub ReferenceCalc()
  Hs = source_height_meas_cond
  Hr = mic_height_meas_cond
  R = m_meas_distance
  Sigma = SigmaRef
  Em2 = Em2Ref
  Ingard  'ground_effect
  AnsiHumidity  ' atmos_absorption
  For I = 0 To 23
    ground_effect_reference(I) = ground_effect(I)
    AtmAbsRef(I) = atmos_absorption(I)
  Next I
End Sub
Function TargetdBA()
  Dim E As Double
  E = 0
  For I = 0 To 23
    E = E + Exp(Log10Div10 * (target_spectrum(I) + A_weight_levels(I)))
  Next I
  TargetdBA = TenDivLog10 * Log(E)
End Function
Function ListenerdBA()
  Dim E As Double
  E = 0
  For I = 0 To 23
    E = E + Exp(Log10Div10 * (target_spectrum(I) + prop_loss_cum(I) + A_weight_levels(I)))
  Next I
  ListenerdBA = TenDivLog10 * Log(E)
End Function
Function BkgdBA()
  Dim E As Double
  E = 0
  For I = 0 To 23
    E = E + Exp(Log10Div10 * (background_noise_spec(I) + A_weight_levels(I)))
  Next I
  BkgdBA = TenDivLog10 * Log(E)
End Function
Sub InitMacros()
  Dim Rg1 As Range
  DeleteAllCharts
  Trg = Sheets("Model").Range("B6").Value
  Bkg = Sheets("Model").Range("C6").Value
  Hth = Sheets("Model").Range("D6").Value
  For I = 0 To 23
    freq(I) = Sheets("Data").Range("A2").Offset(I, 0).Value
    Sheets("Model").Range("A8").Offset(I, 0).Value = freq(I)
  Next I
  If Trg <> 0 Then
    For I = 0 To 23
      target_spectrum(I) = Sheets("Data").Range("A2").Offset(I, Trg).Value
      Sheets("Model").Range("A8").Offset(I, 1).Value = target_spectrum(I)
    Next I
    m_meas_distance = Sheets("Data").Range("A2").Offset(24, Trg).Value
    Sheets("Model").Range("C3").Value = m_meas_distance
  Else
    Set Rg1 = Sheets("Model").Range("C3")
    If IsNumeric(Rg1) And Len(Rg1) > 0 Then
      m_meas_distance = Rg1.Value
    Else
      MsgBox ("Enter a Distance in meters in Cell C3 for Microphone from Source")
      End
    End If
  End If
  If Bkg <> 0 Then
    For I = 0 To 23
      background_noise_spec(I) = Sheets("Data").Range("A2").Offset(I, Bkg).Value
      Sheets("Model").Range("A8").Offset(I, 2).Value = background_noise_spec(I)
    Next I
  End If
  If Hth <> 0 Then
    For I = 0 To 23
      hearing_threshold(I) = Sheets("Data").Range("A2").Offset(I, Hth).Value
      Sheets("Model").Range("A8").Offset(I, 3).Value = hearing_threshold(I)
    Next I
  End If
  For I = 0 To 23: A_weight_levels(I) = Sheets("Data").Range("A2").Offset(I, 35).Value:  Next I  ' 35 A weight levels
  For I = 0 To 23: AiWt(I) = Sheets("Data").Range("A2").Offset(I, 37).Value: Next I ' 37 AI weight levels
  ColOffset = 29
  For I = 0 To 23
    For J = 0 To 4
      one_third_oct_band_weight(I, J) = Sheets("Data").Range("A2").Offset(I, J + ColOffset).Value
    Next J
    J = 5
    B3(I) = Sheets("Data").Range("A2").Offset(I, J + ColOffset).Value
  Next I
  Sheets("Model").Range("A39").Select
  Pi = 4 * Atn(1)
  Log10Div10 = 0.230258509
  TenDivLog10 = 1 / Log10Div10
  Log10 = 2.302585093
  Co = 331.32
  zerDegCelsius_inKelvin = 273.15
  precision_fraction = 0.001 ' precision fraction
  source_height_meas_cond = Range("A3").Value ' m Source Height meas
  mic_height_meas_cond = Range("B3").Value ' m Mic Height meas
  ' m_meas_distance m in "C3"  Mic distance from source
  Tref = Range("m_meas_distance").Value     ' Deg C meas
  Href = Range("E3").Value     ' % r.h. meas
  SigmaRef = Range("F3").Value ' m Flow resistivity meas
  Em2Ref = Range("G3").Value   ' Em2 turbulence factor meas
  celsius_degrees_det = Tref
  relative_humid_percent = Href
  
  ReferenceCalc
  
  source_height_det = Range("H3").Value ' m Source Height (detection conditions)
  listener_height_det = Range("I3").Value ' m Listener Height det
  celsius_degrees_det = Range("J3").Value    ' Deg C det
  relative_humid_percent = Range("K3").Value    ' % relative humidity det
  SigmaDet = Range("L3").Value ' m Flow resistivity det
  Em2Det = Range("M3").Value   ' Em2 turbulence factor det
  WindSpeed = Range("N3").Value ' wind_refrac_loss speed det
  observer_efficiency = Range("r3").Value  ' Observer efficiency
  hit_prob = Range("s3").Value  ' Hit prob
  false_alarm_rate = Range("t3").Value  ' False alarm prop
  calculate_d = Dprime             ' Calculate d' statistic
  Range("u3").Value = calculate_d
  WindFlag = 0
  WindDir = " Upwind"
  barrier_number = Range("H5").Value ' barrier? 0 or 1
  distance_from_source = Range("I5").Value ' distance from source m
  barrier_height_det = Range("J5").Value ' barrier height detec m
  foliated_zone_nums = Range("K5").Value ' foliage? 0 or 1
  distance_frmSource_to_fol = Range("L5").Value ' distance in meters from source to near edge of foliage
  extent_of_foliage = Range("M5").Value ' depth (extent) of foliage in meters
  leaf_area_per_unit = Range("N5").Value ' leaf area per unit vol dense hardwood brush in m^-1
  avg_leaf_width = Range("O5").Value ' average leaf width in cm.
  type_of_surface = 1   ' Type of surface
  Iwthr1 = 0
  Surface = "Grass"
  deflt_bkg_num = 11 ' Default background number 12=G.C. 11=Low EPA
  typical_vehicle = 3   ' Typical vehicle
  Hnumber = 2  ' ISO Hearing Threshold for Pure tones

End Sub

Sub Detection_Distance_Macro()
  InitMacros

  BinarySearch

  Range("A32:D37").Select
  Range("D37").Activate
  Selection.HorizontalAlignment = xlLeft
  Range("A32").Select
  ActiveCell.Offset(0, 0).Value = "Target Level (dBA)"
  ActiveCell.Offset(1, 0).Value = "Detect Distance (m)"
  ActiveCell.Offset(2, 0).Value = "Detect Frequency (Hz)"
  ActiveCell.Offset(3, 0).Value = "Detect 1/3 Oct Band Level"
  ActiveCell.Offset(4, 0).Value = "Detect Level (dBA)"
  ActiveCell.Offset(5, 0).Value = "Background Noise Level (dBA)"
  ActiveCell.Offset(0, 4).Value = TargetdBA
  ActiveCell.Offset(1, 4).Value = detection_dist
  ActiveCell.Offset(2, 4).Value = freq(B1)
  ActiveCell.Offset(3, 4).Value = target_spectrum(B1) + prop_loss_cum(B1)
  ActiveCell.Offset(4, 4).Value = ListenerdBA
  ActiveCell.Offset(5, 4).Value = BkgdBA
  Range("G32:J37").Select
  Range("J37").Activate
  Selection.HorizontalAlignment = xlLeft
  Range("G32").Select
  ActiveCell.Offset(0, 0).Value = ""
  ActiveCell.Offset(1, 0).Value = ""
  ActiveCell.Offset(2, 0).Value = ""
  ActiveCell.Offset(3, 0).Value = ""
  ActiveCell.Offset(4, 0).Value = ""
  ActiveCell.Offset(5, 0).Value = ""
  ActiveCell.Offset(0, 4).Value = ""
  ActiveCell.Offset(1, 4).Value = ""
  ActiveCell.Offset(2, 4).Value = ""
  ActiveCell.Offset(3, 4).Value = ""
  ActiveCell.Offset(4, 4).Value = ""
  ActiveCell.Offset(5, 4).Value = ""
  Range("A8").Select

  InverseDistance
  K = 4: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  Atmosphere
  K = 5: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  GroundEffect
  K = 7: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  Barrier
  K = 8: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  Foliage
  K = 9:  For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  K = 10: For I = 0 To 23: ActiveCell.Offset(I, K).Value = target_spectrum(I) + prop_loss_cum(I): Next I
  K = 11: For I = 0 To 23: ActiveCell.Offset(I, K).Value = listener_spectrum(I): Next I
  K = 12: For I = 0 To 23: ActiveCell.Offset(I, K).Value = mod_background_noise(I): Next I
  K = 13: For I = 0 To 23: ActiveCell.Offset(I, K).Value = max_modBackground_thresh(I): Next I
  Range("B8:B31").Select:  Selection.Interior.ColorIndex = 34
  Range("c8:c31").Select:  Selection.Interior.ColorIndex = 38
  Range("k8:k31").Select:  Selection.Interior.ColorIndex = 35
  DrawLayers
  Range("B8:K37").Select
  Range("K37").Activate
  Selection.NumberFormat = "0.0"
  With Selection
    .HorizontalAlignment = xlCenter
  End With
  Range("A38").Activate
End Sub
Sub NonDetect_Spectrum_Macro()
  Dim Rg2 As Range
  InitMacros
  Set Rg2 = Sheets("Model").Range("C4")
  If IsNumeric(Rg2) And Len(Rg2) > 0 Then
    detection_dist = Rg2.Value
  Else
    MsgBox ("Enter a Distance in meters in Cell C4 for Listener from Source")
    End
  End If
  Reverse

  Range("A32:D37").Select
  Range("D37").Activate
  Selection.HorizontalAlignment = xlLeft
  Range("A32").Select
  ActiveCell.Offset(0, 0).Value = "Target Level (dBA)"
  ActiveCell.Offset(1, 0).Value = "Detect Distance (m)"
  ActiveCell.Offset(2, 0).Value = "Detect Frequency (Hz)"
  ActiveCell.Offset(3, 0).Value = "Detect 1/3 Oct Band Level"
  ActiveCell.Offset(4, 0).Value = "Detect Level (dBA)"
  ActiveCell.Offset(5, 0).Value = "Background Noise Level (dBA)"
  ActiveCell.Offset(0, 4).Value = 0
  ActiveCell.Offset(1, 4).Value = 0
  ActiveCell.Offset(2, 4).Value = 0
  ActiveCell.Offset(3, 4).Value = 0
  ActiveCell.Offset(4, 4).Value = 0
  ActiveCell.Offset(5, 4).Value = 0
  Range("G32:J37").Select
  Range("J37").Activate
  Selection.HorizontalAlignment = xlLeft
  Range("G32").Select
  ActiveCell.Offset(0, 0).Value = ""
  ActiveCell.Offset(1, 0).Value = ""
  ActiveCell.Offset(2, 0).Value = ""
  ActiveCell.Offset(3, 0).Value = ""
  ActiveCell.Offset(4, 0).Value = ""
  ActiveCell.Offset(5, 0).Value = ""
  ActiveCell.Offset(0, 4).Value = ""
  ActiveCell.Offset(1, 4).Value = ""
  ActiveCell.Offset(2, 4).Value = ""
  ActiveCell.Offset(3, 4).Value = ""
  ActiveCell.Offset(4, 4).Value = ""
  ActiveCell.Offset(5, 4).Value = ""
  Range("A8").Select
  K = 1: For I = 0 To 23
    ActiveCell.Offset(I, K).Value = mod_background_noise_2(I) - prop_loss_cum(I)
    target_spectrum(I) = mod_background_noise_2(I) - prop_loss_cum(I)
  Next I
  InverseDistance
  K = 4: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  Atmosphere
  K = 5: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  GroundEffect
  K = 7: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  Barrier
  K = 8: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  Foliage
  K = 9:  For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  K = 10: For I = 0 To 23: ActiveCell.Offset(I, K).Value = target_spectrum(I) + prop_loss_cum(I): Next I
  K = 11: For I = 0 To 23: ActiveCell.Offset(I, K).Value = listener_spectrum(I): Next I
  K = 12: For I = 0 To 23: ActiveCell.Offset(I, K).Value = mod_background_noise(I): Next I
  K = 13: For I = 0 To 23: ActiveCell.Offset(I, K).Value = max_modBackground_thresh(I): Next I
  
  Range("B8:B31").Select:  Selection.Interior.ColorIndex = 34
  Range("c8:c31").Select:  Selection.Interior.ColorIndex = 38
  Range("k8:k31").Select:  Selection.Interior.ColorIndex = 35
  DeleteAllCharts
'  DrawLayers
  Range("B8:K37").Select
  Range("K37").Activate
  Selection.NumberFormat = "0.0"
  With Selection
    .HorizontalAlignment = xlCenter
  End With
  Range("A38").Activate
End Sub
Sub Level_Losses_Distance_Macro()
  Dim Rg2 As Range
  InitMacros
  Set Rg2 = Sheets("Model").Range("C4")
  If IsNumeric(Rg2) And Len(Rg2) > 0 Then
    detection_dist = Rg2.Value
  Else
    MsgBox ("Enter a Distance in meters in Cell C4 for Listener from Source")
    End
  End If
  Propagate
  Range("A32:D37").Select
  Range("D37").Activate
  Selection.HorizontalAlignment = xlLeft
  Range("A32").Select
  ActiveCell.Offset(0, 0).Value = "Target Level (dBA)"
  ActiveCell.Offset(1, 0).Value = "Detect Distance (m)"
  ActiveCell.Offset(2, 0).Value = "Detect Frequency (Hz)"
  ActiveCell.Offset(3, 0).Value = "Detect 1/3 Oct Band Level"
  ActiveCell.Offset(4, 0).Value = "Detect Level (dBA)"
  ActiveCell.Offset(5, 0).Value = "Background Noise Level (dBA)"
  ActiveCell.Offset(0, 4).Value = TargetdBA
  ActiveCell.Offset(1, 4).Value = detection_dist
  ActiveCell.Offset(2, 4).Value = freq(B1)
  ActiveCell.Offset(3, 4).Value = target_spectrum(B1) + prop_loss_cum(B1)
  ActiveCell.Offset(4, 4).Value = ListenerdBA
  ActiveCell.Offset(5, 4).Value = BkgdBA
  Range("G32:J37").Select
  Range("J37").Activate
  Selection.HorizontalAlignment = xlLeft
  Range("G32").Select
  ActiveCell.Offset(0, 0).Value = ""
  ActiveCell.Offset(1, 0).Value = ""
  ActiveCell.Offset(2, 0).Value = ""
  ActiveCell.Offset(3, 0).Value = ""
  ActiveCell.Offset(4, 0).Value = ""
  ActiveCell.Offset(5, 0).Value = ""
  ActiveCell.Offset(0, 4).Value = ""
  ActiveCell.Offset(1, 4).Value = ""
  ActiveCell.Offset(2, 4).Value = ""
  ActiveCell.Offset(3, 4).Value = ""
  ActiveCell.Offset(4, 4).Value = ""
  ActiveCell.Offset(5, 4).Value = ""
  Range("A8").Select

  InverseDistance
  K = 4: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  Atmosphere
  K = 5: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  GroundEffect
  K = 7: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  Barrier
  K = 8: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  Foliage
  K = 9:  For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  K = 10: For I = 0 To 23: ActiveCell.Offset(I, K).Value = target_spectrum(I) + prop_loss_cum(I): Next I
  K = 11: For I = 0 To 23: ActiveCell.Offset(I, K).Value = listener_spectrum(I): Next I
  K = 12: For I = 0 To 23: ActiveCell.Offset(I, K).Value = mod_background_noise(I): Next I
  K = 13: For I = 0 To 23: ActiveCell.Offset(I, K).Value = max_modBackground_thresh(I): Next I
  Range("B8:B31").Select:  Selection.Interior.ColorIndex = 34
  Range("c8:c31").Select:  Selection.Interior.ColorIndex = 38
  Range("k8:k31").Select:  Selection.Interior.ColorIndex = 35
  DrawLayers
  
End Sub
Sub dBA_Distance_Macro()
  Dim Rg2 As Range
  InitMacros
  Set Rg2 = Sheets("Model").Range("G4")
  If IsNumeric(Rg2) And Len(Rg2) > 0 Then
    dBa = Rg2.Value
  Else
    MsgBox ("Enter a dBA value in Cell G4")
    End
  End If
  BinarySearchA
  Range("A32:D37").Select
  Range("D37").Activate
  Selection.HorizontalAlignment = xlLeft
  Range("A32").Select
  ActiveCell.Offset(0, 0).Value = "Target Level (dBA)"
  ActiveCell.Offset(1, 0).Value = "A-wt Distance (m)"
  ActiveCell.Offset(2, 0).Value = ""
  ActiveCell.Offset(3, 0).Value = ""
  ActiveCell.Offset(4, 0).Value = "Detect Level (dBA)"
  ActiveCell.Offset(5, 0).Value = "Background Noise Level (dBA)"
  ActiveCell.Offset(0, 4).Value = TargetdBA
  ActiveCell.Offset(1, 4).Value = detection_dist
  ActiveCell.Offset(2, 4).Value = ""
  ActiveCell.Offset(3, 4).Value = ""
  ActiveCell.Offset(4, 4).Value = ListenerdBA
  ActiveCell.Offset(5, 4).Value = BkgdBA
  Range("G32:J37").Select
  Range("J37").Activate
  Selection.HorizontalAlignment = xlLeft
  Range("G32").Select
  ActiveCell.Offset(0, 0).Value = ""
  ActiveCell.Offset(1, 0).Value = ""
  ActiveCell.Offset(2, 0).Value = ""
  ActiveCell.Offset(3, 0).Value = ""
  ActiveCell.Offset(4, 0).Value = ""
  ActiveCell.Offset(5, 0).Value = ""
  ActiveCell.Offset(0, 4).Value = ""
  ActiveCell.Offset(1, 4).Value = ""
  ActiveCell.Offset(2, 4).Value = ""
  ActiveCell.Offset(3, 4).Value = ""
  ActiveCell.Offset(4, 4).Value = ""
  ActiveCell.Offset(5, 4).Value = ""
  Range("A8").Select

  InverseDistance
  K = 4: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  Atmosphere
  K = 5: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  GroundEffect
  K = 7: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  Barrier
  K = 8: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  Foliage
  K = 9:  For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  K = 10: For I = 0 To 23: ActiveCell.Offset(I, K).Value = target_spectrum(I) + prop_loss_cum(I): Next I
  K = 11: For I = 0 To 23: ActiveCell.Offset(I, K).Value = listener_spectrum(I): Next I
  K = 12: For I = 0 To 23: ActiveCell.Offset(I, K).Value = mod_background_noise(I): Next I
  K = 13: For I = 0 To 23: ActiveCell.Offset(I, K).Value = max_modBackground_thresh(I): Next I
  Range("B8:B31").Select:  Selection.Interior.ColorIndex = 34
  Range("c8:c31").Select:  Selection.Interior.ColorIndex = 38
  Range("k8:k31").Select:  Selection.Interior.ColorIndex = 35
  DrawLayers
  Range("B8:K37").Select
  Range("K37").Activate
  Selection.NumberFormat = "0.0"
  With Selection
    .HorizontalAlignment = xlCenter
  End With
  Range("A38").Activate



End Sub





Sub Speech_Intelligibility_Macro()
  Dim Rg2 As Range
  InitMacros
  Set Rg2 = Sheets("Model").Range("C4")
  If IsNumeric(Rg2) And Len(Rg2) > 0 Then
    detection_dist = Rg2.Value
  Else
    MsgBox ("Enter a Distance in meters in Cell C4 for Listener from Source")
    End
  End If
  Propagate
  Ai = 0
  For I = 0 To 23 ' max_modBackground_thresh=max(hearing_threshold,background_noise_spec)
    If hearing_threshold(I) > background_noise_spec(I) Then
       max_modBackground_thresh(I) = hearing_threshold(I)
    Else
      max_modBackground_thresh(I) = background_noise_spec(I)
    End If
    signal_noise_ratio(I) = target_spectrum(I) + prop_loss_cum(I) - max_modBackground_thresh(I)
    If signal_noise_ratio(I) > 30 Then signal_noise_ratio(I) = 30
    If signal_noise_ratio(I) < 0 Then signal_noise_ratio(I) = 0
    Ai = Ai + signal_noise_ratio(I) * AiWt(I)
  Next I
  PCor = Ai * (108.9 + Ai * (370.2 + Ai * (-2527.1 + Ai * 4485.3)))
  PCor = PCor / (1 + Ai * (-5.839 + Ai * (21.89 + Ai * (-45.01 + Ai * 52.34))))
  Pc1 = Ai * (58.2 + Ai * (-9.21 + Ai * 328.8))
  Pc1 = Pc1 / (1 + Ai * (-2.596 + Ai * (4.27 + Ai * 1.145)))
  Pc2 = Ai * (442.22 + Ai * (-9090 + Ai * 189940))
  Pc2 = Pc2 / (1 + Ai * (-1.4901 + Ai * (-80.722 + Ai * 1894.1)))
  Pc3 = Ai * (-3.1869 + Ai * 266.79)
  Pc3 = Pc3 / (1 + Ai * (-2.2109 + Ai * 3.9026))
  Pc4 = Ai * (88.676 + Ai * (147.51 + Ai * 1335.1))
  Pc4 = Pc4 / (1 + Ai * (-3.0167 + Ai * (6.67 + Ai * 11.147)))
  
  Range("A32:D37").Select
  Range("D37").Activate
  Selection.HorizontalAlignment = xlLeft
  Range("A32").Select
  ActiveCell.Offset(0, 0).Value = "Target Level (dBA)"
  ActiveCell.Offset(1, 0).Value = "Detect Distance (m)"
  ActiveCell.Offset(2, 0).Value = "Detect Frequency (Hz)"
  ActiveCell.Offset(3, 0).Value = "Detect 1/3 Oct Band Level"
  ActiveCell.Offset(4, 0).Value = "Detect Level (dBA)"
  ActiveCell.Offset(5, 0).Value = "Background Noise Level (dBA)"
  ActiveCell.Offset(0, 4).Value = 0
  ActiveCell.Offset(1, 4).Value = 0
  ActiveCell.Offset(2, 4).Value = 0
  ActiveCell.Offset(3, 4).Value = 0
  ActiveCell.Offset(4, 4).Value = 0
  ActiveCell.Offset(5, 4).Value = 0
  Range("G32:J37").Select
  Range("J37").Activate
  Selection.HorizontalAlignment = xlLeft
  Range("G32").Select
  ActiveCell.Offset(0, 0).Value = "Articulation( AI )"
  ActiveCell.Offset(1, 0).Value = "Digits or Ciphers %"
  ActiveCell.Offset(2, 0).Value = "Modified Rhyme Test %"
  ActiveCell.Offset(3, 0).Value = "UnKnown Sentences %"
  ActiveCell.Offset(4, 0).Value = "PB50 1000 Words %"
  ActiveCell.Offset(5, 0).Value = "CVC Syllables %"
  
  ActiveCell.Offset(0, 4).Value = Ai
  ActiveCell.Offset(1, 4).Value = Pc2
  ActiveCell.Offset(2, 4).Value = PCor
  ActiveCell.Offset(3, 4).Value = Pc4
  ActiveCell.Offset(4, 4).Value = Pc1
  ActiveCell.Offset(5, 4).Value = Pc3

  Range("A8").Select
  K = 1: For I = 0 To 23: ActiveCell.Offset(I, K).Value = mod_background_noise_2(I) - prop_loss_cum(I): Next I
  InverseDistance
  K = 4: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  Atmosphere
  K = 5: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  GroundEffect
  K = 7: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  Barrier
  K = 8: For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  Foliage
  K = 9:  For I = 0 To 23: ActiveCell.Offset(I, K).Value = prop_loss_indiv(I): Next I
  K = 10: For I = 0 To 23: ActiveCell.Offset(I, K).Value = target_spectrum(I) + prop_loss_cum(I): Next I
  Range("B8:B31").Select:  Selection.Interior.ColorIndex = 34
  Range("c8:c31").Select:  Selection.Interior.ColorIndex = 38
  Range("k8:k31").Select:  Selection.Interior.ColorIndex = 35


  Range("K32:K37").Select
  Range("K37").Activate
  Selection.NumberFormat = "0.00"
  With Selection
    .HorizontalAlignment = xlCenter
  End With

  Sheets("Chart 2").Select
  For Each cht In ActiveSheet.ChartObjects
    cht.Copy
  Next cht
  Sheets("Model").Select
  Range("L6").Offset(0).Select
  ActiveSheet.Paste
  ShapeNo = 7: For I = 0 To 23: prop_loss_cum_4graph(I) = max_modBackground_thresh(I): Next I: TraceXY  ' thick solid red
  ShapeNo = 10: For I = 0 To 23: prop_loss_cum_4graph(I) = target_spectrum(I) + prop_loss_cum(I): Next I: TraceXY ' thick solid green
  ShapeNo = 1: For I = 0 To 23: prop_loss_cum_4graph(I) = target_spectrum(I): Next I:  TraceXY   ' thick solid cyan
End Sub
Sub AddChartSheet()
  Dim chtChart As Chart
  'Create a new chart
  Set chtChart = Charts.Add
  With chtChart
    .Name = "Chart1"
    .ChartType = xlXYScatterLinesNoMarkers
    'Link to the source data range
    .SetSourceData source:=Sheets("Model").Range("A38:B61"), PlotBy:=xlColumns
    .HasTitle = True
    .ChartTitle.Text = "=Model!R37C1"
    .Axes(xlCategory, xlPrimary).MinimumScale = 0
    .Axes(xlCategory, xlPrimary).MaximumScale = 25
    .Axes(xlCategory, xlPrimary).HasTitle = True
    .Axes(xlCategory, xlPrimary).AxisTitle.Characters.Text = "Frequency"
    .Axes(xlValue, xlPrimary).MinimumScale = 0
    .Axes(xlValue, xlPrimary).MaximumScale = 90
    .Axes(xlValue, xlPrimary).TickLabels.NumberFormat = "0"
    .Axes(xlValue, xlPrimary).HasTitle = True
    .Axes(xlValue, xlPrimary).AxisTitle.Characters.Text = "Sound Pressure Level"
  End With
End Sub

Sub DeleteAllCharts()
  Dim chtObj As ChartObject
  For Each chtObj In ActiveSheet.ChartObjects
    chtObj.Delete
  Next chtObj
End Sub

Sub CopyCharts()
  Dim chtObj As ChartObject
  'First delete all charts from target sheet
  DeleteEmbeddedCharts ("Model")
  'Some delay
  'Application.Wait Now + TimeSerial(0, 0, 1)
  
  For Each chtObject In Sheets("Chart 2").ChartObjects
    With Sheets("Model")
      .Activate
      chtObj.Copy
      'Paste in row L7+i
      I = 0
      Range("L7").Offset(I).Select
      .Activate
      'Application.Wait Now + TimeSerial(0, 0, 1)
      .Paste
      'Application.Wait Now + TimeSerial(0, 0, 1)
      .Activate
    End With
  Next chtObj
  ' Set the data references to target sheet
  SetChartRef ("Model")
End Sub
Sub TraceXY()
  Dim myCht As Chart
  Dim mySrs As Series
  Dim Npts As Integer, Ipts As Integer
  Dim myBuilder As FreeformBuilder
  Dim myShape As Shape
  Dim Xnode As Double, Ynode As Double
  Dim Xmin As Double, Xmax As Double
  Dim Ymin As Double, Ymax As Double
  Dim Xleft As Double, Ytop As Double
  Dim Xwidth As Double, Yheight As Double
  Dim Y As Double
  
  Set myCht = ActiveChart
  Xleft = myCht.PlotArea.InsideLeft
  Xwidth = myCht.PlotArea.InsideWidth
  Ytop = myCht.PlotArea.InsideTop
  Yheight = myCht.PlotArea.InsideHeight
  Xmin = myCht.Axes(1).MinimumScale
  Xmax = myCht.Axes(1).MaximumScale
  Ymin = myCht.Axes(2).MinimumScale
  Ymax = myCht.Axes(2).MaximumScale
  
  If ShapeNo = 2 Then Npts = 11 Else Npts = 23
  'First point
  Y = prop_loss_cum_4graph(0): If Y < Ymin Then Y = Ymin: If Y > Ymax Then Y = Ymax
  Xnode = Xleft + (1 - Xmin) * Xwidth / (Xmax - Xmin)
  Ynode = Ytop + (Ymax - Y) * Yheight / (Ymax - Ymin)
  Set myBuilder = myCht.Shapes.BuildFreeform(msoEditingAuto, Xnode, Ynode)
  For Ipts = 1 To Npts
    Y = prop_loss_cum_4graph(Ipts): If Y < Ymin Then Y = Ymin: If Y > Ymax Then Y = Ymax
    Xnode = Xleft + (Ipts + 1 - Xmin) * Xwidth / (Xmax - Xmin)
    Ynode = Ytop + (Ymax - Y) * Yheight / (Ymax - Ymin)
    myBuilder.AddNodes msoSegmentLine, msoEditingAuto, Xnode, Ynode
  Next Ipts
  Set myShape = myBuilder.ConvertToShape
  With myShape
    Select Case ShapeNo
    Case 1
      .Line.ForeColor.RGB = RGB(0, 255, 255) ' Cyan
      .Line.DashStyle = msoLineSolid
      .Line.Weight = 2
    Case 2
      .Line.ForeColor.RGB = RGB(0, 255, 0) ' Green
      .Line.DashStyle = msoLineSolid
      .Line.Weight = 2
    Case 3
      .Line.ForeColor.RGB = RGB(0, 255, 0) ' Green
      .Line.DashStyle = msoLineDash
      .Line.Weight = 3
    Case 4
      .Line.ForeColor.RGB = RGB(255, 0, 0) ' Red
      .Line.DashStyle = msoLineDash
      .Line.Weight = 3
    Case 5
      .Line.ForeColor.RGB = RGB(255, 0, 0) ' Red
      .Line.DashStyle = msoLineSolid
      .Line.Weight = 0.75
    Case 6
      .Line.ForeColor.RGB = RGB(255, 0, 0) 'Red
      .Line.DashStyle = msoLineLongDash
      .Line.Weight = 0.75
    Case 7
      .Line.ForeColor.RGB = RGB(255, 0, 0) 'Red
      .Line.DashStyle = msoLineSolid
      .Line.Weight = 2
    Case 8
      .Line.ForeColor.RGB = RGB(192, 192, 192) 'Silver
      .Line.DashStyle = msoLineSolid
      .Line.Weight = 2
    Case 9
      .Line.ForeColor.SchemeColor = 16 ' Brown
      .Line.DashStyle = msoLineSolid
      .Line.Weight = 0.75
    Case 10
      .Line.ForeColor.RGB = RGB(0, 255, 0) ' Green
      .Line.DashStyle = msoLineSolid
      .Line.Weight = 3
    Case Else
      .Line.ForeColor.SchemeColor = 16 ' Brown
      .Line.DashStyle = msoLineSolid
      .Line.Weight = 0.75
    End Select
  End With
End Sub

Sub DrawFilledPolygon()
  Dim myCht As Chart
  Dim mySrs As Series
  Dim Npts As Integer, Ipts As Integer
  Dim myBuilder As FreeformBuilder
  Dim myShape As Shape
  Dim Xnode As Double, Ynode As Double
  Dim Xmin As Double, Xmax As Double
  Dim Ymin As Double, Ymax As Double
  Dim Xleft As Double, Ytop As Double
  Dim Xwidth As Double, Yheight As Double
  Dim Y As Double
  
  Set myCht = ActiveChart
  Xleft = myCht.PlotArea.InsideLeft
  Xwidth = myCht.PlotArea.InsideWidth
  Ytop = myCht.PlotArea.InsideTop
  Yheight = myCht.PlotArea.InsideHeight
  Xmin = myCht.Axes(1).MinimumScale
  Xmax = myCht.Axes(1).MaximumScale
  Ymin = myCht.Axes(2).MinimumScale
  Ymax = myCht.Axes(2).MaximumScale

  Y = prop_loss_cum_4graph(0) + prop_loss_indiv(0): If Y < Ymin Then Y = Ymin: If Y > Ymax Then Y = Ymax
  'First point
  Xnode = Xleft + (1 - Xmin) * Xwidth / (Xmax - Xmin)
  Ynode = Ytop + (Ymax - Y) * Yheight / (Ymax - Ymin)
  Set myBuilder = myCht.Shapes.BuildFreeform(msoEditingAuto, Xnode, Ynode)
  For Ipts = 0 To 23
    Y = prop_loss_cum_4graph(Ipts): If Y < Ymin Then Y = Ymin: If Y > Ymax Then Y = Ymax
    Xnode = Xleft + (Ipts + 1 - Xmin) * Xwidth / (Xmax - Xmin)
    Ynode = Ytop + (Ymax - Y) * Yheight / (Ymax - Ymin)
    myBuilder.AddNodes msoSegmentLine, msoEditingAuto, Xnode, Ynode
  Next Ipts
  For Ipts = 23 To 0 Step -1
    Y = prop_loss_cum_4graph(Ipts) + prop_loss_indiv(Ipts): If Y < Ymin Then Y = Ymin: If Y > Ymax Then Y = Ymax
    Xnode = Xleft + (Ipts + 1 - Xmin) * Xwidth / (Xmax - Xmin)
    Ynode = Ytop + (Ymax - Y) * Yheight / (Ymax - Ymin)
    myBuilder.AddNodes msoSegmentLine, msoEditingAuto, Xnode, Ynode
  Next Ipts
  For Ipts = 0 To 23
    prop_loss_cum_4graph(Ipts) = prop_loss_cum_4graph(Ipts) + prop_loss_indiv(Ipts)
  Next Ipts
  
  Set myShape = myBuilder.ConvertToShape
  With myShape
    .Line.ForeColor.SchemeColor = 16  ' Brown
    .Line.Weight = 0.75
    Select Case ShapeNo
    Case 1
      .Fill.Patterned msoPattern10Percent
      .Fill.ForeColor.RGB = RGB(0, 255, 255) ' Cyan
      .Fill.BackColor.RGB = RGB(160, 160, 160) ' Gray
  
    Case 2
      .Fill.Patterned msoPatternSmallCheckerBoard
      .Fill.ForeColor.SchemeColor = 4 ' Blue
      '.Fill.ForeColor.RGB = RGB(0, 0, 255) ' Blue
    Case 3
      .Fill.Patterned msoPatternNarrowHorizontal
      .Fill.ForeColor.RGB = RGB(0, 127, 0) ' DarkGreen
    Case 4
      .Fill.Patterned msoPatternDarkVertical
      .Fill.ForeColor.RGB = RGB(0, 127, 0) ' DarkGreen
    Case 5
      .Fill.Patterned msoPatternHorizontalBrick
      .Fill.ForeColor.RGB = RGB(190, 190, 190) ' Gray
      .Fill.BackColor.RGB = RGB(127, 0, 0) ' Brown
    Case 6
      .Fill.Patterned msoPattern25Percent
      .Fill.ForeColor.SchemeColor = 7 ' LightCyan
    Case Else
      .Fill.Patterned msoPattern25Percent
      .Fill.ForeColor.SchemeColor = 7 ' LightCyan
    End Select
  End With
End Sub
Sub DrawLayers()
  Dim cht As ChartObject
  
  DeleteAllCharts
  Sheets("Chart 2").Select
  For Each cht In ActiveSheet.ChartObjects
    cht.Copy
  Next cht
  
  Sheets("Model").Select
  Range("L6").Offset(0).Select
  ActiveSheet.Paste

  For I = 0 To 23: prop_loss_cum_4graph(I) = target_spectrum(I): prop_loss_cum(I) = target_spectrum(I):  Next I
  InverseDistance ' prop_loss_cum=prop_loss_cum+Pinv
  ShapeNo = 1
  DrawFilledPolygon
  Atmosphere ' prop_loss_cum=prop_loss_cum+atmos_absorption-AtmAbsRef
  ShapeNo = 2
  DrawFilledPolygon
  GroundEffect ' prop_loss_cum=prop_loss_cum+GeAbs-ground_effect_reference
  ShapeNo = 3
  DrawFilledPolygon
  Foliage  '  prop_loss_cum=prop_loss_cum+foliage_attenuation
  ShapeNo = 4
  DrawFilledPolygon
  Barrier   '  prop_loss_cum=prop_loss_cum+barrier_attenuation
  ShapeNo = 5
  DrawFilledPolygon
  ShapeNo = 1: For I = 0 To 23: prop_loss_cum_4graph(I) = target_spectrum(I): Next I:  TraceXY   ' thick solid cyan
  ShapeNo = 2: For I = 0 To 23: prop_loss_cum_4graph(I) = target_spectrum(I) + prop_loss_cum(I): Next I: TraceXY ' thick solid green
  ShapeNo = 9: For I = 0 To 23: prop_loss_cum_4graph(I) = listener_spectrum(I): Next I: TraceXY ' thin solid brown
  ShapeNo = 3: For I = 0 To 23: prop_loss_cum_4graph(I) = listener_spectrum(I): Next I: TraceXY ' thick dash green
  ShapeNo = 7: For I = 0 To 23: prop_loss_cum_4graph(I) = background_noise_spec(I): Next I: TraceXY ' thick solid red
  ShapeNo = 6: For I = 0 To 23: prop_loss_cum_4graph(I) = hearing_threshold(I): Next I: TraceXY ' thin dash red
  ShapeNo = 6: For I = 0 To 23: prop_loss_cum_4graph(I) = mod_background_noise(I): Next I: TraceXY ' thin dash red
  ShapeNo = 9: For I = 0 To 23: prop_loss_cum_4graph(I) = max_modBackground_thresh(I): Next I: TraceXY  ' thick dash red
  ShapeNo = 4: For I = 0 To 23: prop_loss_cum_4graph(I) = max_modBackground_thresh(I): Next I: TraceXY  ' thick dash red
End Sub

