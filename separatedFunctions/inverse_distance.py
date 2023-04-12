def inverse_distance(detection_distance, prop_loss_cum, prop_loss_indiv):
    m_measure_distance = data_dict.given_cons['m_measure_distance']

    Log10Div10 = 0.230258509
    TenDivLog10 = 1 / Log10Div10
    #print(m_measure_distance)
    #print(detection_distance)
    #D4 should be mic distance source m. changes depending on what macro is used. Can also be detect distance
    p_inv = 2 * TenDivLog10 * math.log(m_measure_distance / detection_distance) #Inverse distance loss from D3 to D4. D3 is the mic distance from the target from excel.
    #print(p_inv)
    for i in range(24):
        if (p_inv == 0):
            prop_loss_indiv[i] = -0.001
        prop_loss_cum[i] = prop_loss_indiv[i]
    print('Inverse_distance', prop_loss_cum)
    
    return prop_loss_cum, prop_loss_indiv
