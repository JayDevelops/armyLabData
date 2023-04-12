def ground_effect(freq, prop_loss_cum, prop_loss_indiv, ground_effect_ref):
    windspeed = data_dict.det_cons['wind_speed'] #this looks a bit different than the original, not sure if thats ok
    Iwthr1 = 0
    
    
    if (Iwthr1 == 0):
        ground_effect = Ingard(freq)
        if (windspeed >= 0):
            for i in range(24):
                prop_loss_indiv[i] = ground_effect[i] - ground_effect_ref[i]
                if (prop_loss_indiv[i] == 0):
                    prop_loss_indiv[i] = -0.01
                prop_loss_cum[i] = prop_loss_cum[i] + prop_loss_indiv[i]
        else:
            for i in range(24):
                if (prop_loss_indiv[i] == 0):
                    prop_loss_indiv[i] = -0.001
                prop_loss_indiv[i] = ground_effect[i] - ground_effect_ref[i]
    print('ground_effect', prop_loss_cum)
    return prop_loss_cum, prop_loss_indiv
