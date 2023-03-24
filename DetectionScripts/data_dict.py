#Dictionaries to store parameters referenced in 'Model' sheet
#For now it is set at constants, these will be adjusted based on user/program input in the future

meas_cons = {
    'mic_height_meas':  1.2,
    'celsius_degrees_meas': 15,
    'relative_humid_percent_meas': 70,
    'sigma_meas': 200,
    'em2_meas': 0.0000006,
}

det_cons = {
    'source_height_det': 1.2,
    'listener_height_det': 1.2,
    'celsius_degrees_det': 15,
    'relative_humid_percent_det': 70,
    'sigma_det': 200,
    'em2_det' : 0.0000006,
    'wind_speed': 0,
    'wind_flag': False,
    'wind_direction': 'UP',
    'barrier_on': False,
    'barrier_dist': 1,
    'barrier_height': 2,
    'foliage_on': False,
    'foliage_dist': 1,
    'foliage_depth': 100,
    'leaf_areapervol': 0.5,
    'leaf_width': 3.2,
}


lis_cons = {
    'observer_efficiency': 0.40,
    'hit_prob': 0.50,
    'false_alarm_rate': 0.01,
    'd_stat': 2.33, #this is calculated using Dprime
    'is_expert': True,
    'is_binaural': False,
}

given_cons = {
    'm_measure_distance': 250,
    'lis_dist_given': 30,
    'dba_given': 55,
    'trg_name': 'M60 Tank idling at 30 meters',
    'bkg_name': 'Urban',
    'hth_name': 'ISO Std',
}
