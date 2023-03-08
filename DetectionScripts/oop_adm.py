import math


def inverse_distance(measure_dist: float, detection_dist: float) -> list:
    # 10 divided by log 10 static, and log of measure distance divided by detection distance log in base 10
    ten_divided_by_log_10, log_m_dist_by_d_dist = (10 / math.log(10, 10)), math.log(measure_dist / detection_dist, 10)

    inverse_distances = []

    for x in range(0, 23):
        inv_dist = 2 * ten_divided_by_log_10 * log_m_dist_by_d_dist
        if inv_dist == 0:
            inverse_distances.append(inv_dist)
        else:
            inverse_distances.append(-0.001)

    return inverse_distances

def reference_calc(sigma: float, source_height_condition: float,
                   measure_dist: float, mic_height_coord: float,
                   em2_ref: float, ground_effect_ref: float) -> list:
    arr = []
    for x in range(0, 23):
        new_num = Ingard() + ansi_humidity()
        arr.append(Ingard() + AnsiHumidity())

    return arr
