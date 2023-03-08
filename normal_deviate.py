import math

def normal_deviate(p):
    t = math.sqrt(-2 * math.log(p))
    return t - (2.30753 + t * 0.27061) / (1 + t * (0.99229 + t * 0.04481))