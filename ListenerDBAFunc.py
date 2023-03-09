# adds the previous E value in the function with the exp function of Log10Div10 * 3 separate array values
def ListenerdBA():
    E = float()
    E = 0.0
    for I in range(0, 23):
        E = E + math.exp(Log10Div10 * (S1[I] + S3[I] + S10[I]))
    listener_dba_return = TenDivLog10 * math.log(E)
    return listener_dba_return
