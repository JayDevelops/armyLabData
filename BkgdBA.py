# Changed my orginal function to follow jonathan's listenerBDA template  for consistency
def BkgdBA():
     E = float()
     E = 0.0
     
     for I in range(0,23):
         E = E + math.exp(Log10Div10 * (S6[I] + S10[I]))
         
     bkgdba_return = TenDivLog10 * math.log(E)
     
     return bkgdba_return    