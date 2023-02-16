

#Runs parameters through set formulas and changes GEF values
def ComplexDiv(A, B, C, D, E, F):
    G = 1 / (C * C + D * D)
    E = (A * C + B * D) * G
    F = (B * C - A * D) * G
    return 1

def main():
    #Test Values; Remove if going to merge with main
    ComplexDiv(5, 6, 4, 7, 9, 2)
    

if __name__=='__main__':
    main()