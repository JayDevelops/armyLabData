

#Runs parameters through set formulas and changes EF values
def ComplexMul(A, B, C, D, E, F):
    E = A * C - B * D
    F = A * D + B * C
    return 1

def main():
    #Test Values; Remove if going to merge with main
    ComplexMul(5, 6, 4, 7, 9, 2)
    

if __name__=='__main__':
    main()