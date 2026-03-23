import csv
from helper_functions.matrix_solving import gauss_elim_p_piv

file  = "SV&SH_data.csv"

def read_csv(file):
    with open(file,'r') as file:    #
        reader = csv.reader(file)
        sv = []
        sh = []
        for SV,SH in reader: 
            sv.append(float(SV))
            sh.append(float(SH))
    return sv, sh

def spline(sv, sh, 
sv, sh = read_csv(file)


#print(sv)
#print(sh)
