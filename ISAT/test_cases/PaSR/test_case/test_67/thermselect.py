species = ["C2H4", "O2", "CO2", "H2O", "OH", "H", "CO", "H2", "O", "HO2", "CH3", "C4H2", "CH4", "C2H2", "HCCO", "CH2CHO", "CH2", "C2H5", "CH3O", "C2H3", "C2H6", "H2O2", "C3H6", "C2O", "CH", "aC3H5", "C2H", "CH2*", "H2CC", "C", "CH2O", "iC4H3", "C4H6", "HCO", "CH2OCH2", "CH2CHCO", "AR", "N2"]

thermfile = open('therm.dat', 'r')

newtherm = ''
ln = 0
sample = False
for line in thermfile:
    words = line.split()
    if words[0] in species:
        sample = True
    if sample:
        ln += 1
        newtherm += line
    if ln == 4:
        sample = False
        ln = 0


newthermfile = open('newtherm.dat', 'w')
newthermfile.write(newtherm)
newthermfile.close()
thermfile.close()
