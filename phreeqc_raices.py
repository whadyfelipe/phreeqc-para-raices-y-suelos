#function phreeqq for raices
#https://www.phreeqpy.com/examples.html
#https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/html/final-54.html

import phreeqpy.iphreeqc.phreeqc_dll as phreeqc_mod

species={'pH':7.0,'pe':12.5,'Ca':0.6,'Cl':1.2,'Na':0,'K':0,'N(5)':0,\
'CB':-1.2e-5,'TOT_H':1.1e2,'TOT_O':5.5e1,'CaX2':5.32e-4,'KX':1.76e-5,\
'NaX':1.76e-5}

def phreeq_raices(species):
    textinput='TITLE Example 11.--Transport and ion exchange. \n'
    textinput+='SOLUTION           1 \n'
    textinput+='units            mmol/kgw \n'
    textinput+='temp             25.0 \n'
    textinput+='pH               '+str(species['pH'])+'    charge \n'
    textinput+='pe               '+str(species['pe'])+'    O2(g)   -0.68 \n'
    textinput+='Ca               '+str(species['Ca'])+'\n'
    textinput+='Cl               '+str(species['Cl'])+'\n'
    textinput+='Na               '+str(species['Na'])+'\n'
    textinput+='K               '+str(species['K'])+'\n'
    textinput+='N(5)              '+str(species['N(5)'])+'\n'
    textinput+='CB               '+str(species['CB'])+'\n'
    textinput+='TOT_H            '+str(species['TOT_H'])+'\n'
    textinput+='TOT_O            '+str(species['TOT_O'])+'\n'
    textinput+='EXCHANGE             1'+'\n'
    textinput+='equilibrate 1'+'\n'
    textinput+='X  0.0011'+'\n'
    textinput+='CaX2          '+str(species['CaX2'])+'\n'
    textinput+='KX            '+str(species['KX'])+'\n'
    textinput+='NaX           '+str(species['NaX'])+'\n'
    textinput+='SELECTED_OUTPUT'+'\n'
    textinput+='-reset           false'+'\n'
    textinput+='-step'+'\n'
    textinput+='-high_precision          true'+'\n'
    textinput+='-ph          true'+'\n'
    textinput+='USER_PUNCH '+'\n'
    textinput+='-headings  charge    H   O  Ca  Cl  Na  K N(5) pH CaX2 KX NaX'+'\n'
    textinput+='10 PUNCH charge_balance'+'\n'
    textinput+='20 PUNCH TOTMOLE("H"), TOTMOLE("O"),TOTMOLE("Ca")'+'\n'
    textinput+='30 PUNCH TOTMOLE("Cl"),TOTMOLE("Na"),TOTMOLE("K"),TOTMOLE("N(5)")'+'\n'
    textinput+='40 PUNCH -LA("H+"),MOL("CaX2"),MOL("KX"),MOL("NaX")'+'\n'
    textinput+='END'+'\n'

    #Initialize IPhreeqc module
    phreeqc = phreeqc_mod.IPhreeqc()
    phreeqc.load_database(r"phreeqc.dat")
    phreeqc.run_string(textinput)
    components = phreeqc.get_component_list()    
#    phc_string = 'RUN_CELLS; -cells 1\n'
#    phreeqc.run_string(phc_string)
    output = phreeqc.get_selected_output_array()
    
    return output

output=phreeq_raices(species)