import pandas
import numpy as np

def main():
    bestandLijst = ['glucose-receptor-homosapiens', 'GPCR-homosapiens', 'HIV1 env-eiwit.fasta', 'HIV1-ivg-eiwit.txt', 'HIV2-gp-eiwit', 'HIV2-ivg-eiwit.txt', 'homosapiens-huishoudgen.fasta', 'receptor-kinase-homosapiens', 'SIV-gp-eiwit', 'SIV-ivg-eiwit.txt', 'SIVmnd2-gp-eiwit.txt', 'SIVmnd2-ivg-eiwit.txt'] 
    alleData = []                                                                                                       #maak eel lege lijst om sublijsten aan toe te voegen
    for i in range(len(bestandLijst)):                                                                                  #zo lang i in de range van de lengte van bestandLijst zit... 
        file = open(bestandLijst[i])                                                                                    #open het bestand met de bestandsnaam op positie met waarde i
        print(bestandLijst[i])                                                              
        completeSeq = lezen(file)                                                                                       #roep de lezen functie aan, geef het bestand mee en verkrijg hieruit alle coderende sequentie uit het bestand
        codonLijst, codonGebLijst = tellen(completeSeq)                                                                 #roep de tellen functie aan, geef de sequentie mee en verkrijg hieruit een lijst met aminozuren en de frequentie hiervan
        totaal = 0                                                                                                      #geef totaal waarde 0
        huidigBestandData = []                                                                                          #maak een lege lijst voor de data in het huidige bestand
        for i in codonGebLijst:                                                                                         #zo lang i in de codonGebLijst lijst zit...
            totaal+=i                                                                                                   #verhoog totaal met wat er op positie i dstaat in codonGebLijst
        huidigBestandData.append(totaal)                                                                                #voeg het totaal toe aan de huidigBestandDataLijst
        huidigBestandData.append(round(float(codonGebLijst[4]/(totaal+1)*100),3))                                       #voegt frequentie cysteine toe aan de huidigBestandDataLijst
        huidigBestandData.append(round(float(codonGebLijst[17]/(totaal+1)*100),3))                                      #" tryptophan, " hydrofoob, " hydrofiel
        huidigBestandData.append(round(float((codonGebLijst[17]+codonGebLijst[0]+codonGebLijst[10]+codonGebLijst[19]+codonGebLijst[11]+codonGebLijst[5]+codonGebLijst[18])/(totaal+1)*100),3))
        huidigBestandData.append(round(float((codonGebLijst[2]+codonGebLijst[15]+codonGebLijst[6]+codonGebLijst[18]+codonGebLijst[16]+codonGebLijst[1]+codonGebLijst[3]+codonGebLijst[12]+codonGebLijst[9])/(totaal+1)*100),3))
        codonGebLijstMax, codonGebLijstMin = lijsten(codonGebLijst)
        tel = 0                                                                                         #zet een teller op 0
        for i in range(0,3):                                                                            #zo lang i in de range van 0 tot 3 zit...
            huidigBestandData.append(codonLijst[codonGebLijstMax.index((max(codonGebLijstMax)))])       #voeg de hoogste waarde toe aan de huidigBestandData lijst
            codonGebLijstMax[codonGebLijstMax.index((max(codonGebLijstMax)))] = 0                       #overschrijf de hoogste waarde met 0
            tel+=1                                                                                      #verhoog de teller met 1
        tel = 0                                                                                         #zet een teller op 0
        for i in range(0,3):                                                                            #zo lang i in de range va 0 tot 3 zit...
            huidigBestandData.append(codonLijst[codonGebLijstMin.index((min(codonGebLijstMin)))])       #voeg de laagste waarde toe aan de huidigBestandData lijst
            codonGebLijstMin[codonGebLijstMin.index((min(codonGebLijstMin)))] = max(codonGebLijstMin)+1 #overschrijf de kleinste waarde met een hoge waarde 
            tel+=1                                                                                      #verhoog de teller met 1
        alleData.append(huidigBestandData)                                                              #voeg de lijst huidigBestandData toe aa alleData
    tabel(alleData,bestandLijst)                                                                        #roep de tabel functie aan den geef alleData en bestandLijst mee

def lijsten(codonGebLijst): #deze functie maakt 2 identieke lijsten waar verschillende posities mogen worden gewijzigd
    codonGebLijstMax = []                           #maak een lege lijst waat de minimale waarden mogen worden overschreven
    codonGebLijstMin = []                           #maak een lege lijst waar de maximale waarden mogen worden overschreven
    for i in codonGebLijst:                         #zo lang i in codonGebLijst zit...
        codonGebLijstMax.append(i)                  #voeg de waarde van i toe aan codonGebLijstMin
        codonGebLijstMin.append(i)                  #voeg de waarde van i toe aan codonGebLijstMin
    return codonGebLijstMax, codonGebLijstMin 

def tabel(alleData, bestandLijst): #deze functie geeft de data weer
    dataVolg = ['               totale lengte', 'Cysteine frequentie','Tryptofaan frequentie','Percentage hydrofoob','Percentage hydrofiel','Meest voork. 1','Meest voork. 2','Meest voork. 3','Minst voork. 1','Minst voork. 2','Minst voork. 3' ]
    alledata = np.array(alleData)                                                   #gebruik panda en numpy om data te weergeven
    tabel = pandas.DataFrame(alleData, bestandLijst, dataVolg)                      #geef pandas de juiste data
    print(tabel)
#- - - - - - - - - - - - - - - - - - - 
#hydrofiel = N, S, Q, Y, T, R, D, K, H
#hydrofoob = W, A, I, V, L, F, Y
#- - - - - - - - - - - - - - - - - - -
        
def lezen(file): #deze functie leest het bestand in en voegt alle coderende stukken samen
    completeSeq = ''                            #maak een lege string waar de sequentie in gaat worden samengesteld
    for line in file.readlines():               #loop door het bestand...
        if line[0] != '>':                      #als het eerste teken in de line geen > is...
            line = line.rstrip()                #strip de line
            completeSeq+=line                   #voeg de line toe aan completeSeq
    return completeSeq                      

def tellen(completeSeq): #deze functie telt de frequentie van alle aminozuren 
    codonLijst = ['A', 'R', 'N', 'D', 'C', 'F', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'P', 'S', 'T', 'W', 'Y', 'V']   #codonLijst is een lijst met de letters van alle aminozuren evenredig aan codonGebLijst 
    codonGebLijst = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]                                        #maak een lijst waar de frequenties in komen te staan
    for i in completeSeq:                                                                                               #zo lang i in de sequentie zit...
        codonGebLijst[codonLijst.index(i)]+=1                                                                           #verhoog de waarde die staat op de positie van waar in codonLijst deze letter staan met 1
    print(codonGebLijst)                                                                                                
    return codonLijst, codonGebLijst                                    


main() #roep de main functie aan 
