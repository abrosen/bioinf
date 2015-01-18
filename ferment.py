#!/usr/bin/python
class Malt(object):
    def __init__(self,name,constant,ferment):
        self.name =name 
        self.constant = constant
        self.ferment =  ferment
    
    def getWeight(self,density,volume):
        return ((density -1000)*volume)/self.constant
    
    def getABV(self,density):
        return (density - 1000)*self.ferment/7.45
    
    def getFinalDensity(self,density):
        return (density - 1000)*(1.-self.ferment)+1000
    def __str__(self):
        return self.name

def main():
    sugars = []
    sugars.append(Malt("malt extract", 303, 0.62))
    sugars.append(Malt("table sugar", 375,1.00))
    sugars.append(Malt("chocolate malt", 200, 0.30))
    density = input("Enter the desired density: ")
    liters = input("How many liters would you like to make? ")
    try:
        density = float(density)
        liters = float(liters)
    except:
        print "Your input is bad and you should feel bad!"
        exit(0)
    for s in sugars:
        print"For", s.name, "You need", s.getWeight(density,liters), "kg which yields a final denisity of",s.getFinalDensity(density),"and",s.getABV(density), "abv"

if __name__ == "__main__":
    main()

