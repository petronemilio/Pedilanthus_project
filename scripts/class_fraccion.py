
def gcd(m,n):
    while m%n != 0:
        oldm = m
        oldn = n

        m = oldn
        n = oldm%oldn
    return n


class Fraction:
    
    def __init__(self,top,bottom):
	self.num = top
        self.den = bottom
    
    def __str__(self):
        return str(self.num)+"/"+str(self.den)
 
    
    def __add__(self,otherfraction):
        
        newnum = self.num*otherfraction.den + otherfraction.num*self.den
        newden = self.den * otherfraction.den
        common = gcd(newnum,newden)
        return Fraction(newnum//common,newden//common)

    def __sub__(self,otherfraction):
         newnum = self.num*otherfraction.den - self.den*otherfraction.num
         newden = self.den * otherfraction.den
         common = gcd(newnum,newden)
         return Fraction(newnum//common,newden//common)

    def __mul__(self,otherfraction):
         newnum = self.num*otherfraction.num 
         newden = self.den * otherfraction.den
         common = gcd(newnum,newden)
         return Fraction(newnum//common,newden//common)
    def __truediv__(self,otherfraction):
         newnum = self.num*otherfraction.num 
         newden = self.den * otherfraction.den
         common = gcd(newnum,newden)
         return Fraction(newnum//common,newden//common) 
  
    def __eq__(self, other):
        firstnum = self.num * other.den
        secondnum = other.num * self.den
        return firstnum == secondnum



myfraction = Fraction(3,5)
