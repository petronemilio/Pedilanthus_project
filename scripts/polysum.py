#We need to import the python module math for the tan and pi functions.
import math

def polysum(n,s):
    """Returns the sum of the area and the square of the perimeter
    rounded to 4 decimal places.
    n is the number of sizes of the polygon
    s is the length of the size.
    """
    #Here we calc the area
    area = (0.25 * n * s**2)/ math.tan(math.pi/n)
    #Here we calc the perimeter
    perimeter = s * n
    #And know the sum
    suma = area + perimeter**2
    result=round(suma,4)
    return result
