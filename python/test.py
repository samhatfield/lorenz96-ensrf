from numpy.random import randn
from math import sqrt

def ar(last, phi, sigma):
    z = sigma * randn()

    return phi * last + sqrt(1-phi**2) * z

series = [0.0] * 3
params = [0.75, 0.997, 0.9997]

for _ in range(200):
    print('{0}\t{1}\t{2}'.format(*tuple(series)))
    series = [ar(l, p, 0.126) for l, p in zip(series, params)]
