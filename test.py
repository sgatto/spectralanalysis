def timeFilter(t,tmax):
    return math.exp(-(((t-tmax*0.5)/(0.35*tmax))**4))
