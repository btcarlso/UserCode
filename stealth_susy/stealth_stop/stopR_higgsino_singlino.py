#!/usr/bin/env python

from math import *
import subprocess

# basic mass parameters
mtop = 174.0
mZ = 91.2
mh = 125.0

#phase space function
def qPS(x,y):
  return 1 + x*x + y*y - 2.0*x - 2.0*y - 2.0*x*y

# left-handed stop branching ratios
# returns [H_1 t_R, H_2 t_R, H+ b_R] branching ratios
def stopLbr(mstop, mu):
  if(mstop > mtop + mu):
    return [0.5, 0.5, 0.0]
  return [0.0, 0.0, 1.0]

#left-handed sbottom branching ratios
# returns [H- t_R, H_1 b_R, H_2 b_R]
def sbottomLbr(msbottom, mu):
  if(msbottom > mtop + mu):
    return [1.0, 0.0, 0.0]
  return [0.0, 0.5, 0.5]

# right-handed stop branching ratios
# returns [H_1 t_L, H_2 t_L, H+ b_L] branching ratios
def stopRHbr(mstop, mu):
  if(mstop > mtop + mu):
    bA = qPS(mtop**2.0 / mstop**2.0, mu**2.0 / mstop**2.0)**0.5 * (mstop**2.0 - mtop**2.0 - mu**2.0)
    bB = qPS(0.0, mu**2.0 / mstop**2.0)**0.5 * (mstop**2.0 - mu**2.0)
    return [0.5*bA/(bA + bB), 0.5*bA/(bA+bB), bB/(bA+bB)]
  return [0.0, 0.0, 1.0]

# higgsino branching ratios
# assumes higgsino -> Z + singlino is open
# mS denotes m_singlino
def higgsinobr(mu, mS, tanbeta, alpha):
  beta = atan(tanbeta)
  if(mu > mS + mh):
    qZ = qPS(mZ**2.0/mu**2.0, mS**2.0/mu**2.0)
    qh = qPS(mh**2.0/mu**2.0, mS**2.0/mu**2.0)
    GamplusZ = (cos(beta) + sin(beta))**2.0 * qZ**0.5/(mu*(mS - mu)**2.0) * (0.5 * (mu**4.0 * qZ + 3*mZ**2.0 * (mu**2.0 + mS**2.0 - mZ**2.0)) - 3 * mu * mS * mZ**2.0)
    GamminusZ = (cos(beta) - sin(beta))**2.0 * qZ**0.5/(mu*(mS + mu)**2.0) * (0.5 * (mu**4.0 * qZ + 3*mZ**2.0 * (mu**2.0 + mS**2.0 - mZ**2.0)) + 3 * mu * mS * mZ**2.0)
    GamplusH = 0.5 * (cos(alpha) - sin(alpha))**2.0 * qh**0.5 / mu * (mu**2.0 - 2.0 * mu * mS + mS**2.0 - mh**2.0)
    GamminusH = 0.5 * (cos(alpha) + sin(alpha))**2.0 * qh**0.5 / mu * (mu**2.0 + 2.0 * mu * mS + mS**2.0 - mh**2.0)
    return [GamminusZ / (GamminusZ + GamminusH), GamminusH / (GamminusZ + GamminusH), GamplusZ/(GamplusZ + GamplusH), GamplusH/(GamplusZ + GamplusH)]
  # else
  return [1.0, 0.0, 1.0, 0.0]

# write the SLHA file 
stopmasses = [150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0] 
higgsinomasses = [200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0]

counter = 0
maxtowrite = 10
for mst in stopmasses:
  for mhi in higgsinomasses:
    if(mst - mhi < 5.0):
      continue
    counter += 1
    slhaname = "StHiS_model_stop"+str(mst)+"_higgsino"+str(mhi)+".slha"
    subprocess.call("cp slha_top.txt "+slhaname, shell=True)
    slhaout = open(slhaname, 'a')
    slhaout.write("      1000006    "+str(mst)+"  # m(stop)\n")
    slhaout.write("      1000022    "+str(mhi)+"  # m(chi0_1)\n")
    slhaout.write("      1000023    "+str(mhi)+"  # m(chi0_2)\n")
    slhaout.write("      1000024    "+str(mhi)+"  # m(chi+_1)\n")
    slhaout.write("\n")
    # write the stop decay table
    stopbrs = stopRHbr(mst, mhi)
    slhaout.write("#         ID            Width\n")
    slhaout.write("DECAY   1000006     1.00000000E+00   # t1 decays \n")
    strchi01t = "%.8E" % stopbrs[0]  # format the branching ratio strings nicely
    strchi02t = "%.8E" % stopbrs[1]
    strchipb  = "%.8E" % stopbrs[2]
    slhaout.write("     "+strchi01t+"    2     1000022         6        # BR(t1 -> chi0_1 t)\n")
    slhaout.write("     "+strchi02t+"    2     1000023         6        # BR(t1 -> chi0_2 t)\n")
    slhaout.write("     "+strchipb +"    2     1000024         5        # BR(t1 -> chi+_1 b)\n")
    # write the gluino decay table
    #slhaout.write("#         ID            Width\n")
    #slhaout.write("DECAY   1000021     1.00000000E+00   # gluino decays \n")
    #if(mgl - mst > mtop):
    #  slhaout.write("     5.00000000E-01    2    -1000006         6        # BR(gluino -> t t1~)\n")
    #  slhaout.write("     5.00000000E-01    2     1000006        -6        # BR(gluino -> t~ t1)\n")
    #else:
    #  # 3-body; insert check if this fails?
    #  slhaout.write("     5.00000000E-01    3    -1000006         5        24        # BR(gluino -> w+ b t1~)\n")
    #  slhaout.write("     5.00000000E-01    3     1000006        -5       -24        # BR(gluino -> w- b~ t1)\n")
    # write the chi0_1 decay table
    # fixing singlino mass = 100 GeV, tan beta = 10.0, alpha = pi/2 - beta
    hibr = higgsinobr(mhi, 100.0, 10.0, pi/2.0 - atan(10.0))
    h1zs = "%.8E" % hibr[0]
    h1hs = "%.8E" % hibr[1]
    h2zs = "%.8E" % hibr[2]
    h2hs = "%.8E" % hibr[3]
    slhaout.write("#         ID            Width\n")
    slhaout.write("DECAY   1000022     1.00000000E-04   # chi0_1 decays \n")
    slhaout.write("     "+h1zs+"    2     5000001        23        # BR(chi0_1 -> Z singlino)\n")
    slhaout.write("     "+h1hs+"    2     5000001        25        # BR(chi0_1 -> h singlino)\n")
    # write the chi0_2 decay table
    slhaout.write("#         ID            Width\n")
    slhaout.write("DECAY   1000023     1.00000000E-04   # chi0_2 decays \n")
    slhaout.write("     "+h2zs+"    2     5000001        23        # BR(chi0_2 -> Z singlino)\n")
    slhaout.write("     "+h2hs+"    2     5000001        25        # BR(chi0_2 -> h singlino)\n")
    # write the chi+_1 decay table
    slhaout.write("#         ID            Width\n")
    slhaout.write("DECAY   1000024     1.00000000E-04   # chi+_1 decays \n")
    slhaout.write("     1.00000000E-00    2     5000001        24        # BR(chi+_1 -> w+ singlino)\n")
    slhaout.close()
    subprocess.call("cat slha_bottom.txt >> "+slhaname, shell=True)
    if(counter >= maxtowrite):
      raise SystemExit
