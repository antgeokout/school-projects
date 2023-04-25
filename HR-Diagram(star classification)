# Import necessary packages (numpy, pyplot, astropy)
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

# Create an Astropy Table with data from Vizier
# mytable is the Gaia catologue cross matched with 2mass
mytable = Table.read("gaiaX2mass-2.vot", format='votable')
# Gaia catalogue of all stars 50pc from the Earth
orgtable = Table.read("Gaia-50pc-all.vot", format='votable')

# Create variables for your magnitudes, parallax, Temperature and calculate color (R - K mag)
kmag = mytable['Kmag']
gmag = mytable['Gmag']
Rmag = mytable['RPmag']
color = Rmag-kmag
plx = mytable['Plx']
e_plx = mytable['e_Plx']
Teff = mytable['Tefftemp']

# more magnitudes to check against cross matched values using 2mass
gmagorg = orgtable['Gmag']
Rmagorg = orgtable['RPmag']
colororg = gmagorg - Rmagorg
plx_org = orgtable['Plx']
e_plx_org = orgtable['e_Plx']

# absolute magnitude formula
absmag = gmag + (5*(np.log10(plx/1000))) + 5
absmag_org = gmagorg + (5*(np.log10(plx_org/1000))) + 5

# configure plot of abs mag vs. color
fig1, ax = plt.subplots(figsize=(12,6))
ax.plot(colororg, absmag_org, marker='.', linestyle='none')
ax.set_ylim(25,-5)
plt.show()

# error in parallax
fig1, ax = plt.subplots(figsize=(12,6))
ax.plot(e_plx_org/plx_org, plx_org, marker='.', linestyle='none')
plt.show()

# Create empty lists to add values with low err in plx
plxcorr = []
e_plx_cor = []
gmag_cor = []
kmag_cor = []
Rmag_cor = []

# for each of our variables, only keep ones with err less than 0.01
ratio = e_plx / plx
for i in range(len(plx)):
    if ratio[i] < 0.01:
        plxcorr.append(plx[i])
        e_plx_cor.append(e_plx[i])
        gmag_cor.append(gmag[i])
        kmag_cor.append(kmag[i])
        Rmag_cor.append(Rmag[i])

# convert our new lists into numpy arrays
plxcorr = np.array(plxcorr)
e_plx_cor= np.array(e_plx_cor)
gmag_cor = np.array(gmag_cor) 
kmag_cor = np.array(kmag_cor)
Rmag_cor = np.array(Rmag_cor)

# calculate abs mag with new values
ratio_cor = e_plx_cor / plxcorr
color_cor = Rmag_cor-kmag_cor
absmag_cor = gmag_cor + (5*np.log10(plxcorr/1000)) + 5 

# plot of err
fig1, ax = plt.subplots(figsize=(12,6))
ax.plot(ratio_cor, plxcorr, marker='.', linestyle='none')
plt.show()

# HR diaram with updated values
fig1, ax = plt.subplots(figsize=(12,6))
ax.plot(color_cor, absmag_cor, marker='.', linestyle='none')
ax.set_ylim(25,-5)
plt.show()

'''
Now do the same with our original non-cross-matched values
'''
plxcorr_org = []
e_plx_cor_org = []
gmag_cor_org = []
Rmag_cor_org = []

ratio_org = e_plx_org / plx_org
for i in range(len(plx_org)):
    if ratio_org[i] < 0.01:
        plxcorr_org.append(plx_org[i])
        e_plx_cor_org.append(e_plx_org[i])
        gmag_cor_org.append(gmagorg[i])
        Rmag_cor_org.append(Rmagorg[i])

plxcorr_org = np.array(plxcorr_org)
e_plx_cor_org= np.array(e_plx_cor_org)
gmag_cor_org = np.array(gmag_cor_org) 
Rmag_cor_org = np.array(Rmag_cor_org)


ratio_cor_org = e_plx_cor_org / plxcorr_org
color_cor_org = gmag_cor_org-Rmag_cor_org
absmag_cor_org = gmag_cor_org + (5*np.log10(plxcorr_org/1000)) + 5 

fig1, ax = plt.subplots(figsize=(12,6))
ax.plot(ratio_cor_org, plxcorr_org, marker='.', linestyle='none')
plt.show()

fig1, ax = plt.subplots(figsize=(12,6))
ax.plot(color_cor_org, absmag_cor_org, marker='.', linestyle='none')
ax.set_ylim(25,-5)
plt.show()

'''

loose description of each star's classification by their temperature to magnitude range

B = < -2
temp = 16,400 - 30,000

A = [-2,0.2]
temp = 8,620 - 16,400

F = [0.2,0.5]
temp = 6,540 - 8,620

G = [0.5,1.2]
temp = 5,610 - 6,540

K = [1.2, 2]
temp = 4,410 - 5,610

M = > 2
temp = 2,650 - 4,410
'''

# empty list for each star type by color
bstar_c = []
astar_c = []
fstar_c = []
gstar_c = []
kstar_c = []
mstar_c = []

# empty list for each star type by abs mag
bstar_m = []
astar_m = []
fstar_m = []
gstar_m = []
kstar_m = []
mstar_m = []

# for loop sorting through each data point and classifying it based on it's color
for i in range(len(absmag_cor)):
    if 2.5*color_cor[i] - 2 < absmag_cor[i] < 2.5*color_cor[i] + 5:
        if -2 < color_cor[i]:
            bstar_c.append(color_cor[i])
            bstar_m.append(absmag_cor[i])
        if -2 < color_cor[i] < 0.2 :
            astar_c.append(color_cor[i])
            astar_m.append(absmag_cor[i])
        if 0.2 < color_cor[i] < 0.5 :
            fstar_c.append(color_cor[i])
            fstar_m.append(absmag_cor[i])
        if 0.5 < color_cor[i] < 1.2 :
            gstar_c.append(color_cor[i])
            gstar_m.append(absmag_cor[i])
        if 1.2 < color_cor[i] < 2 :
            kstar_c.append(color_cor[i])
            kstar_m.append(absmag_cor[i])
        if color_cor[i] > 2 :
            mstar_c.append(color_cor[i])
            mstar_m.append(absmag_cor[i])
            
# convert lists into arrays
np.array(bstar_c)
np.array(astar_c)
np.array(fstar_c)
np.array(gstar_c)
np.array(kstar_c)
np.array(mstar_c)

np.array(bstar_m)
np.array(astar_m)
np.array(fstar_m)
np.array(gstar_m)
np.array(kstar_m)
np.array(mstar_m)

# lines indicating the cuttoff point for each star type
x1, y1 = [-2,-2] , [25,-5]
x2, y2 = [-0.5,-0.5] , [25,-5]
x3, y3 = [0.5,0.5] , [25,-5]
x4, y4 = [1,1] , [25,-5]
x5, y5 = [2,2] , [25,-5]

# mag 1 and 2 cuttoff white dwarf and red giants
mag1 = 2.5*color_cor - 2
mag2 = 2.5*color_cor + 5

# cuttoff for each star type
magA = -2*color_cor + 3
magF = -2*color_cor + 5
magG = -2*color_cor + 7
magK = -2*color_cor + 9
magM = -2*color_cor + 11

# plotting each star type
fig1, ax = plt.subplots(figsize=(12,6))
ax.plot(astar_c,astar_m, marker='.', linestyle='none')
ax.plot(fstar_c,fstar_m, marker='.', linestyle='none')
ax.plot(gstar_c,gstar_m, marker='.', linestyle='none')
ax.plot(kstar_c,kstar_m, marker='.', linestyle='none')
ax.plot(mstar_c,mstar_m, marker='.', linestyle='none')
#ax.plot(color_cor,absmag_cor , marker='.', linestyle='none')

#imaginary lines of cuttoff for star types
#ax.plot(color_cor,mag1 ,color='yellow')
#ax.plot(color_cor,mag2,color='green')
#ax.plot(x1,y1)
#ax.plot(x2,y2)
#ax.plot(x3,y3)
#ax.plot(x4,y4)
#ax.plot(x5,y5)

# graph display
ax.yaxis.get_ticklocs(minor=True)
ax.minorticks_on()
ax.set_ylim(25,-5)
ax.set_xlim(-5,20)
ax.set_xlabel('color (RPmag-Kmag)')
ax.set_ylabel('Mg (mag)')
ax.set_title('HR Diagram (MS-Subtype)')
ax.legend(['A-star','F-star','G-star','K-star','M-star'])
plt.show()

fig1.savefig("HR_diagram.png",bbox_inches='tight')

# Import tabulate
from tabulate import tabulate
# create a list of each star type
starmag = [astar_m, fstar_m, gstar_m, kstar_m, mstar_m]

# empty list for total number of MS stars
num = []

# empty lists for adding mass of each star
astar_mass_log = []
fstar_mass_log = []
gstar_mass_log = []
kstar_mass_log = []
mstar_mass_log = []

# list of number of each subtype
masslog = [astar_mass_log,fstar_mass_log,gstar_mass_log,kstar_mass_log,mstar_mass_log]


#Calculate the mass of each star from each subsection
def masscalc(P1,P2):
    for i in range(len(P1)):
        if P1[i] < 9.58:
            P2.append(10**(-0.15*P1[i] + 1.05))
        else:
            P2.append(10**(-0.07*P1[i] + 0.36))
            
masscalc(astar_m,astar_mass_log)
masscalc(fstar_m,fstar_mass_log)
masscalc(gstar_m,gstar_mass_log)
masscalc(kstar_m,kstar_mass_log)
masscalc(mstar_m,mstar_mass_log)

#total mass of stars in each subsection
total_mass = []
for i in range(len(masslog)):
    total_mass.append(sum(masslog[i]))

#total mass of all stars
all_mass = sum(total_mass)

#percentage of mass in each subtype
pmas = []
for i in range(len(total_mass)):
    pmas.append(total_mass[i]/all_mass)

#number of stars in each subsection
for i in range(len(starmag)):
    num.append(len(starmag[i]))

#total number of stars in MS
numsum = sum(num)

#percent of each subtype in the MS 
percent = []
for i in range(len(starmag)):
    percent.append(len(starmag[i])/numsum)

print(all_mass)  # total mass of all stars in MS
print(numsum) # number of stars in MS

#numsum = sum(num)
#print(tabulate(mydata, headers=head, tablefmt="grid"))
print(len(astar_m), len(fstar_m), len(gstar_m), len(kstar_m), len(mstar_m))  # num of stars of each sub type



# import pretty table
from prettytable import PrettyTable
 
# Specify the Column Names while initializing the Table
myTable = PrettyTable(["Subtype", "#", "Percentage", "Total Mass (solar units)","Percent Mass"])
 
# Add rows       type name, amt, percent of each type, total mass of each subtype, percent mass of each subtype
myTable.add_row(["A type", num[0], percent[0], total_mass[0], pmas[0]])
myTable.add_row(["F type", num[1], percent[1], total_mass[1], pmas[1]])
myTable.add_row(["G type", num[2], percent[2], total_mass[2], pmas[2]])
myTable.add_row(["K type", num[3], percent[3], total_mass[3], pmas[3]])
myTable.add_row(["M type", num[4], percent[4], total_mass[4], pmas[4]])
 
print(myTable)


'''
from here count the number of WD and giants. estimate WD to be ~1 solar mass and giants to be ~5 solar masses
'''
# WD 2,070
# Giant 100

# White Dwarf
num_wd = 2070
mass_wd = num_wd

# Red Giants
num_giant = 100
mass_giant = num_giant * 5

tot_mass = all_mass + mass_wd + mass_giant     # total mass of (MS + giants + WD) in solar mass

# solar mass = 2*10**30 kg
volume = (4/3)*(3.14)*(50**3)
print("Volume of a 50pc sphere: ", volume, "cubic parsec")

#density in solarmass per cubic parsec
density = (all_mass + mass_wd + mass_giant) /volume
print("density of stars in a 50pc sphere by mass: ", density, "solar mass per cubic parsec")

#density in solarmass per cubic parsec
density_num = (numsum + num_wd + num_giant) /volume
print("density of stars in a 50pc sphere: ", density_num, "# of stars per cubic parsec")

#in kilogram
solar_to_kg = (2*(10**30))
density_kg = density*solar_to_kg
print("Density of stars in a 50pc sphere in kg: ", density_kg, "kg per cubic parsec")


#volume of our galaxy
volume_gal = (200)*(3.14)*(15000**2)
print("Volume of our Galaxy assuming cylindrical shape: ", volume_gal,"cubic parsec")

# mass of milky way in solar masses
mass_gal_solar = density*volume_gal
print("Mass of Galaxy: ",mass_of_gal, "Solar Masses")

# mass of milky wasy in kilograms
mass_gal_kg = solar_to_kg*mass_gal_solar
print("Mass of Galaxy: ", mass_gal_kg, "kilograms")

