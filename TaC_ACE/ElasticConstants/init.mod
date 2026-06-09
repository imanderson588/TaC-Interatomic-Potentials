# NOTE: This script can be modified for different atomic structures, 
# units, etc. See in.elastic for more info.
#

# Define the finite deformation size. Try several values of this
# variable to verify that results do not depend on it.
variable up equal 1.0e-6
 
# Define the amount of random jiggle for atoms
# This prevents atoms from staying on saddle points
variable atomjiggle equal 1.0e-5

# Uncomment one of these blocks, depending on what units
# you are using in LAMMPS and for output

# metal units, elastic constants in eV/A^3
#units		metal
#variable cfac equal 6.2414e-7
#variable cunits string eV/A^3

# metal units, elastic constants in GPa
units		metal
variable cfac equal 1.0e-4
variable cunits string GPa

# real units, elastic constants in GPa
#units		real
#variable cfac equal 1.01325e-4
#variable cunits string GPa

# Define minimization parameters
variable etol equal 1e-3
variable ftol equal 1.0e-25
variable maxiter equal 1000
variable maxeval equal 1000
variable dmax equal 1.0e-2

# generate the box and atom positions using a diamond lattice
variable a   equal 2
boundary p p p


lattice         fcc $a 
region          box prism -4 4 -4 4 -4 4 0 0 0
create_box      2 box
create_atoms    1 box
lattice         fcc $a origin 0.5 0.5 0.5
create_atoms    2 box


mass 1 180.9479  # Atomic mass of Ta
mass 2 12.011    # Atomic mass of C
