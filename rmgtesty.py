from rmgpy.species import Species
from rmgpy.molecule import Molecule
from rmgpy.rmg.main import RMG

# 1. Create the molecule
molecule = Molecule().from_smiles('COO')

# 2. Create a species
ch3ooh = Species(label='CH3OOH', molecule=[molecule])

# 3. Create the RMG object
rmg = RMG()

# 4. Add a simple reactor and initial species
rmg.add_simple_reactor(
    temperature=800.0,              # Kelvin
    pressure=1.0,                   # atm
    initial_mole_fractions={ch3ooh: 1.0},
    termination_conversion={ch3ooh: 0.99}
)

# 5. Execute RMG
rmg.execute()

