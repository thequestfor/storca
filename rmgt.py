from rmgpy.species import Species
from rmgpy.thermo.thermoengine import ThermoEngine
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.reaction import Reaction

# -------------------------------
# Initialize the RMG object
# -------------------------------
rmg = RMG()
rmg.database.load(
    thermoLibraries=['primaryThermoLibrary', 'BurkeH2O2'],
    kineticsFamilies='default',
    kineticsDepositories=['training'],
)

# -------------------------------
# Define HOOOH species
# -------------------------------
hoooH = Species().fromSMILES('OOO[OH]')

# Make it reactive
hoooH.generate_resonance_structures()
hoooH.reactive = True

# -------------------------------
# Set up a simple model
# -------------------------------
model = CoreEdgeReactionModel()
model.core.species.append(hoooH)

# -------------------------------
# Use RMG's ThermoEngine to get thermochemistry
# -------------------------------
thermo_engine = ThermoEngine(database=rmg.database)
thermo_engine.apply_thermo(hoooH)

# -------------------------------
# Estimate decomposition reactions
# -------------------------------
# RMG uses kinetics families to enumerate reactions
from rmgpy.data.kinetics.family import KineticsFamily

family = rmg.database.kinetics.families['R_Recombination']  # example; peroxide reactions
reactions = family.generate_reactions([hoooH], rmg.database.kinetics, resonance=True)

print("Generated Reactions for HOOOH:")
for rxn in reactions:
    print(rxn)
    print("Arrhenius params:", rxn.kinetics)

# -------------------------------
# Estimate half-life (simplest)
# -------------------------------
import math

T = 298  # K
R = 8.314  # J/mol·K

# Pick the fastest pathway (lowest Ea)
k_values = []
for rxn in reactions:
    if hasattr(rxn.kinetics, 'T0') and hasattr(rxn.kinetics, 'Ea'):
        A = rxn.kinetics.A.value_si
        n = rxn.kinetics.n
        Ea = rxn.kinetics.Ea.value_si
        k = A * (T / rxn.kinetics.T0.value_si) ** n * math.exp(-Ea / (R * T))
        k_values.append(k)

if k_values:
    k_fastest = max(k_values)
    t_half = math.log(2) / k_fastest
    print(f"Estimated half-life at 298 K: {t_half:.2f} s")
else:
    print("No reactions with kinetics found.")

