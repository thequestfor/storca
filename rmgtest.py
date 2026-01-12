# rmg_hoooh.py
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import CoreEdgeReactionModel

# --- 1. Define the molecule ---
# HOOOOH: simplified as SMILES for RMG
hooooh = Species().fromSMILES('OOOO')  # 4 O chain, implicit H

# --- 2. Create an RMG instance ---
rmg = RMG()
rmg.species = [hooooh]
rmg.temperature = (298,)  # K
rmg.pressure = (1e5,)     # Pa
rmg.time = (0, 1)         # seconds
rmg.maxCoreSize = 10      # max species in network
rmg.verbose = 3

# --- 3. Build the reaction network ---
print("Generating possible reactions for HOOOOH...")
rmg.generate()

# --- 4. Extract reactions involving HOOOOH ---
hoooh_reactions = []
for rxn in rmg.reactions:
    if hooooh in rxn.reactants:
        hooooh_reactions.append(rxn)

# --- 5. Sort reactions by estimated rate at 298 K ---
def rate_at_298(rxn):
    try:
        return rxn.estimateRateCoefficient(T=298.0)
    except Exception:
        return 0.0  # if no rate available

hooooh_reactions.sort(key=rate_at_298, reverse=True)

# --- 6. Print results ---
print("\nHOOOOH decomposition pathways sorted by estimated rate constant (s^-1):\n")
for i, rxn in enumerate(hoooh_reactions, 1):
    k = rate_at_298(rxn)
    reactants = ' + '.join([sp.label for sp in rxn.reactants])
    products = ' + '.join([sp.label for sp in rxn.products])
    print(f"{i:02d}: {reactants} -> {products} | k(298 K) ≈ {k:.3e} s^-1")

