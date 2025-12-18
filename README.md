# ORCA + Multiwfn Automation Workflow

This software automates a workflow for molecular electronic structure calculations using **ORCA**, with support for organic and organometallic compounds.  

It can:

1. Input a molecular coordinate file (`.xyz`) or a **SMILES string**
2. Use **ORCA** to optimize geometry and generate wavefunction (`.gbw`) files
3. Perform single-point calculations for neutral, cationic, and anionic states
4. Convert ORCA `.gbw` files to `.wfn` files
5. Parse HOMO, LUMO, and orbital energies automatically
6. Generate practical hazard estimates and chemical descriptors (SMILES, molecular weight, TPSA, XLogP)

With this workflow, users can:

- Quickly obtain optimized structures and wavefunctions
- Examine HOMO/LUMO energies without manually parsing ORCA output
- Generate WFN files for further analysis in programs like **Multiwfn** (optional)

---

## Dependencies

- **Python 3.10+**  
  - Requires standard scientific Python libraries (RDKit, PubChemPy, typer, rich)
- **ORCA**  
  - Download via the [official ORCA forum](https://orcaforum.kofo.mpg.de/)  
  - Academic users must register and accept the ORCA license  
  - Commercial use requires a commercial ORCA license  
  - **Note:** This repository does **not** include the ORCA executable  

> ⚠️ **Optional:** If you wish to analyze WFN files for Fukui functions or other electronic properties, you can use **Multiwfn** separately, but it is not required for the core workflow.

---

## Usage Notes

- Ensure ORCA is installed and accessible from your system PATH
- You can provide either:
  - A `.xyz` coordinate file, or
  - A SMILES string (the program will generate the XYZ file automatically)
- The workflow produces:
  - Optimized XYZ structures
  - Single-point calculation GBW files for neutral, +1, and -1 charge states
  - Converted WFN files (`N.wfn`, `N+1.wfn`, `N-1.wfn`)
  - Orbital energies and HOMO/LUMO information
  - Estimated practical hazards for lab handling

---

## Citation

If you publish results obtained using this workflow, please cite **ORCA**:

- F. Neese, *WIREs Comput. Mol. Sci.* **8**, e1327 (2018)

> Optional: If Multiwfn is used for post-processing:
> - Tian Lu & Feiwu Chen, *J. Comput. Chem.* **33**, 580–592 (2012)

