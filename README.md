# ORCA + Multiwfn Automation Workflow

This software automates a workflow to:

1. Input a molecular coordinate file (`.xyz`)  
2. Use **ORCA** to optimize the geometry and generate wavefunction (`.gbw`) files  
3. Generate wavefunctions for adjacent positively and negatively charged ions  
4. Convert ORCA `.gbw` files to `.wfn` files  
5. Use **Multiwfn** to calculate the Fukui function of the molecule  
6. Export Gaussian cube files for visualization in software like **VMD** or **Avogadro 2**  

With this workflow, users can:

- Visualize nucleophilic and electrophilic regions of any molecule  
- Easily view HOMO, LUMO, and other orbitals without manually parsing ORCA output  
- Perform additional orbital analyses efficiently  

---

## Dependencies

- **ORCA**  
  - Download via the [official ORCA forum](https://orcaforum.kofo.mpg.de/)  
  - Academic users must register and accept the ORCA license  
  - Commercial use requires a commercial ORCA license  
  - **Note:** This repository does **not** include the ORCA executable  

- **Multiwfn**  
  - Free and open-source under its license  
  - Download from the [official website](http://sobereva.com/multiwfn/)  
  - Users must agree to the Multiwfn license terms  

---

## Usage Notes

- Ensure ORCA and Multiwfn executables are correctly installed and accessible from your system PATH  
- The workflow handles generating wavefunctions, calculating Fukui functions, and producing visualizable cube files automatically  

---

## Citation

If you publish results obtained using this workflow, please cite **Multiwfn**:

1. Tian Lu & Feiwu Chen, *J. Comput. Chem.* **33**, 580–592 (2012)  
2. Tian Lu, *J. Chem. Phys.* **161**, 082503 (2024)  

