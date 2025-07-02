# Molecular Modelling Project

## Features

### 2D/3D molecular visualisation

#### Overview

An interactive molecular visualisation and modeling environment for manipulation and viewing of molecules in 2D/3D representations

#### Development steps

- Integration with open-source software (e.g., [Avogadro](https://avogadro.cc/), [MolView](https://molview.org/))
- Integration in GUI with PyQt5

### Ligand-guided pharmacophore generation

#### Overview

A webserver and/or application for ligand- and receptor-based drug design strategies based on a ligand-guided method for generating pharmacophore hypotheses by using ligand–receptor interactions

#### Development steps

- Ligand to pharmacophore hypotheses
  - From ligand in MOL format
  - From ligand in ligand–receptor complex crystal structure in PDB format
- Generation of molecules with molecular docking features based on pharmacophore hypotheses using [PGMG Server](https://www.csuligroup.com/PGMG)
- Molecular docking tests of generated molecules using [AutoDock Vina](https://vina.scripps.edu/)
- Ranking of generated molecules based on molecular docking affinity
- Ranking of pharmacophore hypotheses based on ranking of generated molecules
- Integration in GUI with PyQt5

## References

### Articles

- [*J. Chem. Inf. Model.* **2024**, *64*(10), 4263–4276](https://pubs.acs.org/doi/10.1021/acs.jcim.3c01920)
- [*J. Cheminform.* **2023**, *15*, 117](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-023-00786-w)
- [*J. Cheminform.* **2014**, *6*, 14](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-6-14)

### Tutorials

- [Ligand-based pharmacophores](https://projects.volkamerlab.org/teachopencadd/talktorials/T009_compound_ensemble_pharmacophores.html)
- [Lists of molecules in RDKit](https://asteeves.github.io/blog/2015/01/12/lists-in-rdkit/)

### Resources

- [AutoDock Vina](https://vina.scripps.edu/)
- [Avogadro](https://avogadro.cc/)
- [BIOVIA Discovery Studio Visualizer](https://discover.3ds.com/discovery-studio-visualizer-download)
- [LigandScout](http://www.inteligand.com/ligandscout/)
- [Maestro](https://www.schrodinger.com/platform/products/maestro/)
- [Molecular Operating Environment](https://www.chemcomp.com/en/Products.htm)
- [MolView](https://molview.org/)
- [PDBe CCDUtils](https://github.com/PDBeurope/ccdutils)
- [Pharmacophore-Guided Molecular Generation Server](https://www.csuligroup.com/PGMG)
- [PharmaCore](https://computorgchemunisa.org/pharmacore/)
- [Phase](https://www.schrodinger.com/platform/products/phase/)
