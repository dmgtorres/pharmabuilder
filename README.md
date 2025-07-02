# PharmaBuilder

## Description

PharmaBuilder is an interactive molecular pharmacophore modeling application. This program enables users to load small-molecule structures, compute 3D pharmacophore features, and visualise them in a user-friendly GUI. Key functionalities include dynamic feature filtering, interactive 3D rendering, and export of pharmacophore data and snapshots.

Further improvements for this program are under development.

## Setup

### Requirements

- Operating system: Windows (version 10 or later), macOS (version 10.15 or later), or Linux
- Python (version 3.8 or later)
- Python module dependencies:
  - PyQt5
  - PyQtWebEngine
  - RDKit
  - py3Dmol
  - Pillow
  - JSON (built-in)

### Installation

1. **Clone the repository**

   ```bash
   git clone https://github.com/yourusername/PharmaBuilder.git
   cd PharmaBuilder
   ```

2. **Create and activate a virtual environment**

   ```bash
   python -m venv venv
   source venv/bin/activate     # macOS/Linux
   venv\\Scripts\\activate.bat  # Windows
   ```

3. **Install dependencies**

   ```bash
   pip install -r requirements.txt
   ```

4. **Run the application**

   ```bash
   python src/main.py
   ```

## Tutorial

### Input

PharmaBuilder supports the following input formats:
- **Manual input**: Copy and paste the SMILES string of the molecule directly into the text field to parse it using RDKit
- **File upload**: Upload the MDL Molfile (.mol), SDF file (.sdf), or PDB file (.pdb) of the molecule to extract atom positions and bonds

### Output

PharmaBuilder provides the following export options:
- **Filtered pharmacophore export**: Download a JSON file (.json) or plain-text file (.txt) containing pharmacophore features (type, atom identifier, and centroid 3D coordinates)
- **Snapshot export**: Export a PNG image (.png) of the 2D molecular structure with pharmacophore features highlighted and a legend of feature types

### Features

- Automatic hydrogen addition, 3D embedding, and force-field optimisation
- Molecular structure 3D visualisation:
  - Wire model
  - Stick model
  - Ball and stick model
  - Sphere model
  - Cartoon model
- Pharmacophore feature detection, filtering, and 3D visualisation:
  - Hydrogen bond donors (HBD)
  - Hydrogen bond acceptors (HBA)
  - Positive ionisable groups (PI)
  - Negative ionisable groups (NI)
  - Hydrophobic groups (H)
  - Aromatic centers (AR)
- Interactive legend to toggle pharmacophore feature visibility

### Utilities

- **Clear filters**: Reset all active feature filters
- **Snapshot**: Save the current 3D display with pharmacophore feature legend (with or without background) as `pharmacophore_snapshot.png`
- **Export features**: Save the current filtered features to `filtered_features.json`
- **Export centroids**: Save centroid coordinates to `centroids.txt`

## Bugs

Please send any bug reports to the authors with a complete test case.

## References

1. Vuorinen, A.; Schuster, D. *Methods* **2015**, *71*, 113–134. https://doi.org/10.1016/j.ymeth.2014.10.013
1. Giordano, D. *et al.* *Pharmaceuticals* **2022**, *15*(5), 646. https://doi.org/10.3390/ph15050646

## Acknowledgements

This work was developed under the supervision of Prof. João Aires de Sousa (NOVA University Lisbon) and Prof. Gilles Marcou (University of Strasbourg) for the Molecular Modelling Project, as part of the Erasmus Mundus Joint Master in ChEMoinformatics+, co-funded by the European Union (Programme ERASMUS2027, ERASMUS-EDU-2021-PEX-EMJM-MOB, Project number 101050809).

## Changelog

Updates to this document will be described here.

### 0.1.0 (2025-06-26)

Initial version.

## Credits and licensing

Copyright (c) 2025 André Maio, Diogo Torres, Pedro Lourenço (NOVA University Lisbon). All rights reserved.

This program is licensed under the GNU General Public License (version 3.0).

**Authors**

- André Maio <<a.maio@campus.fct.unl.pt>>
- Diogo Torres <<dm.torres@campus.fct.unl.pt>>
- Pedro Lourenço <<pg.lourenco@campus.fct.unl.pt>>
