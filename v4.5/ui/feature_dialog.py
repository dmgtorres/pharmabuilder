import json
import io
from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, QTextEdit,
                             QLineEdit, QPushButton, QCheckBox, QListWidget, 
                             QGroupBox, QGridLayout, QFileDialog
)
from PIL import ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont

from core.mol_to_features import get_structured_features

class FeatureDialog(QDialog):
    def __init__(self, mol, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Pharmacophore Features")
        self.setMinimumWidth(600)

        self.features = get_structured_features(mol)

        layout = QVBoxLayout()

        top_layout = QHBoxLayout()
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("Search group name...")
        self.export_btn = QPushButton("Export to JSON")
        self.export_txt_btn = QPushButton("Export Centroids (.txt)")
        self.export_img_btn = QPushButton("Export 2D Snapshot (.png)")
        top_layout.addWidget(QLabel("Search:"))
        top_layout.addWidget(self.search_input)
        top_layout.addWidget(self.export_btn)
        top_layout.addWidget(self.export_txt_btn)
        top_layout.addWidget(self.export_img_btn)
        layout.addLayout(top_layout)

        options_box = QGroupBox("Snapshot Export Options")
        options_layout = QVBoxLayout()
        self.transparent_bg_chk = QCheckBox("Non-Transparent Background")
        options_layout.addWidget(self.transparent_bg_chk)
        options_box.setLayout(options_layout)
        layout.addWidget(options_box)

        self.filter_checkboxes = {}
        feature_types = {"HBD", "HBA", "AR", "H"}
        filter_box = QGroupBox("Filter by Pharmacophore Feature Type")
        filter_layout = QGridLayout()
        for i, t in enumerate(sorted(feature_types)):
            cb = QCheckBox(t)
            cb.setChecked(True)
            self.filter_checkboxes[t] = cb
            filter_layout.addWidget(cb, i // 2, i % 2)
        filter_box.setLayout(filter_layout)
        layout.addWidget(filter_box)

        self.feature_list = QListWidget()
        self.feature_list.currentRowChanged.connect(self.show_feature_details)
        layout.addWidget(self.feature_list)

        self.textbox = QTextEdit()
        self.textbox.setReadOnly(True)
        layout.addWidget(self.textbox)

        self.setLayout(layout)

        self.search_input.textChanged.connect(self.apply_filters)
        for cb in self.filter_checkboxes.values():
            cb.stateChanged.connect(self.apply_filters)
        self.export_btn.clicked.connect(self.export_to_json)

        self.export_txt_btn.clicked.connect(self.export_to_txt)
        self.export_img_btn.clicked.connect(self.export_to_png)

        self.apply_filters()

    def apply_filters(self):
        text = self.search_input.text().lower()
        active = [k for k, cb in self.filter_checkboxes.items() if cb.isChecked()]
        self.filtered_features = []
        self.feature_list.clear()

        for feat in self.features:
            if text and text not in feat["functional_group"].lower():
                continue
            if not any(f in active for f in feat["pharmacophore_features"]):
                continue
            self.filtered_features.append(feat)
            group_name = feat.get("functional_group", "Unknown")
            centroid = feat.get("centroid", ())
            label = f"{group_name} at {centroid}" if centroid else group_name
            self.feature_list.addItem(label)

        self.feature_list.setCurrentRow(0)
        self.show_feature_details(0)

    def show_feature_details(self, index):
        if 0 <= index < len(self.filtered_features):
            feat = self.filtered_features[index]
            self.textbox.setText(json.dumps(feat, indent=2))
        else:
            self.textbox.clear()

    def display_features(self):
        output = ""
        for i, feat in enumerate(self.filtered_features, 1):
            output += f"Feature {i}:"
            output += f"  Group: {feat['functional_group']}"
            output += f"  Types: {', '.join(feat['pharmacophore_features'])}"
            output += f"  Atom Indices: {feat['atom_indices']}"
            output += f"  Coordinates:"
            for coord in feat['atom_coordinates']:
                output += f"    {coord}"
            if 'centroid' in feat:
                output += f"  Centroid: {feat['centroid']}"
            output += ""
        self.textbox.setText(output)

    def export_to_json(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save Features", "filtered_features.json", "JSON Files (*.json)")
        if path:
            with open(path, 'w') as f:
                json.dump(self.filtered_features, f, indent=2)

    def export_to_txt(self):
        if not self.filtered_features:
            return

        path, _ = QFileDialog.getSaveFileName(self, "Save Centroid Data", "centroids.txt", "Text Files (*.txt)")
        if path:
            with open(path, "w") as f:
                f.write("Pharmacophore Centroid Coordinates\n")
                f.write("==================================\n\n")
                for i, feat in enumerate(self.filtered_features, 1):
                    if "centroid" not in feat:
                        continue
                    group = feat.get("functional_group", "Unknown")
                    types = ", ".join(feat.get("pharmacophore_features", []))
                    x, y, z = feat["centroid"]
                    f.write(f"Feature {i}:\n")
                    f.write(f"  Types: {types} {x:.3f} {y:.3f} {z:.3f}\n\n")

    def export_to_png(self):
        mol = self.parent().molecule_store.get_rdkit_mol()
        if mol is None:
            return

        highlight_atoms = set()
        atom_colors = {}
        present_tags = set()

        COLOR_MAP = {
            "HBD": (230/255, 50/255, 50/255),
            "HBA": (50/255, 100/255, 230/255),
            "AR": (30/255, 150/255, 30/255),
            "H": (120/255, 120/255, 120/255),
            "PI": (255/255, 127/255, 20/255),
            "NI": (153/255, 51/255, 204/255),
            "ZnB": (255/255, 215/255, 0/255),
            "Unknown": (0.2, 0.2, 0.2)
        }

        features = self.filtered_features
        for feat in features:
            tag = feat.get("pharmacophore_features", ["Unknown"])[0]
            present_tags.add(tag)
            for idx in feat["atom_ids"]:
                highlight_atoms.add(idx)
                atom_colors[idx] = COLOR_MAP.get(tag, (0.2, 0.2, 0.2))

        mol = Chem.RemoveHs(mol)
        AllChem.Compute2DCoords(mol)
        rdMolDraw2D.PrepareMolForDrawing(mol)

        drawer = rdMolDraw2D.MolDraw2DCairo(600, 600)
        drawer.drawOptions().legendFontSize = 16
        drawer.drawOptions().clearBackground = self.transparent_bg_chk.isChecked()

        drawer.DrawMolecule(
            mol,
            highlightAtoms=list(highlight_atoms),
            highlightAtomColors=atom_colors,
            highlightAtomRadii={idx: 0.45 for idx in highlight_atoms},
            confId=-1
        )
        drawer.FinishDrawing()
        img_bytes = drawer.GetDrawingText()

        # Convert RDKit image to PIL
        img = Image.open(io.BytesIO(img_bytes))

        # Draw legend
        draw = ImageDraw.Draw(img)
        font = ImageFont.load_default()

        legend_x = 420
        legend_y = 20
        draw.text((legend_x, legend_y), "Legend:", fill="black", font=font)
        legend_y += 15

        for tag in sorted(present_tags):
            color = tuple(int(c * 255) for c in COLOR_MAP.get(tag, (50, 50, 50)))
            draw.rectangle([legend_x, legend_y, legend_x + 12, legend_y + 12], fill=color)
            draw.text((legend_x + 18, legend_y), tag, fill="black", font=font)
            legend_y += 16

        # Save
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Snapshot",
            "pharmacophore_snapshot.png",
            "PNG Files (*.png)"
        )

        if path:
            img.save(path)
