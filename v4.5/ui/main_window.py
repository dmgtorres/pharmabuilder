import py3Dmol
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLineEdit, 
                            QLabel, QFileDialog, QComboBox, QCheckBox, QMessageBox, 
                            QGroupBox, QScrollArea, QHBoxLayout, QVBoxLayout)
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont

from core.molecule_store import MoleculeStore
from ui.feature_dialog import FeatureDialog
from core.mol_to_features import get_structured_features

CENTROID_COLORS = {
            "HBD": "red",                   # Hydrogen Bond Donor
            "HBA": "blue",                  # Hydrogen Bond Acceptor
            "AR": "green",                  # Aromatic
            "H": "gray",                    # Hydrophobe
            "PI": "orange",                 # Positively Ionizable
            "NI": "purple",                 # Negatively Ionizable
            "ZnB": "gold",                  # Zinc Binder
            "Unknown": "black"              # Fallback for untagged types
}

COLOR_MAP = {
            "HBD": "#e63232",            # red
            "HBA": "#3264e6",            # blue
            "AR": "#1e961e",             # green
            "H":   "#787878",            # gray
            "PI": "#ff7f14",             # orange
            "NI": "#9933cc",             # purple
            "ZnB": "#ffd700",            # gold
            "Unknown": "#333333"         # black
        }

class MoleculeAppWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Pharmacophore Molecule Viewer")
        self.setMinimumSize(900, 700)
        self.setWindowState(Qt.WindowMaximized)

        self.molecule_store = MoleculeStore()
        self.active_feature_filters = set()

        self._init_ui()

    def _init_ui(self):
        main_layout = QVBoxLayout()
        controls = QHBoxLayout()

        self.smiles_input = QLineEdit()
        self.load_smiles_btn = QPushButton("Render SMILES")
        self.load_file_btn = QPushButton("Open File")
        self.center_btn = QPushButton("Center")
        self.feature_btn = QPushButton("Pharmacophore Features")
        self.clear_filters_btn = QPushButton("Clear Filters")

        self.viewer = QWebEngineView()
        self.style_combo = QComboBox()
        self.style_combo.addItems(["Stick", "Sphere", "Line", "Ball and Stick"])

        self.label_toggle = QCheckBox("Show Atom Labels")
        self.legend_toggle = QCheckBox("Show Legend")
        self.legend_toggle.setChecked(True)

        controls.addWidget(self.load_smiles_btn)
        controls.addWidget(self.load_file_btn)
        controls.addWidget(self.center_btn)
        controls.addWidget(self.feature_btn)
        controls.addWidget(self.clear_filters_btn)
        controls.addWidget(QLabel("Style:"))
        controls.addWidget(self.style_combo)
        controls.addWidget(self.label_toggle)
        controls.addWidget(QLabel("SMILES:"))
        controls.addWidget(self.smiles_input)
        controls.addWidget(self.legend_toggle)

        main_layout.addLayout(controls)

        self.legend_layout = QVBoxLayout()
        self.legend_box = QGroupBox("Pharmacophore Legend")
        self.legend_box.setLayout(self.legend_layout)

        self.legend_scroll = QScrollArea()
        self.legend_scroll.setWidget(self.legend_box)
        self.legend_scroll.setWidgetResizable(True)
        self.legend_scroll.setFixedWidth(200)

        viewer_layout = QHBoxLayout()
        viewer_layout.addWidget(self.viewer, stretch=1)
        viewer_layout.addWidget(self.legend_scroll)
        main_layout.addLayout(viewer_layout)

        self.setLayout(main_layout)

        self.load_smiles_btn.clicked.connect(self._handle_smiles)
        self.load_file_btn.clicked.connect(self._handle_file)
        self.center_btn.clicked.connect(self._display_molecule)
        self.feature_btn.clicked.connect(self._show_pharmacophores)
        self.clear_filters_btn.clicked.connect(self._clear_filters)
        self.style_combo.currentTextChanged.connect(self._display_molecule)
        self.label_toggle.stateChanged.connect(self._display_molecule)
        self.legend_toggle.stateChanged.connect(self._toggle_legend)

    def _toggle_legend(self, state):
        self.legend_scroll.setVisible(state == Qt.Checked)

    def _handle_smiles(self):
        smiles = self.smiles_input.text().strip()
        if not smiles:
            QMessageBox.warning(self, "Missing SMILES", "Please enter a SMILES string.")
            return
        self.molecule_store.load_from_smiles(smiles)
        self._display_molecule()

    def _handle_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Open Molecule File", "", "Molecule Files (*.mol *.sdf *.pdb)")
        if file_path:
            self.molecule_store.load_from_file(file_path)
            self._display_molecule()

    def _display_molecule(self):
        mol_data = self.molecule_store.get_viewable_data()
        if not mol_data:
            return

        mol_type, mol_block, rdkit_mol = mol_data
        view = py3Dmol.view(width=self.viewer.width(), height=self.viewer.height())
        view.addModel(mol_block, mol_type)

        style = self.style_combo.currentText().lower()
        if style == "ball and stick":
            view.setStyle({'stick': {"linewidth": 2.0}, 'sphere': {'scale': 0.3}})
        else:
            view.setStyle({style: {"linewidth": 2.0}})

        if self.label_toggle.isChecked() and rdkit_mol:
            conf = rdkit_mol.GetConformer()
            for atom in rdkit_mol.GetAtoms():
                idx = atom.GetIdx()
                pos = conf.GetAtomPosition(idx)
                view.addLabel(f"{atom.GetSymbol()}{idx}", {
                    "position": {"x": pos.x, "y": pos.y, "z": pos.z},
                    "backgroundColor": "rgba(255,255,255,0.8)",
                    "fontColor": "black",
                    "fontSize": 11,
                    "borderThickness": 0.5,
                    "inFront": True
                })

        CENTROID_COLORS = {
            "HBD": "red", "HBA": "blue", "AR": "green", "H": "gray",
            "PI": "orange", "NI": "purple", "ZnB": "gold", "Unknown": "black"
        }
        feats = get_structured_features(rdkit_mol)

        self._update_legend(feats)

        for feat in feats:
            types = feat.get("pharmacophore_features", [])
            if self.active_feature_filters and not any(t in self.active_feature_filters for t in types):
                continue
            if 'centroid' in feat:
                centroid = feat['centroid']
                color = next((CENTROID_COLORS[t] for t in types if t in CENTROID_COLORS), "gray")
                view.addSphere({
                    "center": {"x": centroid[0], "y": centroid[1], "z": centroid[2]},
                    "radius": 1.0,
                    "color": color,
                    "opacity": 0.6,
                    "wireframe": True
                })

        view.zoomTo()
        self.viewer.setHtml(view._make_html())

    def _show_pharmacophores(self):
        mol = self.molecule_store.get_rdkit_mol()
        if mol is None:
            QMessageBox.warning(self, "No Molecule", "Load a molecule before requesting features.")
            return

        dialog = FeatureDialog(mol, parent=self)
        dialog.exec_()

    def _update_legend(self, features):
        while self.legend_layout.count():
            child = self.legend_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

        tag_set = set()
        for feat in features:
            tag_set.update(feat.get("pharmacophore_features", []))

        for tag in sorted(tag_set):
            row = QHBoxLayout()
            color_box = QLabel()
            color_box.setFixedSize(14, 14)
            color_box.setStyleSheet(f"background-color: {COLOR_MAP.get(tag, '#333')}; border: 1px solid black;")

            clickable_label = QLabel(f"{tag}")
            clickable_label.setCursor(Qt.PointingHandCursor)
            font = QFont()
            font.setUnderline(tag in self.active_feature_filters)
            clickable_label.setFont(font)
            clickable_label.mousePressEvent = lambda _, t=tag: self._on_legend_clicked(t)

            row.addWidget(color_box)
            row.addWidget(clickable_label)
            container = QWidget()
            container.setLayout(row)
            self.legend_layout.addWidget(container)

    def _on_legend_clicked(self, tag):
        if tag in self.active_feature_filters:
            self.active_feature_filters.remove(tag)
        else:
            self.active_feature_filters.add(tag)
        self._display_molecule()

    def _clear_filters(self):
        self.active_feature_filters.clear()
        self._display_molecule()
