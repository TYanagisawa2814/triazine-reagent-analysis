
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import rdFreeSASA


TOTAL_CHARGE = 0 # or (+1)

# =========================================
# Bondi van der Waals radii (Å)
# Bondi, J. Phys. Chem. 1964, 68, 441–451
# DOI: 10.1021/j100785a001
# =========================================
VDW_RADII = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52}

xyz_file = r"xtbopt.xyz"
mol = Chem.MolFromXYZFile(xyz_file)

rdDetermineBonds.DetermineBonds(mol, charge=TOTAL_CHARGE)
print("Detected structure:")
print(Chem.MolToSmiles(mol))

radii = []
for atom in mol.GetAtoms():
    symbol = atom.GetSymbol()

    if symbol not in VDW_RADII:
        raise ValueError(f"No Bondi van der Waals radius defined for {symbol}")

    radii.append(VDW_RADII[symbol])

opts = rdFreeSASA.SASAOpts()
opts.algorithm = rdFreeSASA.SASAAlgorithm.ShrakeRupley
opts.probeRadius = 1.4

# Total SASA
total_sasa = rdFreeSASA.CalcSASA(mol, radii, confIdx=0, opts=opts)

# Polar / apolar SASA
# Polar:
# N, O, and H bonded to these atoms
# Apolar:
# Remaining C and H atoms
polar_sasa = 0.0
apolar_sasa = 0.0

for atom in mol.GetAtoms():
    symbol = atom.GetSymbol()
    sasa = float(atom.GetProp("SASA"))

    is_polar = False

    if symbol in ["N", "O"]:
        is_polar = True

    elif symbol == "H":
        neighbors = atom.GetNeighbors()

        if len(neighbors) == 1:
            bonded_symbol = neighbors[0].GetSymbol()

            if bonded_symbol in ["N", "O"]:
                is_polar = True

    if is_polar:
        polar_sasa += sasa
    else:
        apolar_sasa += sasa

polar_fraction = polar_sasa / total_sasa
apolar_fraction = apolar_sasa / total_sasa


# Output
print("----------------------------------------")
print("SASA RESULTS")
print("")

print(f"Total SASA        : {total_sasa:.3f} Å²")
print(f"Polar SASA        : {polar_sasa:.3f} Å²")
print(f"Apolar SASA       : {apolar_sasa:.3f} Å²")
print(f"Polar fraction    : {polar_fraction:.4f}")
print(f"Apolar fraction   : {apolar_fraction:.4f}")
print("")
print("----------------------------------------")
