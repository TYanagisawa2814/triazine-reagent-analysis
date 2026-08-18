
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula


compounds = {"Dimethoxy": "COc1nc(N)nc(OC)n1","Dimorpholino": "Nc1nc(N2CCOCC2)nc(N3CCOCC3)n1",
             "Bis(diethylamino)": "Nc1nc(N(CC)CC)nc(N(CC)CC)n1"}

for name, smiles in compounds.items():
    mol = Chem.MolFromSmiles(smiles)

    canonical_smiles = Chem.MolToSmiles(mol)
    formula = CalcMolFormula(mol)
    molecular_weight = Descriptors.MolWt(mol)
    clogp = Crippen.MolLogP(mol) #using the Wildman–Crippen method

    print(f"\n{name}")
    print("-" * 50)
    print(f"SMILES             : {canonical_smiles}")
    print(f"Molecular formula  : {formula}")
    print(f"Molecular weight   : {molecular_weight:.3f}")
    print(f"cLogP              : {clogp:.3f}")

print("\n" + "-" * 50)
