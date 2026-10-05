import os
from ase.lattice.hexagonal import Graphene
from ase.calculators.lammpsrun import LAMMPS
from ase.optimize import BFGS
from ase.io import write
from ase.visualize import view

if 'ASE_LAMMPSRUN_COMMAND' not in os.environ:
    os.environ['ASE_LAMMPSRUN_COMMAND'] = 'lmp' 

print("--- ÉTAPE 1 : Création du support de graphène 32x32 ---")

a = 2.46  
c = 150.0  # Vide fortement augmenté à 150 Å pour éviter les artefacts périodiques en Z
atoms = Graphene(symbol='C', latticeconstant={'a': a, 'c': c}, size=(32, 32, 1))
atoms.center(axis=2)

print(f"Supercellule créée avec succès : {len(atoms)} atomes de carbone.")

parameters = {
    'pair_style': 'tersoff',
    'pair_coeff': ['* * tersoff.pot C'],
    'units': 'metal',
    'boundary': 'p p f'  # Ouvert (shrink-wrapped) en Z pour empêcher le franchissement
}

atoms.calc = LAMMPS(
    files=['tersoff.pot'],
    **parameters
)

print("Lancement de la relaxation de géométrie (BFGS)...")
opt = BFGS(atoms, trajectory='graphene_32x32_relaxed.traj')
opt.run(fmax=0.01)  

output_filename = "graphene_32x32_relaxed.xyz"
write(output_filename, atoms)
print(f"Relaxation terminée ! Structure enregistrée sous '{output_filename}'.")