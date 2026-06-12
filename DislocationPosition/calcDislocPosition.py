from ovito.io import import_file, export_file
from ovito.modifiers import SelectTypeModifier, CommonNeighborAnalysisModifier, DeleteSelectedModifier, ExpressionSelectionModifier
import warnings
warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')
import ovito._extensions.pyscript


"""
This file uses the Ovito python API to locate and output the position of an edge dislocation within the simulation cell.
This is accomplished by first deleting all carbon atoms, then delecting all atoms in the FCC structure, then selecting atoms
in the center of the simulation cell and calculating the average y-position of the remaining atoms.
"""

pipeline = import_file("dump.652")


pipeline.modifiers.append(SelectTypeModifier(operate_on="particles", property="Particle Type", types = {"Type 2"} ))

pipeline.modifiers.append(DeleteSelectedModifier())

pipeline.modifiers.append(CommonNeighborAnalysisModifier())

pipeline.modifiers.append

data = pipeline.compute()

z_coordinates = data.particles['Position'][:,2]

z_min = min(z_coordinates)
z_max = max(z_coordinates)

surface_cutoff = 5

expression = f"Position.Z < {z_min + surface_cutoff} || Position.Z > {z_max - surface_cutoff}"
pipeline.modifiers.append(ExpressionSelectionModifier(expression=expression))

pipeline.modifiers.append(DeleteSelectedModifier())

final_data = pipeline.compute()


export_file(pipeline, "output_cleaned.dump", format="lammps/dump", columns = ["Position.X", "Position.Y", "Position.Z"])
