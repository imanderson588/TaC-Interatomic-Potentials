from ovito import import_file
from ovito.modifiers import SelectTypeModifier, CommonNeighborAnalysisModifier, DeleteSelectedModifier

"""
This file uses the Ovito python API to locate and output the position of an edge dislocation within the simulation cell.
This is accomplished by first deleting all carbon atoms, then delecting all atoms in the FCC structure, then selecting atoms
in the center of the simulation cell and calculating the average y-position of the remaining atoms.
"""

pipeline = import_file("Dump1")


pipeline.modifiers.append(SelectTypeModifier(operate_on="particles", property="Particle Type", types = {"C"} ))

pipeline.modifiers.append(DeleteSelectedModifier(operate_on="particles"))

pipeline.modifiers.append(CommonNeighborAnalysisModifier())

pipeline.modifiers.append


