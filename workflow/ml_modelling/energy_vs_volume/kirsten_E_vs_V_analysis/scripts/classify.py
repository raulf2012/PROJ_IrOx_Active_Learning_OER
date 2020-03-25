import ase
from protosearch.build_bulk.cell_parameters import CellParameters
from protosearch.build_bulk.classification import PrototypeClassification
from ase.db import connect

db = connect('FinalStructures.db')


for row in db.select():
    #if row.get('p_name', None):
    #    continue

    atoms = row.toatoms()
    print(atoms.info)
    try:
        PC = PrototypeClassification(atoms, tolerance=1e-3)

        prototype = PC.get_classification(include_parameters=False)

    except:
        continue
    p_name = prototype['p_name']
    print(p_name)
    spacegroup = prototype['spacegroup']    
    
    db.update(row.id, p_name=p_name, spacegroup=spacegroup)
