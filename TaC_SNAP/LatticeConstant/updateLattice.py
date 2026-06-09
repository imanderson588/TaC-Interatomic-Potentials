import subprocess


def extract_lammps_lat_const():
    with open('lat_const.out', 'w') as f:
        subprocess.run('mpirun -np 64 lmp -in lat_min.in',
                       stdout=f, shell=True)
    with open("lat_const.out", "r") as file:
        for line in file:
            if "Lattice constant" in line:
                return float(line[30:36])


def update_lat_constant(lat_const):
    filepath = "../ElasticConstants/init.mod"
    with open(filepath, "r") as file:
        lines = file.readlines()

    new_lines = []
    for line in lines:
        if line.strip().startswith("variable a   equal"):
            line = f"variable a   equal {lat_const}\n"
        new_lines.append(line)

    with open(filepath, "w") as file:
        file.writelines(new_lines)

    filepath = "../EdgeRelaxation/in.dislocation"
    with open(filepath, "r") as file:
        lines = file.readlines()

    new_lines = []
    for line in lines:
        if line.strip().startswith("variable a equal"):
            line = f"variable a qual {lat_const}\n"
        new_lines.append(line)

    with open(filepath, "w") as file:
        file.writelines(new_lines)


lat_const = extract_lammps_lat_const()
update_lat_constant(lat_const)
