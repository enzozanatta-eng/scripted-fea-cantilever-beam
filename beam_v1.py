"""
Cantilever beam: scripted FEA with Gmsh + CalculiX + Python
Step 1: build geometry and hex mesh with the Gmsh API, write the CalculiX
input file, run the solver from the command line, read the tip deflection
and compare it with Euler-Bernoulli beam theory.

Units: mm, N, MPa
"""
import os
import subprocess
import gmsh

# ---------------- User settings ----------------
CCX = r"C:\Users\enzoz\Downloads\PrePoMax\PrePoMax v3.0.0\Solver\ccx_dynamic.exe"

L, B, H = 200.0, 20.0, 10.0      # length (x), width (y), height (z) [mm]
P = 100.0                        # tip load in -z [N]
E, NU = 210000.0, 0.3            # steel [MPa]
NX, NY, NZ = 40, 4, 4            # elements in each direction
JOB = "beam"
# ------------------------------------------------


def build_mesh(nx, ny, nz):
    """Structured hexahedral mesh by extrusion. Returns nodes and elements."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("beam")

    p = gmsh.model.geo.addPoint(0, 0, 0)
    line = gmsh.model.geo.extrude([(0, p)], L, 0, 0, [nx])
    surf = gmsh.model.geo.extrude([line[1]], 0, B, 0, [ny], recombine=True)
    gmsh.model.geo.extrude([surf[1]], 0, 0, H, [nz], recombine=True)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    tags, coords, _ = gmsh.model.mesh.getNodes()
    nodes = {int(t): (coords[3 * i], coords[3 * i + 1], coords[3 * i + 2])
             for i, t in enumerate(tags)}

    etypes, etags, enodes = gmsh.model.mesh.getElements(3)
    elements = []
    for et, tg, nd in zip(etypes, etags, enodes):
        if et != 5:  # 5 = 8-node hexahedron
            raise RuntimeError(f"Unexpected element type {et}")
        for i, t in enumerate(tg):
            elements.append((int(t), [int(n) for n in nd[8 * i:8 * i + 8]]))
    gmsh.finalize()
    return nodes, elements


def fix_orientation(nodes, conn):
    """Make sure the hexahedron has a positive Jacobian (CalculiX needs it)."""
    def sub(a, b):
        return [a[k] - b[k] for k in range(3)]

    p0, p1, p3, p4 = (nodes[conn[i]] for i in (0, 1, 3, 4))
    a, b, c = sub(p1, p0), sub(p3, p0), sub(p4, p0)
    cross = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
    vol = sum(cross[k] * c[k] for k in range(3))
    if vol < 0:
        conn = [conn[0], conn[3], conn[2], conn[1], conn[4], conn[7], conn[6], conn[5]]
    return conn


def write_inp(nodes, elements, path):
    tol = 1e-6
    fix = [n for n, (x, y, z) in nodes.items() if abs(x) < tol]
    tip = [n for n, (x, y, z) in nodes.items() if abs(x - L) < tol]

    def write_set(f, ids):
        for i in range(0, len(ids), 10):
            f.write(", ".join(str(n) for n in ids[i:i + 10]) + "\n")

    with open(path, "w") as f:
        f.write("*NODE, NSET=NALL\n")
        for n, (x, y, z) in nodes.items():
            f.write(f"{n}, {x:.6f}, {y:.6f}, {z:.6f}\n")
        # C3D8I: 8-node hex with incompatible modes, accurate in bending
        f.write("*ELEMENT, TYPE=C3D8I, ELSET=EALL\n")
        for e, conn in elements:
            conn = fix_orientation(nodes, conn)
            f.write(f"{e}, " + ", ".join(str(n) for n in conn) + "\n")
        f.write("*NSET, NSET=FIX\n")
        write_set(f, fix)
        f.write("*NSET, NSET=TIP\n")
        write_set(f, tip)
        f.write("*MATERIAL, NAME=STEEL\n*ELASTIC\n")
        f.write(f"{E}, {NU}\n")
        f.write("*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL\n")
        f.write("*STEP\n*STATIC\n")
        f.write("*BOUNDARY\nFIX, 1, 3, 0.0\n")
        f.write("*CLOAD\n")
        for n in tip:
            f.write(f"{n}, 3, {-P / len(tip):.8f}\n")
        f.write("*NODE PRINT, NSET=TIP\nU\n")
        f.write("*NODE FILE\nU\n*EL FILE\nS\n")
        f.write("*END STEP\n")
    return len(fix), len(tip)


def read_tip_deflection(dat_path):
    """Average z-displacement of the tip nodes from the .dat file."""
    uz = []
    reading = False
    with open(dat_path) as f:
        for line in f:
            if "displacements" in line:
                reading = True
                continue
            if reading:
                parts = line.split()
                if len(parts) == 4:
                    try:
                        int(parts[0])
                        uz.append(float(parts[3]))
                    except ValueError:
                        pass
    if not uz:
        raise RuntimeError("No displacements found in the .dat file")
    return sum(uz) / len(uz)


def run(nx=NX, ny=NY, nz=NZ, job=JOB):
    nodes, elements = build_mesh(nx, ny, nz)
    write_inp(nodes, elements, job + ".inp")
    result = subprocess.run([CCX, job], capture_output=True, text=True)
    if result.returncode != 0 or not os.path.exists(job + ".dat"):
        print(result.stdout[-2000:])
        raise RuntimeError("CalculiX failed")
    return read_tip_deflection(job + ".dat"), len(elements)


if __name__ == "__main__":
    I = B * H ** 3 / 12
    delta_theory = P * L ** 3 / (3 * E * I)

    delta_fea, n_el = run()
    error = (abs(delta_fea) - delta_theory) / delta_theory * 100

    print(f"Elements:                 {n_el}")
    print(f"Tip deflection (FEA):     {abs(delta_fea):.4f} mm")
    print(f"Tip deflection (theory):  {delta_theory:.4f} mm")
    print(f"Difference:               {error:+.2f} %")
