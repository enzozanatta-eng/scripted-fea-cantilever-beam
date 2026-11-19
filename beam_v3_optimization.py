"""
Step 3: automated design optimization.
Sweeps the beam height H, runs CalculiX for each design, reads the maximum
bending stress and the tip deflection, checks them against beam theory, and
finds the lightest design that meets a safety factor and a deflection limit.

Units: mm, N, MPa, kg
"""
import os
import subprocess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gmsh

import beam_v1 as beam

# ---------------- Design problem ----------------
L, B = beam.L, beam.B          # length and width are fixed [mm]
P = 300.0                      # tip load [N]
E, NU = beam.E, beam.NU        # steel
RHO = 7.85e-6                  # density [kg/mm^3]
SY = 250.0                     # yield strength [MPa]
SF_MIN = 2.0                   # minimum safety factor
DEFL_MAX = L / 250             # deflection limit [mm]
HEIGHTS = [8, 10, 12, 14, 16, 18, 20]  # candidate heights [mm]
# ------------------------------------------------


def build_mesh(h, nx=40, ny=4, nz=6):  # keep nz = 6 (used in the stress extrapolation)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("beam")
    p = gmsh.model.geo.addPoint(0, 0, 0)
    line = gmsh.model.geo.extrude([(0, p)], L, 0, 0, [nx])
    surf = gmsh.model.geo.extrude([line[1]], 0, B, 0, [ny], recombine=True)
    gmsh.model.geo.extrude([surf[1]], 0, 0, h, [nz], recombine=True)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    nodes = {int(t): (coords[3 * i], coords[3 * i + 1], coords[3 * i + 2])
             for i, t in enumerate(tags)}
    _, etags, enodes = gmsh.model.mesh.getElements(3)
    elements = [(int(t), [int(n) for n in enodes[0][8 * i:8 * i + 8]])
                for i, t in enumerate(etags[0])]
    gmsh.finalize()
    return nodes, elements


def write_inp(nodes, elements, h, path):
    """Clamped end that lets the section contract laterally (no Poisson restraint):
    u_x = 0 on the whole root face, u_y = 0 on the vertical mid-line,
    u_z = 0 on the horizontal mid-line. This matches the beam-theory clamp."""
    tol = 1e-6
    root = [n for n, (x, y, z) in nodes.items() if abs(x) < tol]
    root_y = [n for n in root if abs(nodes[n][1] - B / 2) < tol]
    root_z = [n for n in root if abs(nodes[n][2] - h / 2) < tol]
    tip = [n for n, (x, y, z) in nodes.items() if abs(x - L) < tol]

    def write_set(f, name, ids):
        f.write(f"*NSET, NSET={name}\n")
        for i in range(0, len(ids), 10):
            f.write(", ".join(str(n) for n in ids[i:i + 10]) + "\n")

    with open(path, "w") as f:
        f.write("*NODE, NSET=NALL\n")
        for n, (x, y, z) in nodes.items():
            f.write(f"{n}, {x:.6f}, {y:.6f}, {z:.6f}\n")
        f.write("*ELEMENT, TYPE=C3D8I, ELSET=EALL\n")
        for e, conn in elements:
            conn = beam.fix_orientation(nodes, conn)
            f.write(f"{e}, " + ", ".join(str(n) for n in conn) + "\n")
        write_set(f, "ROOT", root)
        write_set(f, "ROOTY", root_y)
        write_set(f, "ROOTZ", root_z)
        write_set(f, "TIP", tip)
        f.write(f"*MATERIAL, NAME=STEEL\n*ELASTIC\n{E}, {NU}\n")
        f.write("*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL\n")
        f.write("*STEP\n*STATIC\n*BOUNDARY\nROOT, 1, 1, 0.0\nROOTY, 2, 2, 0.0\nROOTZ, 3, 3, 0.0\n")
        f.write("*CLOAD\n")
        for n in tip:
            f.write(f"{n}, 3, {-P / len(tip):.8f}\n")
        f.write("*NODE PRINT, NSET=TIP\nU\n")
        f.write("*EL PRINT, ELSET=EALL\nS\n")
        f.write("*NODE FILE\nU\n*EL FILE\nS\n*END STEP\n")


def read_max_stress(dat_path):
    """Maximum |sigma_xx| (bending stress) over all integration points."""
    smax = 0.0
    reading = False
    with open(dat_path) as f:
        for line in f:
            if "stresses" in line:
                reading = True
                continue
            if reading:
                parts = line.split()
                if len(parts) == 8:
                    try:
                        int(parts[0])
                        smax = max(smax, abs(float(parts[2])))
                    except ValueError:
                        pass
    return smax


def evaluate(h):
    job = f"design_H{h}"
    nodes, elements = build_mesh(h)
    write_inp(nodes, elements, h, job + ".inp")
    res = subprocess.run([beam.CCX, job], capture_output=True, text=True)
    if res.returncode != 0 or not os.path.exists(job + ".dat"):
        print(res.stdout[-2000:])
        raise RuntimeError("CalculiX failed")
    delta = abs(beam.read_tip_deflection(job + ".dat"))
    sigma_ip = read_max_stress(job + ".dat")
    # Integration points are inside the elements, not on the surface. Bending stress
    # varies linearly through the height, so extrapolate the outermost Gauss point
    # value to the top/bottom surface.
    nz = 6
    c_ip = h / 2 - (h / nz) * (1 - 1 / 3 ** 0.5) / 2
    sigma = sigma_ip * (h / 2) / c_ip

    I = B * h ** 3 / 12
    delta_th = P * L ** 3 / (3 * E * I)
    sigma_th = P * L * (h / 2) / I
    mass = RHO * L * B * h
    return dict(h=h, mass=mass, delta=delta, delta_th=delta_th,
                sigma=sigma, sigma_th=sigma_th, sf=SY / sigma)


def main():
    print(f"Load {P:.0f} N | SF >= {SF_MIN} | deflection <= {DEFL_MAX:.2f} mm\n")
    print(f"{'H [mm]':>7} {'Mass [kg]':>10} {'Stress FEA':>11} {'Stress th.':>11} "
          f"{'SF':>6} {'Defl FEA':>9} {'Defl th.':>9} {'OK?':>5}")
    designs = []
    for h in HEIGHTS:
        d = evaluate(h)
        d["ok"] = d["sf"] >= SF_MIN and d["delta"] <= DEFL_MAX
        designs.append(d)
        print(f"{h:>7} {d['mass']:>10.4f} {d['sigma']:>11.1f} {d['sigma_th']:>11.1f} "
              f"{d['sf']:>6.2f} {d['delta']:>9.3f} {d['delta_th']:>9.3f} {'yes' if d['ok'] else 'no':>5}")

    valid = [d for d in designs if d["ok"]]
    if valid:
        best = min(valid, key=lambda d: d["mass"])
        print(f"\nLightest valid design: H = {best['h']} mm, mass = {best['mass']:.4f} kg, "
              f"SF = {best['sf']:.2f}, deflection = {best['delta']:.3f} mm")
    else:
        print("\nNo design meets the requirements. Increase the height range.")

    hs = [d["h"] for d in designs]
    fig, ax1 = plt.subplots(figsize=(8, 4.8))
    ax1.plot(hs, [d["sf"] for d in designs], "o-", color="tab:blue", label="Safety factor (FEA)")
    ax1.axhline(SF_MIN, color="tab:blue", ls="--", alpha=0.6, label=f"SF min = {SF_MIN}")
    ax1.set_xlabel("Beam height H [mm]")
    ax1.set_ylabel("Safety factor", color="tab:blue")
    ax2 = ax1.twinx()
    ax2.plot(hs, [d["delta"] for d in designs], "s-", color="tab:red", label="Tip deflection (FEA)")
    ax2.axhline(DEFL_MAX, color="tab:red", ls="--", alpha=0.6, label=f"Deflection limit = {DEFL_MAX:.2f} mm")
    ax2.set_ylabel("Tip deflection [mm]", color="tab:red")
    if valid:
        ax1.axvline(best["h"], color="green", ls=":", label=f"Lightest valid: H = {best['h']} mm")
    lines = ax1.get_legend_handles_labels()[0] + ax2.get_legend_handles_labels()[0]
    labels = ax1.get_legend_handles_labels()[1] + ax2.get_legend_handles_labels()[1]
    ax1.legend(lines, labels, loc="center right", fontsize=8)
    ax1.grid(alpha=0.3)
    plt.title("Design sweep: cantilever beam height")
    plt.tight_layout()
    plt.savefig("design_sweep.png", dpi=150)
    print("Plot saved: design_sweep.png")


if __name__ == "__main__":
    main()
