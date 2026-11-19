"""
Step 2: mesh independence study.
Runs the cantilever beam from beam_v1.py on four meshes, refined by a constant
ratio of 2, and checks when the tip deflection stops changing.
Also estimates the discretization uncertainty with the Grid Convergence Index (GCI).
"""
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import beam_v1 as beam

MESHES = [(10, 1, 1), (20, 2, 2), (40, 4, 4), (80, 8, 8)]  # (nx, ny, nz), ratio r = 2


def main():
    I = beam.B * beam.H ** 3 / 12
    delta_theory = beam.P * beam.L ** 3 / (3 * beam.E * I)

    results = []
    print(f"{'Mesh':>12} {'Elements':>9} {'Deflection [mm]':>16} {'Change vs previous':>19}")
    previous = None
    for nx, ny, nz in MESHES:
        delta, n_el = beam.run(nx, ny, nz, job=f"beam_{nx}x{ny}x{nz}")
        delta = abs(delta)
        change = "" if previous is None else f"{(delta - previous) / previous * 100:+.3f} %"
        print(f"{f'{nx}x{ny}x{nz}':>12} {n_el:>9} {delta:>16.5f} {change:>19}")
        results.append((n_el, delta))
        previous = delta

    # GCI with the three finest meshes (Roache / ASME V&V 20 style)
    f3, f2, f1 = results[-3][1], results[-2][1], results[-1][1]  # coarse, medium, fine
    r = 2.0
    p = math.log(abs((f3 - f2) / (f2 - f1))) / math.log(r)       # observed order
    f_exact = f1 + (f1 - f2) / (r ** p - 1)                      # Richardson extrapolation
    gci_fine = 1.25 * abs((f1 - f2) / f1) / (r ** p - 1) * 100   # safety factor 1.25

    print()
    print(f"Observed order of accuracy p:   {p:.2f}")
    print(f"Richardson extrapolated value:  {f_exact:.5f} mm")
    print(f"GCI (fine mesh):                {gci_fine:.3f} %")
    print(f"Euler-Bernoulli theory:         {delta_theory:.5f} mm")
    print(f"Fine mesh vs theory:            {(f1 - delta_theory) / delta_theory * 100:+.2f} %")

    n_el = [res[0] for res in results]
    defl = [res[1] for res in results]
    plt.figure(figsize=(7, 4.5))
    plt.semilogx(n_el, defl, "o-", label="CalculiX (C3D8I)")
    plt.axhline(delta_theory, color="gray", ls="--", label="Euler-Bernoulli theory")
    plt.axhline(f_exact, color="green", ls=":", label="Richardson extrapolation")
    plt.xlabel("Number of elements")
    plt.ylabel("Tip deflection [mm]")
    plt.title("Mesh independence study: cantilever beam")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("mesh_study.png", dpi=150)
    print("\nPlot saved: mesh_study.png")


if __name__ == "__main__":
    main()
