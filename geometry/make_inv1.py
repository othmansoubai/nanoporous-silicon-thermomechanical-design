#!/usr/bin/env python3
# Minimal builder for the INV1 inverse-design geometry.
# Uses the included independent visualization constructor for INV1.
import os
import build_structures as bs

# INV1_AR2.5_s37 : Pore A x[0,5] y[0,2], Pore B x[3,8] y[3,5]  (offset +3 = 37.5% stagger, AR=2.5)
bulk  = bs.build_bulk_lattice()
atoms = bs.carve_pores(bulk, 0, 5, 0, 2, 3, 8, 3, 5)

out_dir = os.path.join(os.path.dirname(__file__), "generated_structures")
os.makedirs(out_dir, exist_ok=True)
out = os.path.join(out_dir, "structure_INV1_AR2.5_s37.data")
bs.write_lammps_data(atoms, out, "INV1_AR2.5_s37")

phi = 100.0 * (1.0 - len(atoms) / len(bulk))
print(f"wrote {out}")
print(f"  {len(atoms)} atoms of {len(bulk)}  ->  phi = {phi:.2f}%")
