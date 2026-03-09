"""
/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the 3-Clause BSD License, (the "License").
******************************************************************************/
"""

import math
import os


MATERIAL_FILES = {
    "Graphene": "GrafUnitCell.in",
    "Boron Nitride": "BNUnitCell.in",
    "Phosphorene": "PhosUnitCell.in",
    "MoS\u2082": "MoS2UnitCell.in",
    "Ti\u2083C\u2082": "Ti3C2UnitCell.in",
}

MATERIALS = list(MATERIAL_FILES.keys())


def _get_unit_cell_path(material_name):
    pkg_dir = os.path.dirname(__file__)
    filename = MATERIAL_FILES[material_name]
    return os.path.join(pkg_dir, "UnitCellExamples", filename)


def _setup_chiraltube_globals(ct, unit_cell_path):
    """Load unit cell file and initialise chiraltube module-level globals."""
    nat, A, UnitCell = ct.read_arch(unit_cell_path)
    UnitCell = ct.shift_cell(A, UnitCell)

    ct.a1 = A[0]
    ct.a2 = A[1]
    ct.a3 = A[2]
    ct.a2 = ct.a2 * (-1)
    for i in range(len(UnitCell)):
        UnitCell[i] = ct.a2 + UnitCell[i]
    ct.UnitCell = UnitCell
    ct.phi = math.acos(
        (ct.a1.x * ct.a2.x + ct.a1.y * ct.a2.y + ct.a1.z * ct.a2.z)
        / (ct.a1.mag() * ct.a2.mag())
    )


def generate(opts):
    from . import chiraltube as ct

    options = opts.get("options", {})
    n = int(options.get("n", 5))
    m = int(options.get("m", 5))
    repeats = int(options.get("repeats", 1))
    material = options.get("material", MATERIALS[0])
    unit_cell_path = _get_unit_cell_path(material)

    _setup_chiraltube_globals(ct, unit_cell_path)

    theta, Res = ct.robtainxy(n, m)
    if not Res:
        raise ValueError(
            f"No integer solutions found for (n,m)=({n},{m}). "
            "Try different chiral indices."
        )
    x, y, _, _ = Res[0]
    x, y = round(x), round(y)

    Arr = ct.arr_initial(n, m, x, y)
    Arr = ct.rotate(Arr, theta)
    Arr, disy, disx = ct.eliminate(Arr, n, m, x, y)
    Arr, radio = ct.nanotube(Arr, n, m)

    # Remove atoms that coincide with periodic images on the boundary
    for ele in list(Arr):
        displaced = ct.point3d(0, 0, -disy, "") + ele
        for ele2 in Arr:
            if ct.distance(displaced, ele2) < 0.2:
                Arr.remove(ele)
                break

    # repeats=1 → original unit only; repeats=2 → doubled; etc.
    rep = repeats - 1
    if rep > 0:
        Arr = ct.repeat(Arr, disy, rep)
    total_height = disy * repeats

    lines = [str(len(Arr))]
    lines.append(
        f"NT (n,m)=({n},{m}). z={total_height:.4f} A, r={radio:.4f} A. "
        f"T(x,y)=({x},{y}) material={material}"
    )
    for ele in Arr:
        lines.append(f"{ele.ele}    {ele.x:.6f}    {ele.y:.6f}    {ele.z:.6f}")

    return "\n".join(lines) + "\n"


def run(avo_input):
    return {
        "append": True,
        "moleculeFormat": "xyz",
        "xyz": generate(avo_input),
    }
