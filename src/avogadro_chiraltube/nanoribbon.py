"""
/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the 3-Clause BSD License, (the "License").
******************************************************************************/
"""

from .nanotube import MATERIALS, _get_unit_cell_path, _setup_chiraltube_globals


def generate(opts):
    from . import chiraltube as ct

    options = opts.get("options", {})
    n = int(options.get("n", 5))
    m = int(options.get("m", 5))
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

    # Remove atoms on the periodic boundaries (x, y, and diagonal)
    for ele in list(Arr):
        displaced = ct.point3d(0, -disy, 0, "") + ele
        for ele2 in Arr:
            if ct.distance(displaced, ele2) < 0.2:
                Arr.remove(ele)
                break
    for ele in list(Arr):
        displaced = ct.point3d(-disx, 0, 0, "") + ele
        for ele2 in Arr:
            if ct.distance(displaced, ele2) < 0.2:
                Arr.remove(ele)
                break
    for ele in list(Arr):
        displaced = ct.point3d(-disx, -disy, 0, "") + ele
        for ele2 in Arr:
            if ct.distance(displaced, ele2) < 0.2:
                Arr.remove(ele)
                break

    lines = [str(len(Arr))]
    lines.append(
        f"Nanoribbon (n,m)=({n},{m}). Y={disy:.4f} A, X={disx:.4f} A "
        f"material={material}"
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
