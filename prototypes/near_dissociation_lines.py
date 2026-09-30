"""Build the line list to the B dissociation limit (19 600-20 050 cm-1) and say how far it now reaches.

The default line list stops near v' = 62 (19 955 cm-1): its B box ends at 7 A and every level beyond is a
box state. With the 12 A box of intensity_model(grid=NEAR_DISSOCIATION) every level bound by more than
0.3 cm-1 is converged (docs/design/near-dissociation.md). The result is cached, for the 19 700-20 035 volume.
"""
import time
import numpy as np
from i2spec.intensity import intensity_model, master_line_list, NEAR_DISSOCIATION

if __name__ == "__main__":
    for label, kw in (("default grid", {}), ("12 A grid", dict(grid=NEAR_DISSOCIATION, nlev_b=95))):
        t = time.time()
        m = intensity_model("127I2", **kw)
        ml = master_line_list(m, 19600.0, 20050.0, T_range=(250.0, 350.0), S_min=1e-28)
        ll = ml.at(295.0, 19600.0, 20050.0, S_min=1e-26)
        vu = np.asarray(ll.v_upper)
        print(f"{label}: {time.time() - t:.0f}s; {len(ll.nu)} lines above 1e-26 cm at 295 K; "
              f"v' up to {vu.max()}; highest line {np.max(ll.nu):.2f} cm-1; lines above 19 955: {np.sum(np.asarray(ll.nu) > 19955)}",
              flush=True)
