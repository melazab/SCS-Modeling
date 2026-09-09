"""Chunk 0 smoke test: is NEURON installed and can it compile a mod file?

Two things are checked, because they fail independently:

  1. The Python module imports and can integrate a cable equation. This is what
     `pip install neuron` gives you.
  2. `nrnivmodl` compiles a .mod file into a loadable mechanism. This needs a
     working C compiler, and it is what the MRG axon in the next chunk depends
     on -- MRG ships custom channel mechanisms as .mod files. A NEURON that
     imports fine but cannot compile is useless for this project, and the
     failure would otherwise not show up until much later.

Run:  src/neuron/.venv/bin/python src/neuron/smoke_test.py
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def check_import():
    from neuron import h
    h.load_file("stdrun.hoc")
    return h.nrnversion()


def check_spike():
    """A Hodgkin-Huxley soma should fire when injected with enough current.

    This is NOT a model of anything in the study -- HH is unmyelinated squid at
    6.3 C. It is here only because it is the cheapest end-to-end exercise of
    section, mechanism, stimulus and integrator.
    """
    from neuron import h
    h.load_file("stdrun.hoc")
    soma = h.Section(name="soma")
    soma.L = soma.diam = 20
    soma.insert("hh")

    stim = h.IClamp(soma(0.5))
    stim.delay, stim.dur, stim.amp = 1.0, 1.0, 5.0     # nA, well over threshold

    v = h.Vector().record(soma(0.5)._ref_v)
    t = h.Vector().record(h._ref_t)
    h.finitialize(-65)
    h.continuerun(10.0)

    vmax = max(v)
    peak_t = t[list(v).index(vmax)]
    return vmax, peak_t


MOD = """
NEURON {
    SUFFIX smoketest
    NONSPECIFIC_CURRENT i
    RANGE g, e
}
PARAMETER {
    g = 0.001 (siemens/cm2)
    e = -65   (millivolt)
}
ASSIGNED {
    v (millivolt)
    i (milliamp/cm2)
}
BREAKPOINT {
    i = g * (v - e)
}
"""


def check_nrnivmodl():
    """Compile a trivial passive mechanism and confirm it loads."""
    nrnivmodl = os.path.join(os.path.dirname(sys.executable), "nrnivmodl")
    if not os.path.exists(nrnivmodl):
        nrnivmodl = shutil.which("nrnivmodl")
    if not nrnivmodl:
        return False, "nrnivmodl not found"

    tmp = tempfile.mkdtemp(prefix="nrn_smoke_")
    try:
        with open(os.path.join(tmp, "smoketest.mod"), "w") as fh:
            fh.write(MOD)
        p = subprocess.run([nrnivmodl], cwd=tmp, capture_output=True, text=True)
        if p.returncode != 0:
            return False, (p.stdout + p.stderr)[-500:]

        # Load it in a *subprocess*: a mechanism cannot be unloaded once loaded,
        # and one compiled in a temp dir should not pollute this process.
        #
        # No nrn_load_dll here on purpose. NEURON auto-loads x86_64/libnrnmech.so
        # from the working directory at import, so loading it explicitly as well
        # raises "The user defined name already exists" -- which looks exactly
        # like a compile failure but is in fact proof the compile worked.
        probe = (
            "from neuron import h;"
            "s = h.Section();"
            "s.insert('smoketest');"
            "print('MECH_OK')"
        )
        q = subprocess.run([sys.executable, "-c", probe], cwd=tmp,
                           capture_output=True, text=True)
        return ("MECH_OK" in q.stdout), (q.stdout + q.stderr).strip()[-500:]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    ok = True

    ver = check_import()
    print("NEURON version      : %s" % ver.split("\n")[0])

    vmax, peak_t = check_spike()
    spiked = vmax > 0.0
    print("HH soma peak Vm     : %+.1f mV at t = %.2f ms  -> %s"
          % (vmax, peak_t, "SPIKE" if spiked else "NO SPIKE"))
    ok = ok and spiked

    compiled, detail = check_nrnivmodl()
    print("nrnivmodl compile   : %s" % ("OK" if compiled else "FAILED"))
    if not compiled:
        print("   %s" % detail)
    ok = ok and compiled

    print("\nRESULT: %s" % ("chunk 0 OK" if ok else "chunk 0 FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
