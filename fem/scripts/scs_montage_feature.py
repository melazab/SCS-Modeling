"""App::FeaturePython proxy for `SCS_Montage` -- the interactive per-contact
current control for the "SCS Potential Field Visualizer".

Mohamed edits N float properties (mA, one per lead contact; + = anode/
source, - = cathode/sink, matching assign_and_solve.py's own convention) in
the object's Data tab -- N = however many contacts the CURRENT live lead
has (ensure_contact_properties() below keeps this in sync as leads are
swapped; see its own docstring -- this used to be a fixed 8 regardless of
the actual lead). Recompute re-samples the N-basis-field superposition
(scs_montage.py -- pure numpy/scipy/vtk, no FreeCAD dependency, importable
and unit-tested from a plain terminal) onto the real epidural/dura/CSF/white/grey STL
surfaces and reloads the five sibling `Fem::FemPostPipeline` objects this
module does not own or create (see scs_potential_viz.FCMacro for that half).

Kept in its own importable module, not written inline in the macro, so
FreeCAD can reconstruct this object's Proxy after a document save/reload --
an App::FeaturePython's Proxy is persisted by module+class reference, not by
value.

WHY THE RELOAD LOGIC LIVES IN THE VIEW PROVIDER'S updateData, NOT execute()

`execute()` runs inside the document's own recompute transaction. Calling
`doc.recompute()` again from inside it (needed to make the sibling pipeline
objects' ViewObject.Field enum populate/repaint -- see the stale-colour-bar
comment below) is reentrant into that same transaction, which this repo's
own rule (fem/README.md, HANDOFF.md: rehearse before trusting an FEM-adjacent
change) says not to assume is safe without checking. So `execute()` here
does ONLY the plain-Python half (montage arithmetic, writing the .vtp files,
recording a Status string) and touches no other document object; the
ViewProvider observes a monotonically increasing PlotRevision and schedules
the Gui-facing half with a zero-delay Qt callback, after recompute returns.
This reloads the pipelines even when successive montages have identical
summary extrema, without running the montage arithmetic twice.
"""
import os
import re
import sys
import json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scs_montage as M

# Must match the Fem::FemPostPipeline object names created by
# scs_potential_viz.FCMacro.
PIPE_NAMES = {
    "epidural": "SCS_Epidural_PotentialField",
    "dura": "SCS_Dura_PotentialField",
    "csf": "SCS_CSF_PotentialField",
    "white": "SCS_White_PotentialField",
    "grey": "SCS_Grey_PotentialField",
}

PIPE_NAMES.update({t: 'SCS_%s_PotentialField' % t for t in M.TISSUES if t not in PIPE_NAMES})

def group_lead_pipes(doc, contact_map):
    """Move existing field pipelines into persistent per-lead tree groups."""
    root = doc.getObject('SCS_PotentialFieldViz')
    if root is None:
        return
    mapping = {str(k): v for k, v in contact_map.items()}
    leads = sorted({int(v[0]) for v in mapping.values()})
    groups = {}
    for lead in leads:
        name = 'SCS_Lead%d_PotentialFields' % lead
        group = doc.getObject(name)
        if group is None:
            group = doc.addObject('App::DocumentObjectGroup', name)
        group.Label = 'Lead %d' % lead
        if group not in root.Group:
            root.addObject(group)
        groups[lead] = group
    for obj in list(doc.Objects):
        lead, label = None, None
        for prefix, suffix in (('SCS_contact_', 'contact'), ('SCS_insulator_', 'insulator')):
            if obj.Name.startswith(prefix) and obj.Name.endswith('_PotentialField'):
                index = obj.Name[len(prefix):-len('_PotentialField')]
                if suffix == 'contact' and index in mapping:
                    lead, contact = mapping[index]
                    label = 'Contact %d -- SCS potential field' % contact
                elif suffix == 'insulator' and index.isdigit():
                    lead, label = int(index), 'Insulation -- SCS potential field'
        if obj.Name == 'SCS_insulator_PotentialField' and len(leads) == 1:
            lead, label = leads[0], 'Insulation -- SCS potential field'
        if lead not in groups:
            continue
        visible = obj.ViewObject.Visibility
        for parent in list(obj.InList):
            if parent == root or parent in groups.values():
                if obj in parent.Group and parent != groups[lead]:
                    parent.removeObject(obj)
        if obj not in groups[lead].Group:
            groups[lead].addObject(obj)
        obj.Label = label
        obj.ViewObject.Visibility = visible


_CONTACT_PROP_RE = re.compile(r"^Contact(\d+)_mA$")


def contact_prop_name(i):
    return "Contact%d_mA" % i


def contact_indices(fp):
    """Sorted contact numbers `fp` (an SCS_Montage object) currently has a
    Contact<N>_mA property for -- DISCOVERED by property-name pattern, not
    assumed to be 1..8. 2026-09-14 correction: this object used to always
    create exactly 8 fixed properties (Contact1_mA..Contact8_mA) in
    __init__, regardless of how many contacts the live lead it was built
    for actually had. Every reader of "how many contacts does this montage
    have, and what are they" (execute(), the ViewProvider's updateData())
    now goes through this function instead of a hardcoded range(1, 9), so
    it reflects whatever ensure_contact_properties() last set it up for."""
    out = []
    for name in fp.PropertiesList:
        m = _CONTACT_PROP_RE.match(name)
        if m:
            out.append(int(m.group(1)))
    return sorted(out)


class SCSMontageProxy:
    """Proxy for the App::FeaturePython 'SCS_Montage' object. App-layer only:
    reads the selected lead's currents, runs the montage, writes the .vtp files, records a
    human-readable Status. Never touches another document object or the Gui
    -- see module docstring for why.

    EVERY CONTACT DEFAULTS TO 0 mA -- no hardcoded source/sink. Mohamed's
    own words on an earlier version of this panel: "it shouldn't assume
    what are the cathodes, anodes and inactive electrodes, that's for the
    user to select!" `SCS_Potential_Visualizer.FCMacro`'s panel computes a
    role label (source/sink/floating) live from each contact's CURRENT
    typed value, purely for display -- nothing here or there hardcodes
    config.SOURCE_CONTACT/SINK_CONTACT as a starting state. That panel's
    own "Reset to validated bipolar" convenience button was removed
    2026-09-14 per Mohamed's own review ("I don't like the reset to
    validated bipolar, dude") -- config.SOURCE_CONTACT/SINK_CONTACT/
    CURRENT_A now only mean anything as `solve_basis_fields()`'s own
    combination for computing the stored `V` array (see
    assign_and_solve.py) and as `scs_montage.default_bipolar_mA()`'s CLI
    default; there is no UI shortcut back to that specific case any more --
    typing contact 3 = +1000 mA / contact 5 = -1000 mA reproduces it by
    hand.
    """

    def __init__(self, fp, n_contacts=0):
        """`n_contacts` -- how many Contact<i>_mA properties to create up
        front. Defaults to 0 (none) because this object is normally created
        by SCS_Potential_Visualizer.FCMacro's ensure_pipes() BEFORE the
        panel has polled the live document for a lead (see that macro's
        PotentialVisualizerPanel.__init__ / refresh_gates()) -- the real
        count is set moments later via ensure_contact_properties() once the
        live lead (and its actual contact count) is known, and again
        whenever that count changes (a different lead swapped in)."""
        fp.Proxy = self
        self.ensure_contact_properties(fp, n_contacts)
        fp.addProperty(
            "App::PropertyString",
            "SolutionNpz",
            "Montage",
            "Path to the solution.npz this montage samples -- one lead's "
            "fem/out/lead_runs/<run_key>/solution.npz, keyed by a content "
            "fingerprint of that lead's live geometry (see "
            "fem/scripts/live_lead.py). Set by SCS_Potential_Visualizer."
            "FCMacro's panel once its self-sufficient 'Plot Field' action "
            "has a solved field for the CURRENT live lead (meshing/solving "
            "it first if needed), never edited here directly.",
        )
        fp.SolutionNpz = ""
        fp.addProperty(
            "App::PropertyFloat",
            "ClampVolts",
            "Montage",
            "Bulk-tissue colour-legend saturation, +/-volts (see the "
            "V_clamped_* array note in scs_montage.py). 0 = auto-choose "
            "from this montage's own data (95th percentile |V|).",
        )
        fp.ClampVolts = 0.0
        fp.addProperty(
            "App::PropertyString",
            "Status",
            "Montage",
            "Last recompute's summary (informational, read-only).",
        )
        fp.setEditorMode("Status", 1)  # read-only in the property editor
        self.ensure_revision(fp)
        import scs_waveform
        scs_waveform.ensure_properties(fp)
        fp.WaveformMode = "FAST-style biphasic"

    @staticmethod
    def ensure_revision(fp):
        if 'PlotRevision' not in fp.PropertiesList:
            fp.addProperty('App::PropertyInteger', 'PlotRevision', 'Montage',
                           'Successful plot update counter')
        fp.setEditorMode('PlotRevision', 2)

    @staticmethod
    def validate_binding(fp):
        import artifacts
        import live_lead
        path = fp.SolutionNpz
        if not os.path.isfile(path):
            raise RuntimeError('Saved FEM cache is missing: %s. Restore the cache directory to reuse this solve.' % path)
        run_dir = os.path.dirname(path)
        background = os.path.basename(path) == 'solution_bg.npz'
        if not artifacts.valid(path, artifacts.solution_signature(run_dir, background)):
            raise RuntimeError('Cached solution is missing or stale. Use Solve / Plot Field.')
        info = live_lead.find_live_lead(fp.Document)
        if info['status'] != live_lead.LiveLeadStatus.OK:
            raise RuntimeError('Cannot match the plotted field to the current lead.')
        with open(os.path.join(run_dir, 'params.json')) as f:
            params = json.load(f)
        params.pop('LEAD_DIR', None)
        import scs_document_cache
        if not scs_document_cache.geometry_matches(live_lead.tessellate_lead_arrays(info), run_dir):
            raise RuntimeError('Lead geometry changed. Generate the matching mesh and field.')
        scs_document_cache.bind_mesh(fp.Document, run_dir, params)

    def __getstate__(self):
        # Plot statistics and GUI timers are transient; paths/currents live in
        # document properties and the verified numerical cache stays on disk.
        return None

    def __setstate__(self, state):
        pass

    def onDocumentRestored(self, fp):
        self.ensure_revision(fp)
        import scs_waveform
        scs_waveform.ensure_properties(fp)

    @staticmethod
    def ensure_contact_properties(fp, n_contacts):
        """Add/remove Contact<i>_mA properties so `fp` ends up with EXACTLY
        {1..n_contacts}, whatever it had before. Called from __init__ (a
        fresh object, nothing to remove) AND from
        SCS_Potential_Visualizer.FCMacro's panel whenever the live lead's
        own contact count changes on an ALREADY-EXISTING SCS_Montage object
        -- so swapping a 4-contact lead for a 6-contact one adds
        Contact5_mA/Contact6_mA (starting at 0 mA, per this class's own
        "no contact starts pre-assigned a role" rule) without disturbing
        Contact1_mA..Contact4_mA's current values, and swapping down again
        removes the now-meaningless extras instead of leaving stale
        properties that don't correspond to any real contact on the
        current lead."""
        have = set(contact_indices(fp))
        want = set(range(1, n_contacts + 1))
        if have != want and "SolutionNpz" in fp.PropertiesList:
            # A basis for the previous contact layout cannot serve the new one.
            fp.SolutionNpz = ""
        for i in sorted(want - have):
            fp.addProperty(
                "App::PropertyFloat",
                contact_prop_name(i),
                "Montage",
                "Current into contact %d, mA. Positive = source/anode, "
                "negative = sink/cathode (assign_and_solve.py's own "
                "convention). Starts at 0 -- no contact is assumed to be "
                "a source or sink." % i,
            )
            setattr(fp, contact_prop_name(i), 0.0)
        for i in sorted(have - want):
            fp.removeProperty(contact_prop_name(i))

    @staticmethod
    def ensure_contact_mapping(fp, info):
        mapping = json.dumps(info['contact_map'], sort_keys=True)
        old = getattr(fp, 'ContactMap', '')
        # Changing layout must not silently transfer currents to other leads.
        if old and old != mapping:
            for i in contact_indices(fp):
                setattr(fp, contact_prop_name(i), 0.0)
            fp.SolutionNpz = ''
        SCSMontageProxy.ensure_contact_properties(fp, info['n_contacts'])
        if 'ContactMap' not in fp.PropertiesList:
            fp.addProperty('App::PropertyString', 'ContactMap', 'Montage',
                           'Solver contact ID to [lead number, local contact number]')
        if fp.ContactMap != mapping:
            fp.ContactMap = mapping
        fp.setEditorMode('ContactMap', 1)
        for i in info['contacts']:
            lead, local = info['contact_map'][i]
            fp.setGroupOfProperty(contact_prop_name(i), 'Lead %d currents' % lead)
            fp.setDocumentationOfProperty(contact_prop_name(i),
                'Lead %d / contact %d: first-phase current (mA)' % (lead, local))

    def execute(self, fp):
        # Opening the panel and document recomputes happen before Solve / Plot.
        # An empty binding means idle, never the historical CLI default field.
        if not fp.SolutionNpz:
            fp.Status = "No field selected. Use Solve / Plot Field for the current lead."
            return
        try:
            self.validate_binding(fp)
        except Exception as exc:
            fp.Status = str(exc)
            return
        ids = contact_indices(fp)
        currents_mA = [getattr(fp, contact_prop_name(i)) for i in ids]
        clamp_v = None if fp.ClampVolts == 0.0 else fp.ClampVolts
        sol_npz = fp.SolutionNpz
        out_dir = os.path.join(os.path.dirname(sol_npz), "scs_viz")
        try:
            import scs_waveform as W
            W.ensure_properties(fp)
            factor, phase_name = 1., 'Static'
            if fp.WaveformMode != 'Static':
                W.require_balanced(currents_mA)
                factor, phase_name = W.phase(fp.WaveformTimeMs, fp.FrequencyHz,
                                             fp.PhaseWidthUs, fp.InterphaseGapUs)
                # Fix the scale to the first-phase field, including off frames.
                if clamp_v is None:
                    clamp_v = M.choose_clamp_v(currents_mA, sol_npz=sol_npz)
            source, sink = W.totals(currents_mA)
            stats = M.run_montage([factor*v for v in currents_mA], out_dir=out_dir, clamp_v=clamp_v, sol_npz=sol_npz)
        except Exception as exc:
            fp.Status = 'Field update failed: %s' % exc
            return
        stats.update(phase=phase_name, factor=factor, source_mA=source, sink_mA=sink)
        self._last_stats = stats
        bits = [phase_name, "source=%.4g mA; sink=%.4g mA" % (source, sink),
            "clamp=+/-%.3gV" % stats["clamp_v"],
            "charge_sum=%.4gmA" % stats["charge_sum_mA"],
        ]
        for t in M.TISSUES:
            s = stats["tissues"][t]
            bits.append(
                "%s=[%.4g,%.4g]V" % (t, s["v_min"], s["v_max"])
                + (" NaN=%d" % s["n_nan"] if s["n_nan"] else "")
            )
        fp.Status = "; ".join(bits)
        self.ensure_revision(fp)
        fp.PlotRevision += 1


class SCSMontageViewProxy:
    """ViewProvider for SCS_Montage. Holds the Gui-facing half of a
    recompute: reload the three sibling FemPostPipeline objects and apply
    the stale-colour-bar-range fix (see view_elmer_result.FCMacro, which
    found this bug first: setting ViewObject.Field once immediately after
    a pipeline's data exists captures the -0.5/+0.5 default range and never
    updates again, even though the underlying data is correct -- the fix is
    Field="None", recompute, THEN Field=<real name>, recompute)."""

    def __init__(self, vp):
        vp.Proxy = self

    def attach(self, vp):
        self.Object = vp.Object

    def updateData(self, fp, prop):
        # fp.Status is the LAST property execute() sets, so this fires once
        # per recompute, after the montage/vtp-writing work is already done.
        if prop == 'PickFieldInViewport':
            import scs_field_display
            scs_field_display.apply(fp.Document)
            return
        if prop != "PlotRevision" or not fp.SolutionNpz:
            return
        import FreeCAD

        if not FreeCAD.GuiUp:
            return
        # Defer pipeline mutations until the document recompute has returned.
        from PySide import QtCore
        if getattr(self, '_pending', False):
            return
        self._pending = True
        QtCore.QTimer.singleShot(0, lambda: self.reload_pipes(fp))

    def reload_pipes(self, fp):
        self._pending = False
        stats = getattr(fp.Proxy, '_last_stats', None)
        if not stats or not fp.SolutionNpz:
            return
        # SAME sol_npz/out_dir resolution as execute() above -- this was
        # missing when SolutionNpz was first added (caught live this
        # session): without it, the 3D pipelines silently kept showing the
        # frozen fem/out/solution.npz field forever, no matter which lead
        # was selected, even though the readout stats (computed separately
        # by the panel) correctly reflected the selected lead. Both halves
        # of a recompute MUST resolve the same file the same way.
        sol_npz = fp.SolutionNpz
        out_dir = os.path.join(os.path.dirname(sol_npz), "scs_viz")
        doc = fp.Document
        active_names = {PIPE_NAMES.get(t, 'SCS_%s_PotentialField' % t)
                        for t in stats['tissues']}
        group = doc.getObject('SCS_PotentialFieldViz')
        descendants = list(group.Group)
        for child in group.Group:
            if child.TypeId == 'App::DocumentObjectGroup':
                descendants.extend(child.Group)
        for sibling in descendants:
            if sibling.TypeId == 'Fem::FemPostPipeline' and sibling.Name not in active_names:
                sibling.ViewObject.Visibility = False
        created_lead = False
        for tissue in stats['tissues']:
            pipe_name = PIPE_NAMES.get(tissue, 'SCS_%s_PotentialField' % tissue)
            pipe = doc.getObject(pipe_name)
            if pipe is None:
                pipe = doc.addObject('Fem::FemPostPipeline', pipe_name)
                pipe.Label = (M.TISSUE_LABELS.get(tissue, 'Lead ' + tissue.replace('_', ' ')) + ' -- SCS potential field')
                doc.getObject('SCS_PotentialFieldViz').addObject(pipe)
                created_lead = tissue.startswith('contact_') or tissue.startswith('insulator') or created_lead
            if tissue.startswith('contact_'):
                mapping = json.loads(getattr(fp, 'ContactMap', '') or '{}')
                pair = mapping.get(tissue.split('_')[-1])
                if pair:
                    pipe.Label = 'Lead %d / contact %d -- SCS potential field' % tuple(pair)
            visible = pipe.ViewObject.Visibility
            pipe.read(stats["tissues"][tissue]["path"])
            doc.recompute([pipe])
            # Field's enumeration list is stale here -- still the PREVIOUS
            # vtp's array names (e.g. "V_clamped_75V" after a reload that
            # now has "V_clamped_20V") immediately after recompute()
            # returns, so the fallback to "V_volts" below would silently
            # fire on every single edit. `FreeCADGui.updateGui()` (pumping
            # the Qt event loop) does force the refresh, but tested live
            # -- via the MCP execute_code bridge this was built through --
            # it twice made that bridge's GUI dispatch hang for 90+ seconds.
            # Whether that risk is specific to this RPC bridge or would
            # also bite a real interactive edit in Mohamed's own session is
            # untested; a DisplayMode round-trip forces the same internal
            # rebuild (this is what the "must be Surface before Field is
            # usable" comment in view_elmer_result.FCMacro is really
            # describing) and was verified not to hang either way, so use
            # that instead and avoid the question.
            pipe.ViewObject.DisplayMode = "Wireframe"
            doc.recompute([pipe])
            pipe.ViewObject.DisplayMode = "Surface"
            doc.recompute([pipe])
            clamp_field = "V_clamped_%gV" % stats["clamp_v"]
            enum = pipe.ViewObject.getEnumerationsOfProperty("Field")
            target = clamp_field if clamp_field in enum else ("V_volts" if "V_volts" in enum else None)
            if target is None:
                print("WARNING: no expected Field on %s -- got %s" % (pipe_name, enum))
                continue
            pipe.ViewObject.Field = "None"
            doc.recompute([pipe])
            pipe.ViewObject.Field = target
            doc.recompute([pipe])
            pipe.ViewObject.Visibility = visible
        import scs_field_display
        scs_field_display.apply(doc)
        group_lead_pipes(doc, json.loads(getattr(fp, 'ContactMap', '') or '{}'))
        if created_lead:
            # Coincident CAD surfaces otherwise obscure the coloured lead.
            import live_lead
            info = live_lead.find_live_lead(doc)
            if info['status'] == live_lead.LiveLeadStatus.OK:
                for obj in list(info['contacts'].values()) + list(info['insulators'].values()):
                    obj.ViewObject.Visibility = False
        print("SCS_Montage recompute done: " + fp.Status)

    def getIcon(self):
        return ":/icons/FEM_PostFilterClip.svg"

    def __getstate__(self):
        return None

    def __setstate__(self, state):
        return None
