"""Run with runpy.run_path(..., run_name='__main__') inside FreeCAD's GUI.

Uses temporary empty documents and a sleeping worker, never an anatomical mesh.
"""
from pathlib import Path
import runpy
import tempfile


def run():
    import FreeCAD as App
    macro = Path(__file__).resolve().parents[1] / 'scripts' / 'SCS_Mesh_Generator.FCMacro'
    ui = runpy.run_path(str(macro), run_name='mesh_panel_regression')
    original = App.ActiveDocument
    documents = []
    stack = None
    with tempfile.TemporaryDirectory(prefix='scs-mesh-ui-') as tmp:
        try:
            a = App.newDocument('SCSMeshPanelTestA')
            documents.append(a.Name)
            stack = ui['DocumentMeshPanels']()
            pa = stack.currentWidget()
            pa.h_max.setValue(1.25)
            b = App.newDocument('SCSMeshPanelTestB')
            documents.append(b.Name)
            stack.sync_document()
            pb = stack.currentWidget()
            pb.h_max.setValue(2.0)
            assert all(w.text() == '-' for w in pb.mesh_report_labels.values())
            App.setActiveDocument(a.Name)
            stack.sync_document()
            assert stack.currentWidget() is pa and pa.h_max.value() == 1.25
            assert pb.h_max.value() == 2.0

            worker = Path(tmp) / 'worker.py'
            worker.write_text("import time\nprint('PROGRESS 12', flush=True)\ntime.sleep(30)\n")
            pa.mesh_runner.script = str(worker)
            pa.set_mesh_state('running')
            pa.mesh_runner.start([tmp])
            pa.refresh_gates()
            assert pa.btn_mesh.text() == 'Abort Mesh' and pa.btn_mesh.isEnabled()
            App.setActiveDocument(b.Name)
            stack.sync_document()
            assert pa.mesh_runner.is_running() and stack.other_job_running(pb)
            assert pb.mesh_progress.value() == 0
            App.setActiveDocument(a.Name)
            stack.sync_document()
            pa.btn_mesh.click()
            assert not pa.mesh_runner.is_running()
            assert pa.mesh_state == 'aborted' and pa.btn_mesh.text() == 'Generate Mesh'
            pa.on_mesh_failed('test failure')
            assert pa.mesh_state == 'failed' and pa.btn_mesh.text() == 'Generate Mesh'
            assert pb.mesh_state == 'idle'
            closed_name = a.Name
            App.closeDocument(closed_name)
            stack.sync_document()
            assert closed_name not in stack.panels
            print('PASS: per-document controls, job ownership, abort/failure transitions, document close')
        finally:
            if stack is not None:
                stack.on_dock_closing()
                stack.deleteLater()
            for name in documents:
                if name in App.listDocuments():
                    App.closeDocument(name)
            if original is not None:
                App.setActiveDocument(original.Name)


if __name__ == '__main__':
    run()
