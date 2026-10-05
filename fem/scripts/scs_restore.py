"""Restore lightweight callbacks, without solving or recomputing on file open."""
_observer = None


def repair(doc):
    import FreeCAD
    import scs_montage_feature as F
    obj = doc.getObject('SCS_Montage')
    if obj is None:
        return
    if not isinstance(obj.Proxy, F.SCSMontageProxy):
        obj.Proxy = F.SCSMontageProxy.__new__(F.SCSMontageProxy)
    F.SCSMontageProxy.ensure_revision(obj)
    if FreeCAD.GuiUp:
        import scs_field_display
        scs_field_display.ensure_setting(obj)
        scs_field_display.apply(doc)
    if FreeCAD.GuiUp and not isinstance(obj.ViewObject.Proxy, F.SCSMontageViewProxy):
        F.SCSMontageViewProxy(obj.ViewObject)


def release_closed_document_cache(closed_doc):
    import FreeCAD
    import scs_montage as M
    paths = []
    for doc in FreeCAD.listDocuments().values():
        if doc is closed_doc:
            continue
        montage = doc.getObject('SCS_Montage')
        if montage is not None:
            paths.append(getattr(montage, 'SolutionNpz', ''))
    M.release_unused_basis(paths)


class RestoreObserver:
    def slotDeletedDocument(self, doc):
        release_closed_document_cache(doc)

    def slotFinishRestoreDocument(self, doc):
        repair(doc)


def install():
    import FreeCAD
    global _observer
    if _observer is None:
        _observer = RestoreObserver()
        FreeCAD.addDocumentObserver(_observer)
    for doc in FreeCAD.listDocuments().values():
        repair(doc)
