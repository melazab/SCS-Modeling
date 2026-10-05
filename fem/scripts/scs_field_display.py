"""Lightweight field interaction settings; never resample or modify FEM data."""
PICK_NODE_NAME = 'SCS_FieldPickStyle'


def ensure_setting(montage):
    if 'PickFieldInViewport' not in montage.PropertiesList:
        montage.addProperty('App::PropertyBool', 'PickFieldInViewport', 'Display',
                            'Allow mouse picking on field surfaces; slower for dense models. Tree selection remains available.')
        montage.PickFieldInViewport = False


def apply(doc, enabled=None):
    from pivy import coin
    montage = doc.getObject('SCS_Montage')
    if enabled is None:
        enabled = bool(getattr(montage, 'PickFieldInViewport', False))
    for obj in doc.Objects:
        if obj.TypeId != 'Fem::FemPostPipeline' or not obj.Name.startswith('SCS_'):
            continue
        root = obj.ViewObject.RootNode
        # Reuse our node across repeated updates, never accumulate scene nodes.
        node = next((root.getChild(i) for i in range(root.getNumChildren())
                     if str(root.getChild(i).getName()) == PICK_NODE_NAME), None)
        if node is None:
            node = coin.SoPickStyle()
            node.setName(PICK_NODE_NAME)
            root.insertChild(node, 0)
        node.style = coin.SoPickStyle.SHAPE if enabled else coin.SoPickStyle.UNPICKABLE


def set_enabled(doc, enabled):
    montage = doc.getObject('SCS_Montage')
    if montage is not None:
        ensure_setting(montage)
        montage.PickFieldInViewport = bool(enabled)
    apply(doc, enabled)
