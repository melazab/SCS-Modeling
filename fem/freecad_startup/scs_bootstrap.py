"""Bootstrap from an imported module: Init.py's __file__ can name FreeCAD's loader."""
import os
import sys


def initialize():
    repo = os.environ.get('SCS_MODELING_REPO', os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
    scripts = os.path.join(repo, 'fem', 'scripts')
    if scripts not in sys.path:
        sys.path.insert(0, scripts)
    # Explicitly load our installed addon's callbacks before document restore.
    # FreeCAD's restore policy permits already-loaded modules; do not broaden
    # its trusted directory list or disable its import checks.
    import scs_montage_feature
