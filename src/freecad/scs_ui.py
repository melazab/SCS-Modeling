"""Shared widgets for the SCS panels, so they look and behave alike.

The rules they implement, learned on lead_designer.FCMacro:

* One piece of state, one control. A thing that is either in the document or
  not gets a single checkable button whose text says the state ("Show tracts"
  / "Tracts shown"), green while on, never a Show button beside a Clear button.
* Theme-proof: the off state is painted from the widget's own palette and the
  on state paints its own green and white, so both read on FreeCAD's light and
  dark stylesheets (see lead_designer.FCMacro toggle_style for the history).
* Explanations go in tooltips; the panel itself carries one grey hint line.
* Mutually exclusive choices are a segmented row of buttons, not a combo box
  hiding two options.

Import from a macro with the repo's src/freecad on sys.path.
"""
try:
    from PySide import QtCore, QtGui, QtWidgets
except ImportError:                                     # pragma: no cover
    try:
        from PySide6 import QtCore, QtGui, QtWidgets
    except ImportError:
        from PySide2 import QtCore, QtGui, QtWidgets

# A mid green dark enough for white text to clear 4.5:1 (lead_designer's).
ON_BG, ON_HOVER, ON_BORDER, ON_FG = "#2e7d32", "#388e3c", "#1b5e20", "#ffffff"
WARN_FG = "#c62828"


def _palette(widget):
    pal = widget.palette()
    col = lambda role: pal.color(role).name()
    return (col(QtGui.QPalette.Button), col(QtGui.QPalette.ButtonText), col(QtGui.QPalette.Mid),
            col(QtGui.QPalette.Midlight), col(QtGui.QPalette.Dark))


def toggle_style(widget, radius=3, padding="4px 10px"):
    face, text, edge, lit, off = _palette(widget)
    return ("""
        QPushButton { background: %s; color: %s; border: 1px solid %s;
                      border-radius: %dpx; padding: %s; }
        QPushButton:hover { background: %s; }
        QPushButton:disabled { color: %s; }
        QPushButton:checked { background: %s; color: %s; border: 1px solid %s; font-weight: bold; }
        QPushButton:checked:hover { background: %s; }
    """ % (face, text, edge, radius, padding, lit, off, ON_BG, ON_FG, ON_BORDER, ON_HOVER))


class StateToggle(QtWidgets.QPushButton):
    """A checkable button naming its state: off_text while off, on_text (bold,
    green) while on. set_state() changes it without firing toggled, for
    syncing to the document."""

    def __init__(self, off_text, on_text, tooltip="", parent=None):
        super().__init__(off_text, parent)
        self.off_text, self.on_text = off_text, on_text
        self.setCheckable(True)
        self.setMinimumHeight(28)
        self.setStyleSheet(toggle_style(self))
        self.setToolTip(tooltip)
        self.toggled.connect(self._relabel)

    def _relabel(self, on):
        self.setText(self.on_text if on else self.off_text)

    def set_state(self, on):
        if self.isChecked() != on:
            self.blockSignals(True)
            self.setChecked(on)
            self.blockSignals(False)
        self._relabel(on)


class Segmented(QtWidgets.QWidget):
    """Exclusive choice as joined buttons. changed(int) fires with the index."""
    changed = QtCore.Signal(int)

    def __init__(self, labels, tooltips=None, parent=None):
        super().__init__(parent)
        lay = QtWidgets.QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        self.group = QtWidgets.QButtonGroup(self)
        self.group.setExclusive(True)
        n = len(labels)
        for i, text in enumerate(labels):
            b = QtWidgets.QPushButton(text)
            b.setCheckable(True)
            left = 3 if i == 0 else 0
            right = 3 if i == n - 1 else 0
            b.setStyleSheet(toggle_style(b, radius=0, padding="3px 10px").replace(
                "border-radius: 0px;", "border-top-left-radius: %dpx; border-bottom-left-radius: %dpx; "
                "border-top-right-radius: %dpx; border-bottom-right-radius: %dpx;" % (left, left, right, right)))
            if tooltips:
                b.setToolTip(tooltips[i])
            self.group.addButton(b, i)
            lay.addWidget(b, 1)
        self.group.button(0).setChecked(True)
        self.group.idClicked.connect(self.changed.emit) if hasattr(self.group, "idClicked") \
            else self.group.buttonClicked[int].connect(self.changed.emit)

    def index(self):
        return self.group.checkedId()

    def set_index(self, i):
        self.group.button(i).setChecked(True)


def hint(text=""):
    """The panel's one grey line."""
    lab = QtWidgets.QLabel(text)
    lab.setWordWrap(True)
    lab.setStyleSheet("color: gray;")
    return lab


def swatch(hex_colour, size=12):
    lab = QtWidgets.QLabel()
    pm = QtGui.QPixmap(size, size)
    pm.fill(QtGui.QColor(hex_colour))
    lab.setPixmap(pm)
    lab.setFixedSize(size, size)
    return lab
