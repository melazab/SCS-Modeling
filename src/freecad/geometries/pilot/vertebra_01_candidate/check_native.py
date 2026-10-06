import FreeCAD as App, FreeCADGui as Gui, Part, MeshPart
from PySide import QtCore
from pathlib import Path
import json,time,traceback
HERE=Path('/home/mohamed/Projects/SCS-Modeling/src/freecad/geometries/pilot/vertebra_01_candidate'); start=time.perf_counter(); r={}
def shapeinfo(s):return {'valid':s.isValid(),'closed':s.isClosed(),'solids':len(s.Solids),'faces':len(s.Faces),'volume_mm3':s.Volume,'bbox':str(s.BoundBox)}
def write(): (HERE/'native_checks.json').write_text(json.dumps(r,indent=2)+'\n')
try:
 d=App.openDocument(str(HERE/'Vertebra_01_Assembly_Candidate.FCStd'))
 b=d.getObject('VertebralBody'); p=d.getObject('PosteriorFaceted'); tip=d.getObject('Channels')
 r['body']=shapeinfo(b.Shape);r['posterior']=shapeinfo(p.Shape);r['tip']=shapeinfo(tip.Shape);write()
 # Check meaningful native parameters through final body feature, then restore.
 checks=[]
 for obj,prop,delta in [(d.getObject('BodyHeight'),'Length',1),(d.getObject('EndplateRims'),'Radius',.1)]:
  original=getattr(obj,prop).Value;setattr(obj,prop,original+delta);d.recompute();checks.append({'parameter':obj.Name+'.'+prop,'delta_mm':delta,'result':shapeinfo(b.Shape)});setattr(obj,prop,original);d.recompute()
 h=d.getObject('ChannelCenters');radius=h.Constraints[0].Value;h.setDatum(0,App.Units.Quantity(str(radius+.05)+' mm'));d.recompute();checks.append({'parameter':'ChannelCenters radius constraint 0','delta_mm':.05,'result':shapeinfo(b.Shape)});h.setDatum(0,App.Units.Quantity(str(radius)+' mm'));d.recompute()
 r['parameter_checks']=checks; r['restored_body']=shapeinfo(b.Shape);write()
 # Try a direct fuzzy Boolean; preserve static status if successful.
 fused=b.Shape.fuse(p.Shape,1e-4);r['fuzzy_union']=shapeinfo(fused);r['elapsed_seconds']=time.perf_counter()-start;write()
 if len(fused.Solids)==1 and fused.isValid():
  f=d.addObject('PartDesign::Feature','DiagnosticFusion');f.Label='Hybrid fusion — static diagnostic (not linked to parameters)';f.Shape=fused;d.getObject('Vertebra_01').addObject(f);b.Visibility=False;p.Visibility=False;f.ViewObject.ShapeColor=(.82,.73,.57)
  MeshPart.meshFromShape(Shape=fused,LinearDeflection=.005,AngularDeflection=.10,Relative=False).write(str(HERE/'native_union.stl'))
  Gui.activeDocument().activeView().viewAxonometric();Gui.activeDocument().activeView().fitAll();Gui.activeDocument().activeView().saveImage(str(HERE/'candidate.png'),1280,1024,'White')
  d.recompute();d.saveAs(str(HERE/'Vertebra_01_Hybrid_Static_Candidate.FCStd'))
except Exception:r['error']=traceback.format_exc();write()
finally:
 for n in list(App.listDocuments()):App.closeDocument(n)
 QtCore.QTimer.singleShot(1000,Gui.getMainWindow().close)
