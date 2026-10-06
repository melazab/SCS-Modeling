"""Engineering construction in an isolated FreeCAD GUI; not a lesson macro."""
import json, math, time, traceback
from pathlib import Path
import FreeCAD as App, FreeCADGui as Gui, Part, Sketcher, Mesh, MeshPart
from PySide import QtCore
HERE=Path('/home/mohamed/Projects/SCS-Modeling/src/freecad/geometries/pilot/vertebra_01_candidate')
BASE=HERE.parent/'vertebra_01'
started=time.perf_counter()
(HERE/'build_error.txt').unlink(missing_ok=True)
def mark(s):
    with (HERE/'build_progress.log').open('a') as f: f.write(f'{time.perf_counter()-started:.3f} {s}\n')
def run():
    data=json.loads((BASE/'measurements.json').read_text())
    doc=App.newDocument('Vertebra_01_Candidate')
    ass=doc.addObject('App::Part','Vertebra_01'); ass.Label='Vertebra 01 — hybrid candidate; level unassigned'
    for name,value in {'SourceFilename':data['source_filename'],'SourceSHA256':data['source_sha256'],'Status':'Hybrid: editable body; faceted posterior retained; accuracy pending'}.items():
        ass.addProperty('App::PropertyString',name,'Provenance'); setattr(ass,name,value)
    body=doc.addObject('PartDesign::Body','VertebralBody'); ass.addObject(body)
    body.Label='Vertebral body — fitted editable features'
    cx=sum(p['center'][0] for p in data['channels'])/2; cy=sum(p['center'][1] for p in data['channels'])/2
    z=data['bottom_local_z_mm']; R=data['rotation_global_to_local']
    origin=App.Vector(*[R[0][i]*cx+R[1][i]*cy+R[2][i]*z for i in range(3)])
    angle=-math.degrees(math.atan2(R[2][1],R[2][2]))
    body.Placement=App.Placement(origin,App.Rotation(App.Vector(1,0,0),angle))
    sketch=body.newObject('Sketcher::SketchObject','BodyOutline'); sketch.Label='01 — Body outline: fitted periodic spline'
    curve=Part.BSplineCurve(); knots=data['spline_knots'][3:-3]
    curve.buildFromPolesMultsKnots([App.Vector(x-cx,y-cy,0) for x,y in data['spline_coefficients_xy'][:-3]],[1]*len(knots),knots,True,3)
    sketch.addGeometry(curve,False)
    pad=body.newObject('PartDesign::Pad','BodyHeight'); pad.Profile=sketch; pad.Length=data['height_mm']; pad.Label='02 — Body height'
    doc.recompute(); mark('pad built')
    # Identify the two closed boundary edges of the extrusion, not its seam.
    edgeids=[f'Edge{i+1}' for i,e in enumerate(pad.Shape.Edges) if e.isClosed()]
    fillet=body.newObject('PartDesign::Fillet','EndplateRims'); fillet.Base=(pad,edgeids); fillet.Radius=1.5; fillet.Label='03 — Rounded endplate rims'
    doc.recompute(); assert fillet.Shape.isValid() and len(fillet.Shape.Solids)==1
    mark('fillet built '+str(edgeids))
    holes=body.newObject('Sketcher::SketchObject','ChannelCenters'); holes.Label='04 — Two measured channels'
    for i,p in enumerate(data['channels']):
        x,y=p['center']; j=holes.addGeometry(Part.Circle(App.Vector(x-cx,y-cy,0),App.Vector(0,0,1),p['radius']),False)
        holes.addConstraint(Sketcher.Constraint('Radius',j,p['radius']))
        holes.addConstraint(Sketcher.Constraint('DistanceX',j,3,x-cx))
        # Keep measured tiny transverse asymmetry rather than exactly symmetric assumption.
        holes.addConstraint(Sketcher.Constraint('DistanceY',j,3,y-cy))
    pocket=body.newObject('PartDesign::Pocket','Channels'); pocket.Profile=holes; pocket.Type=1; pocket.Reversed=True; pocket.Label='05 — Through channels'
    doc.recompute(); assert pocket.Shape.isValid() and len(pocket.Shape.Solids)==1
    mark('pocket built')
    for o in (sketch,pad,fillet,holes): o.Visibility=False
    pocket.Visibility=True; pocket.ViewObject.ShapeColor=(0.88,0.7,0.42)
    result={'body_valid':pocket.Shape.isValid(),'body_closed':pocket.Shape.isClosed(),'body_volume_mm3':pocket.Shape.Volume,'body_faces':len(pocket.Shape.Faces),'fillet_edges':edgeids,'pocket_status':pocket.State,'fillet_status':fillet.State}
    # MeshPart explicit absolute deflection gives repeatable comparison tessellation.
    shape=pocket.Shape.copy(); shape.Placement=body.Placement.multiply(shape.Placement)
    mesh=MeshPart.meshFromShape(Shape=shape,LinearDeflection=.005,AngularDeflection=.10,Relative=False)
    mesh.write(str(HERE/'body_candidate.stl')); shape.exportBrep(str(HERE/'body_candidate.brep'))
    mark('body exported')
    # Save intermediate document before attempting expensive posterior conversion.
    doc.recompute(); doc.saveAs(str(HERE/'Vertebra_01_Body_Candidate.FCStd'))
    raw=json.loads((HERE/'posterior_working.json').read_text())
    posterior=Part.Shape(); posterior.makeShapeFromMesh(([App.Vector(*p) for p in raw['vertices']],[tuple(f) for f in raw['faces']]),1e-7)
    mark('posterior shell built')
    posterior=Part.makeSolid(posterior)
    mark('posterior validity '+str((posterior.isValid(),posterior.isClosed(),posterior.Volume)))
    assert posterior.isValid() and posterior.isClosed() and posterior.Volume>0
    rear=doc.addObject('PartDesign::Feature','PosteriorFaceted'); ass.addObject(rear); rear.Shape=posterior; rear.Label='Posterior — faceted reference solid (NOT parametrized)'; rear.ViewObject.ShapeColor=(0.60,0.68,0.77)
    result.update(posterior_valid=posterior.isValid(),posterior_closed=posterior.isClosed(),posterior_faces=len(posterior.Faces),posterior_volume_mm3=posterior.Volume)
    mark('posterior solid validated')
    doc.recompute(); doc.saveAs(str(HERE/'Vertebra_01_Assembly_Candidate.FCStd'))
    # Boolean may be costly: report actual duration and retain editable operands.
    union=doc.addObject('Part::MultiFuse','CompleteHybrid'); ass.addObject(union); union.Shapes=[body,rear]; union.Refine=True; union.Label='Complete vertebra — hybrid union'
    doc.recompute(); mark('union recomputed')
    result.update(union_valid=union.Shape.isValid(),union_closed=union.Shape.isClosed(),union_solids=len(union.Shape.Solids),union_volume_mm3=union.Shape.Volume,union_faces=len(union.Shape.Faces),union_status=union.State)
    mark('union metrics '+str(result))
    (HERE/'construction_result.json').write_text(json.dumps(result,indent=2)+'\n')
    assert result['union_valid'] and result['union_solids']==1
    mesh=MeshPart.meshFromShape(Shape=union.Shape,LinearDeflection=.005,AngularDeflection=.10,Relative=False); mesh.write(str(HERE/'hybrid_candidate.stl'))
    union.Shape.exportBrep(str(HERE/'hybrid_candidate.brep'))
    body.Visibility=False; rear.Visibility=False; union.Visibility=True; union.ViewObject.ShapeColor=(.82,.73,.57)
    ass.Status='Hybrid candidate: editable body, faceted posterior; see measured report'
    Gui.activeDocument().activeView().viewAxonometric(); Gui.activeDocument().activeView().fitAll()
    Gui.activeDocument().activeView().saveImage(str(HERE/'candidate.png'),1280,1024,'White')
    doc.recompute(); doc.saveAs(str(HERE/'Vertebra_01_Hybrid_Candidate.FCStd'))
    result['construction_wall_seconds']=time.perf_counter()-started
    (HERE/'construction_result.json').write_text(json.dumps(result,indent=2)+'\n'); mark('done')
try: run()
except Exception:
    (HERE/'build_error.txt').write_text(traceback.format_exc()); mark('ERROR')
finally:
    for name in list(App.listDocuments()): App.closeDocument(name)
    QtCore.QTimer.singleShot(1000,Gui.getMainWindow().close)
