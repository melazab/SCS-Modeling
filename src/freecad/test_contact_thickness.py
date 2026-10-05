"""Contact-wall checks. Call main() in FreeCAD; creates no document objects."""
import math
import build_lead_config as blc


def main():
    import make_scs_lead as msl
    cfg=blc.parse_defaults_file()
    assert blc.resolve_lead(cfg,{'type':'dorsal'})['contact_thickness']==0
    for value in (-.1,float('nan'),float('inf')):
        try:
            blc.resolve_lead(cfg,{'type':'dorsal','contact_thickness':value})
        except ValueError:
            pass
        else:
            raise AssertionError(f'Invalid thickness accepted: {value}')
    try:
        blc.resolve_lead(cfg,{'type':'drg','contact_thickness':.1})
    except ValueError:
        pass
    else:
        raise AssertionError('Fixed DRG hardware accepted a thickness override')
    segments,_,_=msl.segment_plan(2,3,1,1,0)
    old=msl.sweep_path(segments,lambda z:(0,0,z),.65)
    same=msl.sweep_path(segments,lambda z:(0,0,z),.65,0)
    raised=msl.sweep_path(segments,lambda z:(0,0,z),.65,.1)
    for (_,_,a),(_,_,b) in zip(old,same):
        assert a.exportBrepToString()==b.exportBrepToString()
    assert msl.fuse_insulator(raised).isValid()
    for k,i,s in raised:
        assert s.isValid()
        if k=='contact':
            assert abs(s.Volume-math.pi*.75**2*3)<1e-9
            assert abs(2*s.BoundBox.XMax-1.5)<1e-9
    expected_extra=2*math.pi*(.75**2-.65**2)*3
    assert abs(sum(s.Volume for k,i,s in raised)-sum(s.Volume for k,i,s in old)-expected_extra)<1e-9
    print('PASS: projection validation, legacy geometry, contact diameter and volume')


if __name__=='__main__':main()
