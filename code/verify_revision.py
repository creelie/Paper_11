from pathlib import Path
import json,hashlib,importlib.util
R=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('dyons',R/'siegel_dyons.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def cone(q,p,l):return q>0 and p>0 and q*p>l*l and 3*l<=2*q+p and l<=2*q and l<=p and l>=0
stored=json.loads((R/'siegel_dyons_out.json').read_text())
large=m.coefficient_tables('tor2',2,14,5,6,'R')
out=[]
for q in range(13):
 for p in range(0,9,2):
  for l in range(5):
   if cone(q,p,l):
    d=m.d_index(large,2,q,p//2,l)
    if d<0:out.append([q,p,l,d])
assert out==stored['tor2_negative_single_centred'];assert len(out)==30
assert cone(3,2,2) and not (3*2*2<=2*3+2),'Linear/squared predicate control failed'
assert all(x[3]<0 for x in out)
report={'negative_triples':len(out),'truncation_stability':'larger L=14,K=5,P=6 reproduces entire scan box','linear_cone_negative_control':'passed','full_Lean_coefficient_coverage':False,'scope':'exact finite product computation; no physical horizon identification or proof of all-charge positivity'}
(R/'revision-checks.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
