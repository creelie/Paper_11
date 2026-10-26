from pathlib import Path
import hashlib,json,re
P=Path(__file__).resolve().parents[1]
manifest=json.loads((P/'PACKAGE-MANIFEST.json').read_text(encoding='utf8'))
for item in manifest:
 p=P/item['path'];assert p.is_file(),item['path']
 assert p.stat().st_size==item['bytes'],item['path']
 assert hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256'],item['path']
receipt=json.loads((P/'verification/all-lean-audit.json').read_text(encoding='utf8'))
assert receipt['status']=='passed' and receipt['named_statements']==75
assert not receipt['full_manuscript_formalised'] and not receipt['full_four_dimensional_infrared_closure']
for m in receipt['modules']:
 assert m['compile_exit_code']==m['replay_exit_code']==0
 assert hashlib.sha256((P/'lean'/(m['module']+'.lean')).read_bytes()).hexdigest()==m['source_sha256']
main=(P/'manuscript/main.tex').read_text(encoding='utf8')
assert 'onecolumn' in main and 'twocolumn' not in main
assert not any(x in main for x in ['\\usepackage{lineno}', '\\linenumbers'])
assert not list((P/'manuscript').rglob('*.lean'))
print('Frozen payload hashes, scoped 75-statement receipt and single-column source controls passed.')
