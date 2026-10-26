from pathlib import Path
import subprocess,time,hashlib,re,json,os
P=Path('lean');F=P/'SenConjectures.lean';S=hashlib.sha256(F.read_bytes()).hexdigest()
O=Path('core-audit');O.mkdir(exist_ok=True)
env=dict(os.environ,LEAN_NUM_THREADS='1',LEAN_PATH=str(P.resolve()))
print(subprocess.check_output(['lean','--version'],text=True),flush=True)
t=time.monotonic();r=subprocess.run(['lean','-o','SenConjectures.olean','SenConjectures.lean'],cwd=P,env=env,capture_output=True,text=True)
log=r.stdout+r.stderr;(O/'audit-SenConjectures-compile.log').write_text(log)
assert r.returncode==0,log
statements=re.findall(r"^'([^']+)' (?:depends on axioms: \[([^\]]*)\]|does not depend on any axioms)",log,re.M)
assert len(statements)==57,len(statements)
assert 'sorryAx' not in log and 'Lean.ofReduceBool' not in log
for name,ax in statements:assert set(filter(None,(x.strip() for x in ax.split(','))))<={'propext'},(name,ax)
entry={'module':'SenConjectures','lean':'4.34.1','mathlib':None,'source_sha256':S,'named_statements':57,'compile_exit_code':0,'compile_seconds':round(time.monotonic()-t,3),'compile_log_sha256':hashlib.sha256(log.encode()).hexdigest(),'statements':[{'name':n,'axioms':[x.strip() for x in a.split(',') if x.strip()]} for n,a in statements]}
print('Core compile passed; replaying kernel declarations.',flush=True)
t=time.monotonic();r=subprocess.run(['leanchecker','-v','SenConjectures'],cwd=P,env=env,capture_output=True,text=True)
log=r.stdout+r.stderr;(O/'audit-SenConjectures-replay.log').write_text(log)
assert r.returncode==0,log
assert hashlib.sha256(F.read_bytes()).hexdigest()==S
entry.update(replay_exit_code=0,replay_seconds=round(time.monotonic()-t,3),replay_log_sha256=hashlib.sha256(log.encode()).hexdigest(),olean_sha256=hashlib.sha256((P/'SenConjectures.olean').read_bytes()).hexdigest(),execution='GitHub Actions',run_url=os.environ['AUDIT_RUN_URL'])
(O/'core-audit.json').write_text(json.dumps(entry,indent=2)+'\n')
print(json.dumps({k:v for k,v in entry.items() if k!='statements'},indent=2))
