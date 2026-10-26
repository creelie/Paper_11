from pathlib import Path
import os,subprocess,time,json
R=Path(__file__).resolve().parent;P=R/'revision/lean'
PK=Path(r'D:\CLAUDE PAPERS\research-closure-2026-10-23\NavierStokesAndEuler\.lake\packages')
env=dict(os.environ,LEAN_PATH=os.pathsep.join([str(P)]+[str(p/'.lake/build/lib/lean') for p in PK.iterdir() if p.is_dir()]))
exe=r'C:\Users\user\.elan\toolchains\leanprover--lean4---v4.34.0-rc2\bin\lean.exe'
t=time.monotonic();p=subprocess.run([exe,'-o','AnalyticScopeChecks.olean','AnalyticScopeChecks.lean'],cwd=P,env=env,capture_output=True,encoding='utf8')
(R/'analytic-lean.log').write_text(p.stdout+p.stderr,encoding='utf8');print(p.stdout+p.stderr)
(R/'analytic-run.json').write_text(json.dumps({'exit_code':p.returncode,'wall_seconds':round(time.monotonic()-t,3),'lean':'4.34.0-rc2','mathlib':'85e3a25e006c35636f0e53b0e9296caca2685bc0'},indent=2))
raise SystemExit(p.returncode)
