from pathlib import Path
import os, subprocess, time, json, hashlib

root = Path(__file__).resolve().parent
project = root / 'revision/lean'
packages = Path(r'D:\CLAUDE PAPERS\research-closure-2026-10-23\NavierStokesAndEuler\.lake\packages')
environment = dict(os.environ, LEAN_PATH=os.pathsep.join(
    [str(project)] + [str(p/'.lake/build/lib/lean') for p in packages.iterdir() if p.is_dir()]))
executable = r'C:\Users\user\.elan\toolchains\leanprover--lean4---v4.34.0-rc2\bin\lean.exe'
start = time.monotonic()
result = subprocess.run([executable, '-o', 'InfraredBounds.olean', 'InfraredBounds.lean'],
                        cwd=project, env=environment, capture_output=True, encoding='utf8')
log = result.stdout + result.stderr
(root/'infrared-lean.log').write_text(log, encoding='utf8')
receipt = {'exit_code': result.returncode, 'wall_seconds': round(time.monotonic()-start, 3),
           'lean': '4.34.0-rc2', 'mathlib': '85e3a25e006c35636f0e53b0e9296caca2685bc0',
           'source_sha256': hashlib.sha256((project/'InfraredBounds.lean').read_bytes()).hexdigest(),
           'scope': 'Six elementary complex-norm, corner-algebra and degree-count theorems; no amplitude or integral formalisation.'}
(root/'infrared-run.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf8')
print(log)
print(json.dumps(receipt))
raise SystemExit(result.returncode)
