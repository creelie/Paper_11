"""Compile and replay the three local proof modules, with exact scope receipts.

This script checks source code, not the unformalised mathematical claims in
the paper. Imported Mathlib declarations are trusted dependencies of replay.
"""
from pathlib import Path
import os, subprocess, json, re, hashlib, time, datetime, sys

root = Path(__file__).resolve().parent
project = root/'revision/lean'
packages = Path(r'D:\CLAUDE PAPERS\research-closure-2026-10-23\NavierStokesAndEuler\.lake\packages')
mathlib = '85e3a25e006c35636f0e53b0e9296caca2685bc0'
configurations = [
    ('SenConjectures', '4.34.1', None, 57),
    ('AnalyticScopeChecks', '4.34.0-rc2', mathlib, 12),
    ('InfraredBounds', '4.34.0-rc2', mathlib, 6),
]
supplements_only = '--supplements' in sys.argv
if supplements_only:
    configurations = configurations[1:]
report = {'status': 'running', 'modules': [], 'full_manuscript_formalised': False,
          'full_four_dimensional_infrared_closure': False,
          'replay_scope': "Lean's own kernel replay, relative to module imports; not a separate proof assistant."}
receipt = root/('supplement-audit.json' if supplements_only else 'all-lean-audit.json')


def save():
    receipt.write_text(json.dumps(report, indent=2)+'\n', encoding='utf8')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


save()
for module, version, dependency, expected_count in configurations:
    source = project/(module+'.lean')
    source_hash = sha(source)
    binary = Path(r'C:\Users\user\.elan\toolchains')/('leanprover--lean4---v'+version)/'bin'
    environment = os.environ.copy()
    environment['LEAN_NUM_THREADS'] = '1'
    paths = [str(project)]
    if dependency:
        paths += [str(p/'.lake/build/lib/lean') for p in packages.iterdir() if p.is_dir()]
    environment['LEAN_PATH'] = os.pathsep.join(paths)
    text = source.read_text(encoding='utf8')
    text_without_comments = re.sub(r'/\-.*?\-/', '', text, flags=re.S)
    text_without_comments = re.sub(r'--[^\n]*', '', text_without_comments)
    assert not re.search(r'\b(sorry|admit|axiom|native_decide|ofReduceBool)\b', text_without_comments)
    assert len(re.findall(r'^#print axioms ', text, re.M)) == expected_count
    print('COMPILING '+module+' with Lean '+version, flush=True)
    start = time.monotonic()
    compiled = subprocess.run([str(binary/'lean.exe'), '-o', module+'.olean', module+'.lean'],
                              cwd=project, env=environment, capture_output=True, encoding='utf8',
                              errors='replace', timeout=3600)
    log = compiled.stdout+compiled.stderr
    compile_log = root/('audit-'+module+'-compile.log')
    compile_log.write_text(log, encoding='utf8')
    assert compiled.returncode == 0, log
    assert sha(source) == source_hash, 'Source changed during compilation'
    assert not any(token in log for token in ['sorryAx', 'Lean.ofReduceBool'])
    statements = re.findall(r"^'([^']+)' (?:depends on axioms: \[([^\]]*)\]|does not depend on any axioms)", log, re.M)
    assert len(statements) == expected_count, (module, len(statements))
    allowed = {'propext'} if dependency is None else {'propext', 'Classical.choice', 'Quot.sound'}
    for name, axioms in statements:
        used = set(filter(None, (x.strip() for x in axioms.split(','))))
        assert used <= allowed, (name, used-allowed)
    entry = {'module': module, 'lean': version, 'mathlib': dependency,
             'source_sha256': source_hash, 'named_statements': expected_count,
             'compile_exit_code': compiled.returncode,
             'compile_seconds': round(time.monotonic()-start, 3),
             'compile_log_sha256': sha(compile_log),
             'statements': [{'name': n, 'axioms': [a.strip() for a in ax.split(',') if a.strip()]}
                            for n, ax in statements]}
    report['modules'].append(entry)
    save()
    print('REPLAYING '+module, flush=True)
    start = time.monotonic()
    replay = subprocess.run([str(binary/'leanchecker.exe'), '-v', module], cwd=project,
                            env=environment, capture_output=True, encoding='utf8',
                            errors='replace', timeout=3600)
    replay_log = root/('audit-'+module+'-replay.log')
    replay_log.write_text(replay.stdout+replay.stderr, encoding='utf8')
    assert replay.returncode == 0, replay.stdout+replay.stderr
    entry.update(replay_exit_code=replay.returncode,
                 replay_seconds=round(time.monotonic()-start, 3),
                 replay_log_sha256=sha(replay_log), olean_sha256=sha(project/(module+'.olean')))
    save()
    print('PASSED '+module+' ('+str(expected_count)+' named statements)', flush=True)
report.update(status='passed', named_statements=sum(x['named_statements'] for x in report['modules']),
              completed_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
save()
print(json.dumps({k:v for k,v in report.items() if k != 'modules'}, indent=2), flush=True)
