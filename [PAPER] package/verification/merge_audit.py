from pathlib import Path
import json, hashlib, datetime

root = Path(__file__).resolve().parent
supplement = json.loads((root/'supplement-audit.json').read_text())
core = json.loads((root/'core-audit.json').read_text())
assert supplement['status'] == 'passed' and supplement['named_statements'] == 18
assert core['compile_exit_code'] == core['replay_exit_code'] == 0
modules = [core] + supplement['modules']
for module in modules:
    source = root/'revision/lean'/(module['module']+'.lean')
    assert hashlib.sha256(source.read_bytes()).hexdigest() == module['source_sha256']
    for operation in ['compile', 'replay']:
        log = root/('audit-'+module['module']+'-'+operation+'.log')
        assert hashlib.sha256(log.read_bytes()).hexdigest() == module[operation+'_log_sha256']
        assert module[operation+'_exit_code'] == 0
report = {'status': 'passed', 'modules': modules, 'named_statements': 75,
          'full_manuscript_formalised': False,
          'full_four_dimensional_infrared_closure': False,
          'replay_scope': "Lean's own kernel replay, relative to imports; not a separate proof assistant.",
          'execution': 'Core on GitHub Actions; Mathlib supplements locally. Exact matching source hashes required.',
          'local_duplicate': 'Interrupted after excessive memory paging; its result is not counted as a pass.',
          'completed_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat()}
(root/'all-lean-audit.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf8')
print('All three Lean modules compiled and kernel-replayed: 75 named statements. Full infrared closure: false.')
