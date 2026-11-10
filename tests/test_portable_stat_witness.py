"""Actual metadata helpers: cross-kernel dev only, never other file drift."""
import ast,json,pathlib,types,unittest
from unittest.mock import patch
ROOT=pathlib.Path(__file__).resolve().parents[1]
def helpers(name):
 source=(ROOT/'src/utils/bash'/name).read_text().split("<<'PY'\n",1)[1].split('\nPY\n}',1)[0]
 tree=ast.parse(source);ns=dict(pathlib=pathlib,json=json,root=pathlib.Path('/selected/_ecoda_runs/current'),observed_stats={})
 for node in tree.body:
  if isinstance(node,ast.FunctionDef) and node.name in ('need','stat','selected_beegfs','portable_stat'):
   exec(compile(ast.Module(body=[node],type_ignores=[]),'<actual-portable-witness>','exec'),ns)
 return ns
class Tests(unittest.TestCase):
 def test_both_actual_helpers(self):
  for script in ('legacy_analysis_audit.sh','saved_cluster_screen.sh'):
   with self.subTest(script=script):
    ns=helpers(script);prefix='' if script.startswith('legacy') else 'st_';dev=prefix+'dev'
    expected={prefix+k:v for k,v in dict(dev=68,ino=100,size=20,mode=0o400,mtime_ns=10,ctime_ns=11).items()}
    actual=dict(expected);actual[dev]=173
    good='10 1 0:68 / /selected rw - beegfs beegfs_nodev rw\n'
    with patch.object(pathlib.Path,'read_text',return_value=good):
     self.assertTrue(ns['portable_stat'](actual,expected,'/selected/input'))
     for key in expected:
      if key==dev:continue
      bad=dict(actual);bad[key]+=1
      with self.assertRaisesRegex(SystemExit,'path=/selected/input.*'+key):ns['portable_stat'](bad,expected,'/selected/input')
     with self.assertRaisesRegex(SystemExit,'escapes scratch'):ns['portable_stat'](actual,expected,'/foreign/input')
    for mount in ['10 1 0:68 / /selected rw - ext4 other rw\n',good+'11 10 0:70 / /selected/input rw - tmpfs shadow rw\n']:
     with patch.object(pathlib.Path,'read_text',return_value=mount):
      with self.assertRaisesRegex(SystemExit,'expected BeeGFS'):ns['portable_stat'](actual,expected,'/selected/input')
 def test_device_must_still_match_within_one_invocation(self):
  for script in ('legacy_analysis_audit.sh','saved_cluster_screen.sh'):
   ns=helpers(script);values=[types.SimpleNamespace(st_dev=d,st_ino=100,st_size=20,st_mode=0o100400,st_mtime_ns=10,st_ctime_ns=11) for d in (68,69)]
   class File:
    def stat(self):return values.pop(0)
    def __str__(self):return '/selected/input'
   ns['regular']=lambda path:File();ns['stat']('/selected/input')
   with self.assertRaisesRegex(SystemExit,'in-process stat changed.*dev'):ns['stat']('/selected/input')
if __name__=='__main__':unittest.main()
