"""Preparation errors must be visible in both Mecway UI paths."""
import ast
from pathlib import Path
import unittest

class NoticeTests(unittest.TestCase):
    def test_mesh_failure_shows_reason_once_in_both_ui_modes(self):
        source=Path(__file__).resolve().parents[1]/'Scripts'/'MecwayAsterBridgeMain.py'
        tree=ast.parse(source.read_text(encoding='utf-8-sig'))
        funcs=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ('_notify','_notify_mesh_failure')]
        for popup in (False,True):
            messages=[]
            mw=type('Mecway',(),{'message':staticmethod(messages.append)})()
            scope={'mw':mw,'_NOTES':{'list':[],'popup':popup},'tr':lambda nl,en:en,'log_file_name':lambda:'bridge.log'}
            exec(compile(ast.Module(body=funcs,type_ignores=[]),str(source),'exec'),scope)
            scope['_notify_mesh_failure']({'preparation_error':'Beam offsets: disable modal prestress and complex modes for offset beams.'})
            self.assertEqual(len(messages),1)
            self.assertIn('disable modal prestress',messages[0])
            self.assertIn('was not started',messages[0])
            self.assertEqual(scope['_NOTES']['list'],messages)

if __name__=='__main__': unittest.main()
