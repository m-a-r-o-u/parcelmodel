from unittest import TestCase, skip
from config_tree  import replace_string_in_tree
from config_tree  import set_dict_item_in_tree
from config_tree  import modify_item_in_config_tree

class TestTreeModify(TestCase):
    def test_replace_part_of_string_in_dict(self):
        inp = {'a':'Heeeyyy', 'b':'a/time_stamp/'}
        out = {'a':'Heeeyyy', 'b':'a/2017/'}
        self.assertEqual(replace_string_in_tree(inp, 'time_stamp', '2017'), out)

    def test_replace_part_of_string_in_dict_with_different_result(self):
        inp = {'a':'Heeeyjo', 'b':'a/time_stamp/'}
        out = {'a':'Heeeyjo', 'b':'a/2016/'}
        self.assertEqual(replace_string_in_tree(inp, 'time_stamp', '2016'), out)

    def test_replace_part_of_string_in_dict_containing_list(self):
        inp = {'a':[1,2,3], 'b':'a/time_stamp/'}
        out = {'a':[1,2,3], 'b':'a/2016/'}
        self.assertEqual(replace_string_in_tree(inp, 'time_stamp', '2016'), out)

    def test_replace_part_of_string_in_dict_inside_of_list(self):
        inp = [{'a':[1,2,3], 'b':'a/time_stamp/'}]
        out = [{'a':[1,2,3], 'b':'a/2016/'}]
        self.assertEqual(replace_string_in_tree(inp, 'time_stamp', '2016'), out)

    def test_replace_string_in_list(self):
        inp = [{'a':[1,2,3], 'b':'a/time_stamp/'}, 'time_stamp']
        out = [{'a':[1,2,3], 'b':'a/2016/'}, '2016']
        self.assertEqual(replace_string_in_tree(inp, 'time_stamp', '2016'), out)

    def test_replace_string_in_dict_inside_of_dicts(self):
        inp = {'a':{'b':{'c':'time_stamp'}}}
        out = {'a':{'b':{'c':'2017'}}}
        self.assertEqual(replace_string_in_tree(inp, 'time_stamp', '2017'), out)

    def test_replace_part_of_string_in_list(self):
        inp = [{'a':[1,2,3], 'b':'a/time_stamp/'}, 'b/time_stamp']
        out = [{'a':[1,2,3], 'b':'a/2016/'}, 'b/2016']
        self.assertEqual(replace_string_in_tree(inp, 'time_stamp', '2016'), out)

    def test_replace_string_in_list_inside_of_lists(self):
        inp = ['a',['b',['c','time_stamp']]]
        out = ['a',['b',['c','2017']]]
        self.assertEqual(replace_string_in_tree(inp, 'time_stamp', '2017'), out)

    def test_replace_part_of_string_in_list_inside_dict(self):
        inp = {'b':['a/time_stamp/']}
        out = {'b':['a/2016/']}
        self.assertEqual(replace_string_in_tree(inp, 'time_stamp', '2016'), out)

class TestSetItemInTree(TestCase):
    def test_set_global_path(self):
        inp = {'output':{'file_path':'local_path'}}
        out = {'output':{'file_path':'global_path'}}
        self.assertEqual(set_dict_item_in_tree(inp, 'file_path', 'global_path', ['output']),out)
        self.assertEqual(inp, {'output':{'file_path':'local_path'}})

@skip
class TestSetGlobalsInConfigTree(TestCase):
    def test_replace_dict_item(self):
        inp = {'a':1, 'b':2, 'key':'local_path'}
        out = {'a':1, 'b':2, 'key':'global_path'}
        self.assertEqual(modify_item_in_config_tree(inp, 'key', 'global_path'), out)

    def test_replace_dictitem_in_list_with_input(self):
        inp = [{'a':1, 'b':2, 'key':'local_path'}]
        out = [{'a':1, 'b':2, 'key':'global_path'}]
        self.assertEqual(modify_item_in_config_tree(inp, 'key', 'global_path'), out)

    def test_replace_dictitem_in_list_of_dicts_with_input(self):
        inp = [{'a':1, 'b':2, 'key':'local'}, {'a':1, 'b':2, 'key':'local'}]
        out = [{'a':1, 'b':2, 'key':'global'}, {'a':1, 'b':2, 'key':'global'}]
        self.assertEqual(modify_item_in_config_tree(inp, 'key', 'global'), out)

    def test_replace_dict_item(self):
        inp = {'a':1, 'b':2}
        out = {'a':1, 'b':2, 'key':'global_path'}
        self.assertEqual(modify_item_in_config_tree(inp, 'key', 'global_path'), out)
