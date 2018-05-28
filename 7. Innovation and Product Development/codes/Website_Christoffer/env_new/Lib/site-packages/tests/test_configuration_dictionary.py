from nose.tools import raises
#from tg.configuration.utils import ConfigurationBunch


class TestConfigurationDictionary(object):
    def setup(self):
        self.d = ConfigurationBunch()

"""
    def test_plain_to_nested(self):
        d = self.d
        d['1st.2nd.3rd'] = 'HELLO'
        assert d['1st']['2nd']['3rd'] == 'HELLO'

    def test_nested_to_plain(self):
        d = self.d
        d['1st'] = {}
        d['1st']['2nd'] = {}
        d['1st']['2nd']['3rd'] = 'HELLO'
        assert d['1st.2nd.3rd'] == 'HELLO'

    def test_plain_to_plain(self):
        d = self.d
        d['1st.2nd.3rd'] = 'HELLO'
        assert d['1st.2nd.3rd'] == 'HELLO'

    def test_nested_to_attrs(self):
        d = self.d
        d['st1'] = {}
        d['st1']['nd2'] = {}
        d['st1']['nd2']['rd3'] = 'HELLO'
        assert d.st1.nd2.rd3 == 'HELLO'

    def test_plain_to_attrs(self):
        d = self.d
        d['st1.nd2.rd3'] = 'HELLO'
        assert d.st1.nd2.rd3 == 'HELLO'

    def test_iterate(self):
        d = self.d
        d['simple.sub'] = {'sub1': [3, 4], 'sub2': 'hi', 'sub3': {'subsub1': 7, 'subsub2': 8}}
        d['simple.sub.x'] = 'hi'

        flatkeys = list(iter(d))
        assert flatkeys == ['simple.sub.x', 'simple.sub.sub2', 'simple.sub.sub3.subsub2',
                            'simple.sub.sub3.subsub1', 'simple.sub.sub1'], flatkeys

    @raises(KeyError)
    def test_delete_plain(self):
        d = self.d
        d['st1'] = {}
        d['st1']['nd2'] = {}
        d['st1']['nd2']['rd3'] = 'HELLO'
        d.pop('st1.nd2.rd3', None)

        d['st1']['nd2']['rd3']

    @raises(AttributeError)
    def test_plain_to_attrs_not_found(self):
        d = self.d
        d['st1.nd2.rd3'] = 'HELLO'
        d.st1.nd2.rd4

    def test_set_plain(self):
        self.d['hi'] = 5
        self.d['simple.sub.sub1'] = [3, 4]
        self.d['simple.sub.sub2'] = [3, 4]
        assert sorted(list(self.d.keys())) == sorted(['hi', 'simple']), self.d
        assert 'sub' in self.d['simple'], self.d
        assert len(self.d['simple.sub'].keys()) == 2

    def test_set_subdict(self):
        d = self.d
        d['simple.sub'] = {'sub1': [3, 4], 'sub2': 'hi', 'sub3': {'subsub1': 7, 'subsub2': 8}}
        d['simple.sub.x'] = 'hi'
        assert len(d) == 1, d
        assert d['simple'] == {'sub': {'x': 'hi', 'sub1': [3, 4], 'sub2': 'hi', 'sub3': {'subsub1': 7, 'subsub2': 8}}}

    def test_get_subdict(self):
        self.test_set_subdict()
        assert self.d['simple.sub.sub3'] == {'subsub1': 7, 'subsub2': 8}

    @raises(TypeError)
    def test_odd_parents_and_children1(self):
        self.d['parent'] = 5
        self.d['parent.child'] = 3

    @raises(TypeError)
    def test_odd_parents_and_children2(self):
        self.d['parent.child'] = 3
        self.d['parent'] = 3

        # parent is now an int so it should fail
        self.d['parent.child']
"""