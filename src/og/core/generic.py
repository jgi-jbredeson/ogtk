
class Enumerative(object):
    def __init__(self, iterable=(), *args, **kwargs):
        self._enumerative_values = set()
        self._set_attrs(iterable)
        self._set_attrs(args)
        self._set_attrs(kwargs)

    def _set_attrs(self, values):
        if isinstance(values, dict):
            for attr in values:
                self.add(attr, values[attr])
        else:
            for attr in values:
                self.add(attr)
            
    def __contains__(self, value):
        return value in self._enumerative_values

    def add(self, key, value=None):
        if value is None:
            value = key
        if not hasattr(self, key):
            self._enumerative_values.add(key)
            self._enumerative_values.add(value)
            setattr(self, key, value)
            
    


    
