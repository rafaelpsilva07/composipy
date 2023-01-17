# This file contains the descriptors that validates composipy data
# https://www.youtube.com/watch?v=ovsvGtWD90Y&ab_channel=LivePython

import numbers


class _NumberDescriptor:
    def __init__(self, value, n_min=None, n_max=None, name=None):
        self.value = value
        self.n_min = n_min
        self.n_max = n_max
        self.name = name

    def __get__(self, obj, objtype):
        if not isinstance(self.value, numbers.Real):
            raise ValueError(f'{self.name} Must be a real number')
        if self.n_min is not None and self.value < self.n_min:
            raise ValueError(f'{self.name} must be at least {self.n_min}')
        if self.n_max is not None and self.value > self.n_max:
            raise ValueError(f'{self.name} must not be bigger than {self.n_max}')
        return self.value
    
    def __set__(self, obj, objtype):
        raise AttributeError(f'{self.name} cannot be setted')





class SimpleClass:
    x = _NumberDescriptor(10, n_min=6, n_max=10, name='x')
    y = 2


if __name__ == '__main__':
    s = SimpleClass()
    print(s.x)
    print(type(s.x))
    print(s.y)
    s.x = 9

