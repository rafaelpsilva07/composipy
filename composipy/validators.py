import numbers


class ComposipyValidator:

    def _float_number(self, value, n_min=None, n_max=None, name=None):
        try:
            value = float(value)
        except:
            raise ValueError(f'{name} Must be a real number')
        if n_min is not None and value < n_min:
            raise ValueError(f'{name} must be at least {n_min}')
        if n_max is not None and value > n_max:
            raise ValueError(f'{name} cannot exceed {n_max}')
        return value


    def _int_number(self, value, n_min=None, n_max=None, name=None):
        try:
            value = float(value)
        except:
            raise ValueError(f'{name} Must be a integer number')        
        if not float.is_integer(value):
            raise ValueError(f'{name} Must be a integer number')        
        if n_min is not None and value < n_min:
            raise ValueError(f'{name} must be at least {n_min}')
        if n_max is not None and value > n_max:
            raise ValueError(f'{name} cannot exceed {n_max}')
        return int(value)
    
    