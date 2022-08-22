
__all__ = ['build_sequence', 'build_pcomp']

def _build_pm_unit(angle, times):
    '''Returns a list when there is the plus minus sign '''
    stacking_values = list()
    
    for i in range(times):
        current_angle = (-1)**i*angle
        stacking_values += [current_angle]
    return stacking_values


def _build_stack_unit(stacking):
    '''Returns a list when given a unit for example \\pm45_13'''
    stacking_values = stacking
    is_pm = stacking_values.find('\pm') > -1
    stacking_values = stacking_values.replace('\pm', '')
    stacking_values = stacking_values.strip()   
    stacking_values = stacking_values.split('_')
    len_stacking_values = len(stacking_values)
    is_unit = (len_stacking_values == 2
              or len_stacking_values == 1)

    if not is_unit:
        raise ValueError(f'{stacking} it is not a unit')
    angle = float(stacking_values[0])
    if len_stacking_values == 2:
        times = stacking_values[1].replace('{', '').replace('}', '')
        times = int(times)
    elif len_stacking_values == 1:
        times = 1
    
    if is_pm:
        return _build_pm_unit(angle, times)
    else:
        return [angle] * times
   
    
def build_sequence(stacking):
    '''
    Parameters
    ----------
    stacking : str
        A string in Latex format
    
    Returns
    -------
    unities : list
        A list that contains the stacking sequence
        
    Example
    -------
    >>> stack = '[\pm45_{1}/0_2/90/(0/90)_1]s'
    >>> build_sequence(stack)
    out : 
        [45.0, -45.0, 0.0, 90.0, 0.0, 90.0, 90.0, 0.0, 90.0, 0.0, -45.0, 45.0]
    '''
    stacking = stacking.strip()
    is_symmetric = stacking[-1].lower() == 's'
    stacking = stacking.replace('s', '')
    is_stacking = (stacking[0] == '['
                   and stacking[-1] == ']')
    stacking = stacking.replace('[', '').replace(']', '')
    stacking = stacking + '/'
    
    unities = list()
    unity = list()
    current_block = ''
    is_in_block = False
    is_multiplier = False
    for i, char in enumerate(stacking):
        if char == '(':
            is_in_block = True
        elif char == ')':
            #current_block = ''
            unity += _build_stack_unit(current_block)
            current_block = ''
            is_in_block = False
            is_multiplier = True
        elif (char == '/' 
              and not is_in_block
              and not is_multiplier
             ):
            unity = _build_stack_unit(current_block)
            unities.extend(unity)
            current_block = ''
            unity = list()
        elif (char == '/' 
              and is_in_block):
            unity = _build_stack_unit(current_block)
            current_block = ''
        elif (char == '/' 
             and is_multiplier):
            times = current_block.replace('_', '').replace('{', '').replace('}', '')
            times = int(times)
            unity = times * unity
            unities.extend(unity)
            unity = list()
            current_block = ''
            is_multiplier = False           
        else:
            current_block += char
           
    if is_symmetric:
        symmetric_unities = unities[::-1]
        unities.extend(symmetric_unities)
    return unities


def build_pcomp(sequence):
    '''
    Parameters
    ----------
    sequence : list
        A list of angles.
    
    Returns
    -------
    pcomp : str
        A string in NASTRAN PCOMP format
    '''
    pass