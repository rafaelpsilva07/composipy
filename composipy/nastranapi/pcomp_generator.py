
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
    block = list()
    current_block = ''
    is_in_block = False
    is_multiplier = False
    
    for i, char in enumerate(stacking):
        if char == '(':
            is_in_block = True
        elif char == ')':
            unity = _build_stack_unit(current_block)
            block.extend(unity)
            current_block = ''
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
              and is_in_block
              and not is_multiplier):
            unity = _build_stack_unit(current_block)
            block.extend(unity)
            current_block = ''
        elif (char == '/' 
             and is_multiplier
             and not is_in_block):
            times = current_block.replace('_', '').replace('{', '').replace('}', '')
            times = int(times)
            unity = times * unity
            unities.extend(unity)
            unity = list()
            current_block = ''
            is_multiplier = False
        elif (char == '/'
              and is_multiplier
              and is_in_block):
            times = current_block.replace('_', '').replace('{', '').replace('}', '')
            times = int(times)
            unity = times * block
            unities.extend(unity)
            unity = list()
            block = list()
            current_block = ''
            is_multiplier = False
            is_in_block = False
        else:
            current_block += char
           
    if is_symmetric:
        symmetric_unities = unities[::-1]
        unities.extend(symmetric_unities)
    return unities


def _convert_to_list(value, size):
    '''If it is list return value, if value is a number return a list'''
    try:
        size_value = len(value)
        if size_value == size:
            return value
        else:
            raise ValueError(f'{value} does not fit with sequence size')
    except TypeError:
        return [value for i in range(size)]

def _convert_sout(sout, size):
    if sout == 'FIBER':
        sout = ['NO' for i in range(size)]
        sout[0] = 'YES'
        sout[-1] = 'YES'
        return sout
    elif sout == 'NO':
        sout = ['NO' for i in range(size)]
        return sout
    elif sout == 'YES':
        sout = ['YES' for i in range(size)]
        return sout
    else:
        try:
            size_value = len(sout)
            if size_value == size:
                return sout
            else:
                 raise ValueError
        except:
            raise ValueError(f'{sout} not accepted')
            

def build_pcomp(sequence, midi, ti, pid=1, z0='',sout='FIBER'):
    '''
    Parameters
    ----------
    sequence : list
        Angles of the stacking sequence. A list (or a iterable) containing angles
    midi : list or int
        A list of materials of plies or a material id (int) to be applied to all plies.
    ti : list or float
        A list of thickness of plies or a thickness (float) to be applied to all plies.
    pid : int, default 1
        NASTRAN property id
    z0 : float, default ''
        Laminate offset
    sout : list or str, default 'FIBER'
        A list of output request of plies or a output request to be applied to all plies.
        Options = YES, NO, FIBER
        If FIBER is used, then only the first and the last plies will be set as YES
    
    Returns
    -------
    text : str
        Returns a PCOMP card. 
    '''
    size = len(sequence)
    pid = int(pid)
    ti = _convert_to_list(ti, size)
    midi = _convert_to_list(midi, size)
    sout = _convert_sout(sout, size)
    
    header = f'PCOMP,{pid},{z0},,,,,+\n'
    body = ''
    
    for i in range(0, size, 2):
        i0 = i
        i1 = i + 1
        body += f'+,{midi[i0]},{ti[i0]},{sequence[i0]},{sout[i0]},'
        try:
            body += f'{midi[i1]},{ti[i1]},{sequence[i1]},{sout[i1]}+\n'
        except IndexError:
            body = body[0:-1] # removes last comma
    
    return header + body