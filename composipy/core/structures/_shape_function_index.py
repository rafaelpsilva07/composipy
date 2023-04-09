'''
This module contains functions to treat shape functions indexes and boundary conditions.
'''


from itertools import product


def _generate_index(constraints, m, n, idx_option='ijkl'):
    '''
    Paramenters
    -----------
    constraints : str or dict. Default: \"PINNED\"
        Plate boundary conditions.
        Options: "PINNED", "CLAMPED" or dict
        constraints = {    
            x0 = ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']
            xa = ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']
            y0 = ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']
            yb = ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']}
     m : int, default 10
        Size of shape function along x axis
    n : int, default 10
        Size of shape function along y axis
    plane : str , default 'xy'
        Orientation plane of the plate
    idx_option : str , default 'ijkl'
        'ijkl' -> [0000, 0001, ..., ijkl]
        'ij' -> [00, 01, ..., ij]

    Returns
    -------
    (s1idx, s2idx, s3idx) : tuple
        Each element of the tuple contains a list of index.
        Ex.: [0000, 0001, ..., ijkl]
    '''

    if constraints == 'PINNED':
        x0 = ['TX', 'TY', 'TZ']
        xa = ['TX', 'TY', 'TZ']
        y0 = ['TX', 'TY', 'TZ']
        yb = ['TX', 'TY', 'TZ']
    elif constraints == 'CLAMPED':
        x0 = ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']
        xa = ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']
        y0 = ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']
        yb = ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']
    else:
        x0 = constraints['x0']
        xa = constraints['xa']
        y0 = constraints['y0']
        yb = constraints['yb']

    sm = [i for i in range(m+4)]
    sn = [i for i in range(n+4)]

    um, un = sm.copy(), sn.copy()
    vm, vn = sm.copy(), sn.copy()
    wm, wn = sm.copy(), sn.copy()

    # x0
    if 'TX' in x0:
        um.remove(0)
    if 'TY' in x0:
        vm.remove(0)
    if 'TZ' in x0:
        wm.remove(0)
    if 'RX' in x0:
        um.remove(1)
    if 'RY' in x0:
        vm.remove(1)
    if 'RZ' in x0:
        wm.remove(1)
    #xa
    if 'TX' in xa:
        um.remove(2)
    if 'TY' in xa:
        vm.remove(2)
    if 'TZ' in xa:
        wm.remove(2)
    if 'RX' in xa:
        um.remove(3)
    if 'RY' in xa:
        vm.remove(3)
    if 'RZ' in xa:
        wm.remove(3)
    
    #y0
    if 'TX' in y0:
        un.remove(0)
    if 'TY' in y0:
        vn.remove(0)
    if 'TZ' in y0:
        wn.remove(0)
    if 'RX' in y0:
        un.remove(1)
    if 'RY' in y0:
        vn.remove(1)
    if 'RZ' in y0:
        wn.remove(1)
    
    #yb
    if 'TX' in yb:
        un.remove(2)
    if 'TY' in yb:
        vn.remove(2)
    if 'TZ' in yb:
        wn.remove(2)
    if 'RX' in yb:
        un.remove(3)
    if 'RY' in yb:
        vn.remove(3)
    if 'RZ' in yb:
        wn.remove(3)
    
    um, un = um[0: m], un[0: n]
    vm, vn = vm[0: m], un[0: n]
    wm, wn = wm[0: m], un[0: n]

    if idx_option == 'ijkl':
        uidx = list(product(um, un, um, un))
        vidx = list(product(vm, vn, vm, vn))
        widx = list(product(wm, wn, wm, wn))
    elif idx_option == 'ij':
        uidx = list(product(um, un))
        vidx = list(product(vm, vn))
        widx = list(product(wm, wn))
    else:
        raise ValueError('idx_option must be ijkl or ij')

    return (uidx, vidx, widx)


if __name__ == '__main__':
        uidx, vidx, widx = _generate_index("CLAMPED", 5, 5)
        print(uidx)
