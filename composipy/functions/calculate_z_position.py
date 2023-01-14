'''
This function reproduces the geometry of an N-Layered Laminate.
See Jones figure 4.8
'''


def calculate_z_position(stacking, thickness):
    '''
    Parameters
    ----------
    stacking : list
        A list containing the angles of the stacking sequence.
    thickness : float or list
        A list containing the thickness of each ply or
        a float containing the thickness of all plies.
    
    Returns
    -------
    z_position : list
        A list containing the position of the ply.
    '''
    number_of_plies = len(stacking)
    
    if (isinstance(thickness, float) 
            or isinstance(thickness, int)):
        thickness = [thickness] * number_of_plies
    
    total_thickness = sum(thickness)
    current_z = -total_thickness/2
    z_position = [current_z]
    for t in thickness:
        current_z += t
        z_position.append(current_z)

    return z_position


