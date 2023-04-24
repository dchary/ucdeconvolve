####################################################################################################
# # Copyright (C) 2023-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

from matplotlib.colors import ListedColormap
from enum import Enum

import numpy as np

def create_colormap(
    R : int,
    G : int,
    B : int,
    whitespace : int = 20
) -> ListedColormap:
    """\
    Create Colormap
    
    Creates a transparent colormap for plotting spatial
    data using a base RGB color. Inspired by similar
    function used by cell2location [Kleshchevnikov et. al. 2022]
    
    Params
    ------
    R
        integer for red channel 0 - 255
    G
        integer for green channel 0 - 255
    B
        integer for blue channel 0 - 255
    whitespace
        percent of colormap (expressed as int 0 -255)
        that will be transparent. default is 20%
    
    Returns
    -------
    cmap : ListedColormap
        a listed colormap of the provided colors.

    """
    
    spacing = int(whitespace * 2.55)

    N = 255
    M = 3
    
    alphas = np.concatenate([[0] * spacing * M, np.linspace(0, 1.0, (N - spacing) * M)])
    vals = np.ones((N * M, 4))
    
    for i, color in enumerate([R, G, B]):
        vals[:, i] = color / 255
    vals[:, 3] = alphas

    return ListedColormap(vals)


class CM(Enum):
    """\
    
    Class containing enums representing default colormaps.
    
    """
    
    Yellow = create_colormap(240, 228, 66)
    Orange = create_colormap(213, 94, 0)
    Blue = create_colormap(86, 180, 233)
    Green = create_colormap(0, 158, 115)
    GreenNeon = create_colormap(255, 87, 51)
    Grey = create_colormap(200, 200, 200)
    White = create_colormap(50, 50, 50)
    Purple = create_colormap(90, 20, 165)
    LightPurple = create_colormap(204,153,255)
    LightGreen = create_colormap(189,236,182)
    Yellow2 = create_colormap(254,225,40)
    Red = create_colormap(130,20,49)
    DarkBlue = create_colormap(0,0,139)