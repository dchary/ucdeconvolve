####################################################################################################
# # Copyright (C) 2023-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

# Plotting Utilities
import matplotlib as mpl
import matplotlib.pyplot as plt


def append_custom_cbar(fig, old_cbar, rect, cmap, orientation = 'vertical', extend = 'both', label = ''):
    
    ylims = old_cbar.ax.get_ylim()
    xlims = old_cbar.ax.get_xlim()
    
    new_ax = fig.add_axes(rect)
    
    new_cbar = mpl.colorbar.ColorbarBase(new_ax, 
            orientation=orientation, 
            cmap=cmap,
            norm=mpl.colors.Normalize(
                ylims[0] if orientation == 'vertical' else xlims[0],
                ylims[1] if orientation == 'vertical' else xlims[1]),
            extend=extend,
            label=label)
    
    return new_ax