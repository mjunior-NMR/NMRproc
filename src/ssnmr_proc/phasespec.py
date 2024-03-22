# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 10:31:48 2024

@author: Marcos
"""
import numpy as np

def manual(data):
    """
    Manual Phase correction using matplotlib

    A matplotlib widget is used to manually correct the phase of a Fourier
    transformed dataset. If the dataset has more than 1 dimensions, the first
    trace will be picked up for phase correction.  Clicking the 'Set Phase'
    button will print the current linear phase parameters to the console.
    A ipywidget is provided for use with Jupyter Notebook to avoid changing
    backends. This can be accessed with notebook=True option in this function

    .. note:: Needs matplotlib with an interactive backend.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.    

    Returns
    -------
    p0, p1 : float
        Linear phase correction parameters. Zero and first order phase
        corrections in degrees calculated from pc0, pc1 and pivot displayed
        in the interactive window.

    Examples
    --------
    >>> import nmrglue as ng
    >>> p0, p1 = ng.process.proc_autophase.manual_ps(data)
    >>> # do manual phase correction and close window
    >>> phased_data = ng.proc_base.ps(data, p0=p0, p1=p1)

    """

    if len(data.shape) == 2:
        data = data[0, ...]
    elif len(data.shape) == 3:
        data = data[0, 0, ...]
    elif len(data.shape) == 4:
        data = data[0, 0, 0, ...]

    from matplotlib.widgets import Slider, Button
    import matplotlib.pyplot as plt
    plt.ioff()
    plt.subplots_adjust(left=0.25, bottom=0.35)
    
    interactive, = plt.plot(data.real, lw=1, color='black')
    
    axcolor = 'white'
    axpc0 = plt.axes([0.25, 0.10, 0.65, 0.03], facecolor=axcolor)
    axpc1 = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
    axpiv = plt.axes([0.25, 0.20, 0.65, 0.03], facecolor=axcolor)
    axpst = plt.axes([0.25, 0.25, 0.15, 0.04], facecolor=axcolor)
    axpst1 = plt.axes([0.75, 0.25, 0.15, 0.04], facecolor=axcolor)
    
    spc0 = Slider(axpc0, 'p0', -360, 360, valinit=0)
    spc1 = Slider(axpc1, 'p1', -3600, 3600, valinit=0)
    spiv = Slider(axpiv, 'pivot', 0, data.size, valinit=0)
    axps = Button(axpst, 'Set Phase', color=axcolor)
    done = Button(axpst1, 'Done', color=axcolor)
    
    def update(val):
        pc0 = spc0.val * np.pi / 180
        pc1 = spc1.val * np.pi / 180
        pivot = spiv.val
        interactive.set_ydata((data * np.exp(
            1.0j * (pc0 + (pc1 * np.arange(-pivot, -pivot + data.size) /
                    data.size))).astype(data.dtype)).real)
        plt.draw()
    
    def setphase(val):
        p0 = spc0.val-spc1.val*spiv.val/data.size
        p1 = spc1.val
        print(p0, p1)
    
    def return_ps(event):
        plt.show(block=False)
        return #from this and outer function...

    
    spc0.on_changed(update)
    spc1.on_changed(update)
    spiv.on_changed(update)
    axps.on_clicked(setphase)
    done.on_clicked(return_ps)
    
    plt.show(block=True)
    
    p0 = spc0.val-spc1.val*spiv.val/data.size
    p1 = spc1.val
    return p0, p1