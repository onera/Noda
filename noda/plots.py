# Copyright 2025 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Generate plots."""

from copy import deepcopy

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import matplotlib as mpl

mpl.rcParams['axes.formatter.limits'] = (-2, 3)
mpl.rcParams['axes.grid'] = True
mpl.rcParams['grid.alpha'] = 0.5


def process_xaxis_unit(zunit):
    """
    Convert x-axis unit label to plotting-ready form and make multiplier.

    Parameters
    ----------
    zunit : str
        Unit in which distance is to be plotted on x axis.

    Returns
    -------
    zunit : str
        Unit in which distance is to be plotted on x axis.
    zmult : float
        Factor by which distance (x axis) will be multiplied before plotting.

    """
    zmult = None
    if zunit == 'nm':
        zmult = 1e9
    elif zunit == 'um':
        zmult = 1e6
        zunit = r'$\mu$m'
    return zunit, zmult


def process_yaxis_unit(varname):
    """
    Pre-process info to be used on y axis.

    Convert y-axis variable name to plotting-ready form and make y-axis unit
    label and multiplier.

    Parameters
    ----------
    zunit : str
        Unit in which distance is to be plotted.
    varname : str
        Name of variable to be plotted on y axis.

    Returns
    -------
    varname : str
        Name of variable to be plotted on y axis.
    varunit : float
        Unit of y-axis variable.
    varmult : float
        Factor by which y-axis variable will be multiplied before plotting.

    """
    varmult = 1
    varunit = ''
    if varname.startswith('J'):
        varunit = 'mol m$^{-2}$ s$^{-1}$'
        varname = rf'J^\mathrm{{{varname[1:]}}}'
    elif varname == 'mu':
        varunit = 'kJ/mol'
        varname = r'\mu'
        varmult = 1e-3
    elif varname == 'deformation':
        varunit = '-'
        varname = r'\varepsilon'
    elif varname == 'v':
        varunit = 'm/s'
    elif varname == 'c':
        varunit = 'mol/m3'
    return varname, varunit, varmult


def plot_profile_single(z, var, varname, title, zunit='um'):
    """
    Plot variable as a function of distance.

    Parameters
    ----------
    z : 1D array
        Positions to be used on x axis. Can be either node positions (size
        `nz`) or midpoint positions (size `nz` - 1) depending on the y-axis
        variable.
    var : 1D array or dict of 1D arrays
        Quantity to be plotted on y-axis.
    varname : str
        Name of y-axis variable.
    title : str
        Plot title.
    zunit : str, optional
        x-axis unit. The default is 'um'.

    Returns
    -------
    fig, ax : matplotlib figure and axis

    """
    zunit, zmult = process_xaxis_unit(zunit)
    varname, varunit, varmult = process_yaxis_unit(varname)
    ylabel = None
    fig, ax = plt.subplots()
    if isinstance(var, np.ndarray):
        ax.plot(z*zmult, var*varmult)
        ylabel = f'${varname}$ ({varunit})'
    elif isinstance(var, dict):
        for k in var:
            color = 'k' if k == 'Va' else None
            ax.plot(z*zmult, var[k]*varmult, c=color, label=k)
        ylabel = f'${varname}_i$ ({varunit})'
        ax.legend()
    ax.set_xlabel(f'$z$ ({zunit})')
    ax.set_ylabel(ylabel)
    ax.ticklabel_format(style='sci', scilimits=(-2, 4))
    ax.set_title(title)
    return fig, ax


def add_profile(z, var, varname, zunit='um', extra_legend=None):
    """
    Add set of profiles to existing plot.

    If extra_legend arg is provided, a legend entry is added with a black
    symbol, which applies to all added profiles.

    Parameters
    ----------
    z : 1D array
        Positions to be used on x axis. Can be either node positions (size
        `nz`) or midpoint positions (size `nz` - 1) depending on the y-axis
        variable.
    var : dict of 1D arrays
        Quantities to be added to y-axis.
    varname : str
        Name of y-axis variable.
    zunit : str, optional
        x-axis unit. The default is 'um'.
    extra_legend : str, optional
        Common legend to all added profiles. The default is None.

    """
    zunit, zmult = process_xaxis_unit(zunit)
    varname, _, varmult = process_yaxis_unit(varname)
    ax = plt.gca()
    ax.set_prop_cycle(None)
    for k in var:
        ax.plot(z*zmult, var[k]*varmult, '--k', mfc='none')
    if extra_legend:
        ax.plot([], [], '--k', mfc='none', label=extra_legend)
        ax.legend()


def plot_profile_quartet(result, title, suptitle,
                         zunit='um', ylim=None, exclude_dep=True):
    """
    Plot composition and flux profiles.

    Produce a 2 x 2 grid of subplots with the following variables:

    * Atom fraction
    * Relative difference between simulated and equilibrium vacancy site
      fraction
    * Flux in the lattice-fixed frame
    * Volume fraction of pores.

    Parameters
    ----------
    result : :class:`results.UnitResult`
        Simulation results at given time step.
    title : str
        Axis title.
    suptitle : str
        Plot title.
    zunit : str, optional
        x-axis unit. The default is 'um'.
    ylim : dict, optional
        Lower and upper limits of the y-axis of the four subplots. Four keys
        are recognized: 'x', 'y0', 'J', 'fp'. The values are to be given as a
        tuple or list of two floats. The default is None.
    exclude_dep : bool, optional
        Exclude the dependent constituent from the atom fraction subplots. The
        default is True.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure.
    axes : list of matplotlib.axes._subplots.AxesSubplot
        Axes.
    lines : list of lists of matplotlib.lines.Line2D
        Lines.

    """
    zunit, zmult = process_xaxis_unit(zunit)
    z = result.z*zmult
    x = result.x
    Jlat = result.Jlat
    ryVa = result.ryVa

    fig = plt.figure(figsize=(14, 8), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[2, 1])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[0, 1])
    ax4 = fig.add_subplot(gs[1, 1])

    lines1 = []
    lines2 = []
    lines3 = []
    lines4 = []

    fig.suptitle(suptitle)
    ax1.set_title(title, loc='left')

    x_plot_list = list(x)[:-1] if exclude_dep else list(x)
    for k in x_plot_list:
        line, = ax1.plot(z, x[k], label=k)
        lines1.append(line)
    ax1.legend(loc='upper right')
    ax1.set_xticklabels([])
    ax1.set_ylabel('$x_k$')

    line, = ax2.plot(z, ryVa, c='k')
    lines2.append(line)
    ax2.set_xlabel(f'$z$ ({zunit})')
    ylabel = r'$(y_\mathrm{Va} - y_\mathrm{Va}^\mathrm{eq})'
    ylabel += r' / y_\mathrm{Va}^\mathrm{eq}$'
    ax2.set_ylabel(ylabel)
    ax2.yaxis.set_label_coords(-0.1, 0.5)

    line, = ax3.plot(z, Jlat['Va'], c='k', label='Va')
    lines3.append(line)
    J_plot_list = [k for k in Jlat if k != 'Va']
    for k in J_plot_list:
        line, = ax3.plot(z, Jlat[k], label=k)
        lines3.append(line)
    ax3.set_ylabel('$J_k^{lat}$ (mol m$^{-2}$ s$^{-1}$)')
    ax3.legend(loc='upper right')
    ax3.set_xticklabels([])
    ax3.ticklabel_format(axis='y', style='sci', scilimits=(-2, 4))

    line, = ax4.plot(z, result.fp*100, 'k')
    lines4.append(line)
    ax4.set_xlabel(f'$z$ ({zunit})')
    ax4.set_ylabel('$f_p$ (%)')

    axes = [ax1, ax2, ax3, ax4]
    lines = [lines1, lines2, lines3, lines4]

    if ylim is not None:
        for ax, k in zip([ax1, ax2, ax3, ax4], ['x', 'y0', 'J', 'fp']):
            ylim[k] = ylim.get(k, None)
            if ylim[k]:
                ax.set_ylim(ylim[k])

    return fig, axes, lines


class StaticProfile:
    """
    Static plot of composition and/or flux profile.

    Attributes
    ----------
    res : :class:`results.UnitResult`
        Simulation results at one time step.
    step : int
        Time step.
    th : float
        Time (h).
    title : str
        Default plot title.

    """

    def __init__(self, res, step, th, title):
        """
        Class constructor.

        Parameters
        ----------
        res : :class:`results.UnitResult`
            Simulation results at one time step.
        th : float
            Time (h).
        step : int
            Time step.
        title : str
            Default plot title.

        """
        self.res = res
        self.step = step
        self.th = th
        self.title = title

    def single(self, varname='x', title=None, plot_dep=False, **kwargs):
        """
        Plot profile at given time step.

        Call :func:`plots.plot_profile_single`.
        The plotted variable can be multivalued (dict) or univalued (1D array).

        Parameters
        ----------
        varname : str, optional
            Name of variable to be plotted. The default is `x`.
        plot_dep : bool
            Whether to include dependent constituent in `x` and `y` plots.
        kwargs
            Optional arguments passed to :func:`plot_profile_single`.

        Returns
        -------
        fig, ax : matplotlib figure and axis

        Raises
        ------
        Exception
            If variable is not recognized.

        """
        var = getattr(self.res, varname)
        if plot_dep is False:
            if varname in ['x', 'y']:
                inds = list(var.keys())[:-1]
                var = {k: var[k] for k in inds}

        if title is None:
            title = self.title + f'\nstep {self.step:3}, {self.th:3.1f} h'
        fig, ax = plot_profile_single(self.res.z, var, varname, title, **kwargs)
        return fig, ax

    def quartet(self, title=None, **kwargs):
        """
        Plot composition and flux profiles at given time step.

        Call :func:`plots.plot_profile_quartet`. Plot x, yVa, Jlat, fp.

        Parameters
        ----------
        kwargs
            Optional arguments passed to :func:`plot_profile_quartet`.

        Returns
        -------
        fig, axes, lines : matplotlib figure, axes and lines

        """
        if title is None:
            title = f'step {self.step:3}, {self.th:3.1f} h'
        suptitle = self.title
        fig, axes, lines = plot_profile_quartet(self.res, title, suptitle,
                                                **kwargs)
        return fig, axes, lines


def calculate_view_limits(results, varname):
    """
    Compute y-axis limits adapted to results at all time steps.

    Parameters
    ----------
    results : dict of :class:`results.UnitStep`
        Contains simulation results with time steps as keys.
    varname : str
        Name of variable of interest.

    Returns
    -------
    ylim : dict of lists
        y-axis view limits for variable of interest.

    """
    vals = [getattr(r, varname) for r in results.values()]
    varname, _, varmult = process_yaxis_unit(varname)
    ymin = 0
    ymax = 0

    if isinstance(vals[0], np.ndarray):
        ymin = min(v.min() for v in vals)
        ymax = max(v.max() for v in vals)
    elif isinstance(vals[0], dict):
        ymin = min(np.array(list(v.values())).min() for v in vals)
        ymax = max(np.array(list(v.values())).max() for v in vals)

    ymin *= varmult
    ymax *= varmult

    ymid = (ymin + ymax)/2
    new_span = (ymax - ymin)*1.1
    # set arbitrary span to avoid matplotlib warning if bottom == top == 0
    if new_span == 0:
        new_span = 1
    ymin = ymid - new_span/2
    ymax = ymid + new_span/2

    return ymin, ymax


class InteractivePlot():
    """Interactive plot of simulation results."""

    def __init__(self, varname, comps, results, saved_times, saved_steps,
                 title):
        """
        Class constructor.

        Parameters
        ----------
        varname : str
            Name of variable to be plotted.
        comps : list of str
            System constituents.
        results : dict of :class:`results.UnitResult`
            Simulation results.
        saved_steps : list
            Steps for which simulation results are stored in steps.
        saved_times : 1D array
            Times in h (rounded) that correspond to the saved steps.
        title : str
            Plot title.

        """
        self.comps = comps
        self.saved_times = saved_times
        self.saved_steps = saved_steps
        self.num_out = len(saved_steps)
        self.i = 0
        self.results = deepcopy(results)
        self.varname = varname

        # Remove dependent component for x profile
        if varname == 'x':
            for r in self.results.values():
                r.x = {k: r.x[k] for k in comps[1:-1]}

        xlim = [results[0].z[0]*1e6, results[0].z[-1]*1e6]
        ylim = calculate_view_limits(self.results, varname)

        self.fig, self.ax, self.lines = self.setup_plot(title, xlim, ylim)

        slider_ax = plt.axes([0.4, 0.895, 0.4, 0.02],
                             facecolor='lightgoldenrodyellow')
        self.slider = Slider(slider_ax, '', 0, self.num_out - 1,
                             valinit=0, valstep=1)
        self.slider.valtext.set_visible(False)

        self.fig.canvas.draw()
        self.bg = self.fig.canvas.copy_from_bbox(self.fig.bbox)
        self.update()

        self.fig.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.slider.on_changed(self.update_slider)

    def setup_plot(self, title, xlim, ylim):
        """
        Make figure with empty lines.

        Parameters
        ----------
        title : str
            Figure title.
        xlim : list of floats
            x-axis view limits.
        ylim : dict of lists
            y-axis view limits.

        Returns
        -------
        fig : matplotlib.figure.Figure
            Figure.
        ax : matplotlib.axes._subplots.AxesSubplot
            Axis.
        lines : list of matplotlib.lines.Line2D
            Lines.

        """
        fig, ax = plt.subplots(figsize=(7, 8/1.5), constrained_layout=True)
        fig.suptitle(title)
        ax.set_title(' ', loc='left')

        lines = []

        var = getattr(self.results[0], self.varname)
        varname, varunit, varmult = process_yaxis_unit(self.varname)
        self.varmult = varmult
        ylabel = None

        if isinstance(var, np.ndarray):
            line, = ax.plot([])
            ylabel = f'${varname}$ ({varunit})'
            lines = [line]
        elif isinstance(var, dict):
            for k in var:
                color = 'k' if k == 'Va' else None
                line, = ax.plot([], c=color, label=k)
                lines.append(line)
            ylabel = f'${varname}_i$ ({varunit})'

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.legend(loc='upper right')
        ax.set_xlabel(r'$z$ ($\mu$m)')
        ax.set_ylabel(ylabel)
        ax.ticklabel_format(style='sci', scilimits=(-2, 4))

        return fig, ax, lines

    def on_scroll(self, event):
        """Define callback used on scrolling."""
        if event.button == 'up':
            self.i = (self.i + 1) % self.num_out
        else:
            self.i = (self.i - 1) % self.num_out

        self.slider.set_val(self.i)
        self.slider.ax.figure.canvas.draw_idle()
        self.update()

    def update_slider(self, val):
        """Set slider to given position (time step) and update the figure."""
        self.i = int(val)
        self.update()

    def update(self):
        """Update and draw the figure."""
        th = self.saved_times[self.i]
        n = self.saved_steps[self.i]
        res = self.results[n]
        var = getattr(res, self.varname)
        z = res.z*1e6

        title = self.ax.set_title(f'step {n:3}, {th:3.1f} h', loc='left')

        if isinstance(var, np.ndarray):
            self.lines[0].set_data(z, var*self.varmult)
        elif isinstance(var, dict):
            for k, line in zip(var, self.lines):
                line.set_data(z, var[k]*self.varmult)

        self.fig.canvas.restore_region(self.bg)
        for line in self.lines:
            self.ax.draw_artist(line)
        self.ax.draw_artist(title)
        self.fig.canvas.blit(self.fig.bbox)


class InteractivePlotQuartet():
    """
    Interactive plot of composition and flux profiles.

    The layout and quantities plotted are the same as in
    :func:`plot_profiles_and_fluxes`.

    """

    def __init__(self, comps, results, saved_times, saved_steps, title):
        """
        Class constructor.

        Parameters
        ----------
        comps : list of str
            System constituents.
        results : dict of :class:`results.UnitResult`
            Simulation results.
        saved_steps : list
            Steps for which simulation results are stored in steps.
        saved_times : 1D array
            Times in h (rounded) that correspond to the saved steps.
        title : str
            Plot title.

        """
        self.comps = comps
        self.results = deepcopy(results)
        self.saved_times = saved_times
        self.saved_steps = saved_steps
        self.num_out = len(saved_steps)
        self.i = 0

        # Remove dependent component for x profile
        for r in self.results.values():
            r.x = {k: r.x[k] for k in comps[1:-1]}

        xlim = [results[0].z[0]*1e6, results[0].z[-1]*1e6]
        ylim = {varname: calculate_view_limits(self.results, varname)
                for varname in ['x', 'ryVa', 'fp', 'Jlat']}

        self.fig, self.axes, self.lines = self.setup_plot(title,
                                                          xlim, ylim)
        self.ax1, self.ax2, self.ax3, self.ax4 = self.axes
        self.lines1, self.lines2, self.lines3, self.lines4 = self.lines

        slider_ax = plt.axes([0.25, 0.95, 0.5, 0.02],
                             facecolor='lightgoldenrodyellow')
        self.slider = Slider(slider_ax, '', 0, self.num_out - 1,
                             valinit=0, valstep=1)
        self.slider.valtext.set_visible(False)

        self.fig.canvas.draw()
        self.bg = self.fig.canvas.copy_from_bbox(self.fig.bbox)
        self.update()

        self.fig.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.slider.on_changed(self.update_slider)

    def setup_plot(self, title, xlim, ylim):
        """
        Make figure with empty lines.

        Parameters
        ----------
        title : str
            Figure title.
        xlim : list of floats
            x-axis view limits.
        ylim : dict of lists
            y-axis view limits (see :func:`plot_profiles_and_fluxes` for keys).

        Returns
        -------
        fig : matplotlib.figure.Figure
            Figure.
        axes : list of matplotlib.axes._subplots.AxesSubplot
            Axes.
        lines : list of lists of matplotlib.lines.Line2D
            Lines.

        """
        fig = plt.figure(figsize=(14, 8), constrained_layout=True)
        gs = fig.add_gridspec(2, 2, height_ratios=[2, 1])
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[0, 1])
        ax4 = fig.add_subplot(gs[1, 1])

        fig.suptitle(title)
        ax1.set_title(' ', loc='left')

        lines1 = []
        lines2 = []
        lines3 = []
        lines4 = []

        for k in self.comps[1:-1]:
            line, = ax1.plot([], [], label=k)
            lines1.append(line)
        ax1.set_xlim(xlim)
        ax1.set_ylim(ylim['x'])
        ax1.legend(loc='upper right')
        ax1.set_xticklabels([])
        ax1.set_ylabel('$x_k$')

        line, = ax2.plot([], [], c='k')
        lines2.append(line)
        ax2.set_xlim(xlim)
        ax2.set_ylim(ylim['ryVa'])
        ax2.set_xlabel(r'$z$ ($\mu$m)')
        ylabel = r'$(y_\mathrm{Va} - y_\mathrm{Va}^\mathrm{eq})'
        ylabel += r' / y_\mathrm{Va}^\mathrm{eq}$'
        ax2.set_ylabel(ylabel)
        ax2.yaxis.set_label_coords(-0.1, 0.5)

        line, = ax3.plot([], [], 'k', label='Va')
        lines3.append(line)
        for k in self.comps[1:]:
            line, = ax3.plot([], [], label=k)
            lines3.append(line)
        ax3.set_xlim(xlim)
        ax3.set_ylim(ylim['Jlat'])
        ax3.set_ylabel('$J_k^{lat}$ (mol m$^{-2}$ s$^{-1}$)')
        ax3.legend(loc='upper right')
        ax3.set_xticklabels([])
        ax3.ticklabel_format(axis='y', style='sci', scilimits=(-2, 4))

        line, = ax4.plot([], [], 'k')
        lines4.append(line)
        ax4.set_xlim(xlim)
        ax4.set_ylim([yl*100 for yl in ylim['fp']])
        ax4.set_xlabel(r'$z$ ($\mu$m)')
        ax4.set_ylabel('$f_p$ (%)')

        axes = [ax1, ax2, ax3, ax4]
        lines = [lines1, lines2, lines3, lines4]

        return fig, axes, lines

    def on_scroll(self, event):
        """Define callback used on scrolling."""
        if event.button == 'up':
            self.i = (self.i + 1) % self.num_out
        else:
            self.i = (self.i - 1) % self.num_out

        self.slider.set_val(self.i)
        self.slider.ax.figure.canvas.draw_idle()
        self.update()

    def update_slider(self, val):
        """Set slider to given position (time step) and update the figure."""
        self.i = int(val)
        self.update()

    def update(self):
        """Update and draw the figure."""
        th = self.saved_times[self.i]
        n = self.saved_steps[self.i]
        res = self.results[n]
        z = res.z*1e6
        x = res.x
        Jlat = res.Jlat

        title = self.ax1.set_title(f'step {n:3}, {th:3.1f} h', loc='left')

        for k, line in zip(self.comps[1:], self.lines1):
            line.set_data(z, x[k])

        self.lines2[0].set_data(z, res.ryVa)

        for k, line in zip(Jlat, self.lines3):
            line.set_data(z, Jlat[k])

        self.lines4[0].set_data(z, res.fp*100)

        self.fig.canvas.restore_region(self.bg)
        for ax, lines in zip(self.axes, self.lines):
            for line in lines:
                ax.draw_artist(line)
        self.ax1.draw_artist(title)
        self.fig.canvas.blit(self.fig.bbox)


def plot_mass_balance(results):
    """
    Plot mass balance.

    Parameters
    ----------
    results : :class:`results.SimulationResults`
        Simulation results.

    Returns
    -------
    fig, ax : matplotlib figure and axis

    """
    fig, ax = plt.subplots()
    for i, k in zip(range(1, len(results.comps)), results.comps[1:]):
        rdN = [results.steps[n].rdIN[i] for n in results.saved_steps]
        ax.plot(results.saved_times, rdN, label=k)
    rdNv = [results.steps[n].rdINv for n in results.saved_steps]
    ax.plot(results.saved_times, rdNv, label='Va + pore + void')
    ax.legend()
    ax.set_xlabel('$t$ (h)')
    ax.set_ylabel('$(N_t - N_0)/N_0$')
    return fig, ax
