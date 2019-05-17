import io
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import logging


def rcparam():
    plt.rc('axes', linewidth=4.0, labelsize=18)
    # axes and axes labels properties
    plt.rc('xtick.major', size=10)  # length for x
    plt.rc('ytick.major', size=10)  # length for y
    plt.rc('lines', mew=5)
    plt.rc('lines', lw=4)  # line thickness
    plt.rc('ytick', labelsize=12)  # ytick label size
    plt.rc('xtick', labelsize=12)  # xtick label size
    plt.rc('xtick.major', pad=5)  # xtick padding distance
    plt.rc('ytick.major', pad=5)  # ytick padding distance


def plot_2d_heatmap(data, title=None, axis_labels=None, axis_tick_format='%.0f', force_matching_ticks=True, save_fig=False, fig_name='HeatMap.png', **kwargs):
    rcparam()
    font = {'fontsize': 20,
            'fontweight' : 'bold',
            'verticalalignment': 'baseline'}

    # Figure Creation
    fig, ax = plt.subplots()
    if 'invert_color' in kwargs:
        if kwargs['invert_color']:
            plt.set_cmap('viridis_r')
    if 'absolute_max' in kwargs:
        cax = ax.pcolormesh(data, vmin=0., vmax=kwargs['absolute_max'])
    else:
        cax = ax.pcolormesh(data)
    cbar = fig.colorbar(cax)

    # Adding labels
    if title:
        plt.title(title, fontdict=font)
    if axis_labels:
        plt.xlabel(axis_labels[0])
        plt.ylabel(axis_labels[1])
    if force_matching_ticks:
        xticks = ax.get_xticks()[0:-1]
        ax.set_yticks(xticks)

    ax.xaxis.set_major_formatter(FormatStrFormatter(axis_tick_format))
    ax.yaxis.set_major_formatter(FormatStrFormatter(axis_tick_format))

    if save_fig:
        plt.savefig(fig_name, transparent=True, bbox_inches='tight', pad_inches=0.1)
    plt.close()


def plot_2d_heatmap_buffer(data, title=None, axis_labels=None, axis_tick_format='%.0f', force_matching_ticks=True, save_fig=False, fig_name='HeatMap.png', **kwargs):
    rcparam()
    font = {'fontsize': 20,
            'fontweight' : 'bold',
            'verticalalignment': 'baseline'}

    # Figure Creation
    fig, ax = plt.subplots()
    if kwargs.get('invert_color', False):
        plt.set_cmap('viridis_r')
    else:
        plt.set_cmap('viridis')

    if 'absolute_max' in kwargs:
        cax = ax.pcolormesh(data, vmin=0., vmax=kwargs['absolute_max'])
    else:
        cax = ax.pcolormesh(data)
    cbar = fig.colorbar(cax)

    # Adding labels
    if title:
        plt.title(title, fontdict=font)
    if axis_labels:
        plt.xlabel(axis_labels[0])
        plt.ylabel(axis_labels[1])
    if force_matching_ticks:
        xticks = ax.get_xticks()[0:-1]
        ax.set_yticks(xticks)

    ax.xaxis.set_major_formatter(FormatStrFormatter(axis_tick_format))
    ax.yaxis.set_major_formatter(FormatStrFormatter(axis_tick_format))

    buf = io.BytesIO()
    if save_fig:
        plt.savefig(buf, transparent=True, dpi=75, bbox_inches='tight', pad_inches=0.1, format='png')
    plt.close()
    return buf


def plot_2d_heatmap_dynamic_axis(data, param_dict, plot_group, save_fig=False, fig_name='HeatMap.png'):
    rcparam()
    # Tick/axis label creation and formatting
    label_x = plot_group[1]
    label_y = plot_group[0]
    axis_tick_count = []
    axis_tick_labels = []
    for parameter in plot_group:
        axis_tick_count.append(len(param_dict[parameter]))
        axis_tick_labels.append(['{:.2e}'.format(x) for x in param_dict[parameter]])

    # Figure Creation
    fig, ax = plt.subplots()
    cax = ax.pcolor(data)
    cbar = fig.colorbar(cax)

    # Tick/axis creation and label assignment
    ax.set_xticks(np.arange(0.5, axis_tick_count[1], 1))
    ax.set_yticks(np.arange(0.5, axis_tick_count[0], 1))
    ax.set_xticklabels(axis_tick_labels[1], rotation=45)
    ax.set_yticklabels(axis_tick_labels[0])
    for tick in ax.xaxis.get_majorticklabels():
        tick.set_horizontalalignment("right")
    ax.set_xlabel(label_x)
    ax.set_ylabel(label_y)

    if save_fig:
        plt.savefig(fig_name, transparent=True, bbox_inches='tight', pad_inches=0.1, dpi=300)


if __name__ == "__main__":
    print('Script module used to create 2D heatmaps from simulation data')
