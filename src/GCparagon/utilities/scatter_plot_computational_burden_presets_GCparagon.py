#!/usr/bin/env python3
import pathlib
import plotly.express as px
import pandas as pd
from scipy.stats import linregress as scpst_linregress
from math import floor, log10


REPO_ROOT_PATH = pathlib.Path(__file__).parent.parent.parent.parent
REPO_ROOT_DIR = str(REPO_ROOT_PATH)


color_sequence = ['rgba(51, 98, 255, 0.5)', 'rgba(0, 178, 18, 0.5)', 'rgba(255, 0, 0, 0.5)']
# blue, green, red


# Define function for string formatting of scientific notation; code taken from
# https://stackoverflow.com/questions/18311909/how-do-i-annotate-with-power-of-ten-formatting
def sci_notation(num, prefix=None, postfix=None, decimal_digits=1, precision=None, exponent=None):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    """
    if exponent is None:
        exponent = int(floor(log10(abs(num))))
    coeff = round(num / float(10**exponent), decimal_digits)
    if precision is None:
        precision = decimal_digits
    if prefix is None:
        prefix = ''
    if postfix is None:
        postfix = ''
    return r"${0}{1:.{3}f}\cdot10^{{{2:d}}}{4}$".format(prefix, coeff, exponent, precision, postfix)


if __name__ == '__main__':
    # sample	preset	iteration	memory_consumption	duration    processed_fragments
    out_path = REPO_ROOT_PATH / 'preset_computation/benchmark_results/'
    performance_metrics = pd.read_table(out_path / 'benchmarking_summary.tsv', sep='\t')  # NEEDS TO BE CREATED FIRST!
    # create new preset column
    performance_metrics = performance_metrics.assign(preset_str="preset " + performance_metrics.preset.astype(str))
    # create linear correlation plot
    lin_corr_duration_result = scpst_linregress(performance_metrics['processed_fragments'],
                                                performance_metrics['duration'])
    p_value_duration = float(lin_corr_duration_result.pvalue)
    r_value_duration_str = str(round(lin_corr_duration_result.rvalue, ndigits=3))[:-1]
    lin_corr_memory_result = scpst_linregress(performance_metrics['processed_fragments'],
                                              performance_metrics['memory_consumption'])
    p_value_memory = float(lin_corr_memory_result.pvalue)
    r_value_memory_str = str(round(lin_corr_memory_result.rvalue, ndigits=3))[:-1]
    # computation time figures
    # with trendline (only created to extract trendline)
    trendline_fig = px.scatter(x=performance_metrics['processed_fragments'], y=performance_metrics['duration'],
                               trendline='ols',
                               labels={'x': 'fragments processed / count', 'y': 'computation time / s'})
    trendline_data = trendline_fig.data[1]
    # colored dots with trendline figure
    duration_scatter_fig = px.scatter(x=performance_metrics['processed_fragments'], y=performance_metrics['duration'],
                                      color=performance_metrics['preset_str'].astype(str),
                                      symbol=performance_metrics['sample'], color_discrete_sequence=color_sequence,
                                      labels={'x': 'fragments processed / count', 'y': 'computation time / s'})
    duration_scatter_fig.update_traces(marker=dict(size=16, line=dict(width=2)))  # , opacity=0.5
    duration_scatter_fig.add_trace(trendline_data)
    # switch order: trendline in background
    duration_scatter_fig.data = [duration_scatter_fig.data[-1]] + list(duration_scatter_fig.data[:-1])
    # format trendline
    duration_scatter_fig.data[0].line.color = 'black'
    duration_scatter_fig.data[0].line.width = 4
    duration_scatter_fig.update_layout(width=1200, height=750, template='simple_white',
                                       font={'family': 'Ubuntu', 'size': 22,
                                             'color': 'rgb(20, 20, 20)'},
                                       annotations=[dict(x=int((max(performance_metrics['processed_fragments']) -
                                                                min(performance_metrics['processed_fragments']))
                                                               * 0.25), y=max(performance_metrics['duration']) * 0.8,
                                                         text=sci_notation(
                                                             p_value_duration, decimal_digits=2, postfix=r')',
                                                             prefix=f"Pearson's R={r_value_duration_str} (p="),
                                                         font=dict(family='Ubuntu', size=22, color='rgb(20,20,20)'),
                                                         showarrow=False, bgcolor='white')],
                                       title={'text': 'Dependency of Computation Time on Processed Fragments',
                                              'font': {'family': 'Ubuntu', 'size': 36,
                                                       'color': 'rgb(20, 20, 20)'},
                                              'xanchor': 'center', 'yanchor': 'middle', 'x': 0.5},
                                       legend=dict(orientation='h', xanchor='center', yanchor='top',
                                                   x=0.5, y=-0.18, title='', font=dict(family='Ubuntu', size=18)))
    duration_scatter_fig.update_xaxes(showgrid=True, gridcolor='lightgrey', zerolinecolor='rgb(20,20,20)')
    duration_scatter_fig.update_yaxes(showgrid=True, gridcolor='lightgrey', zerolinecolor='rgb(20,20,20)')
    duration_scatter_fig.write_image(out_path / f'GCparagon_computation_time_presetsAndSamples_linRegress_SSD-ref.png')

    # memory consumption figure
    memory_scatter_fig = px.scatter(x=performance_metrics['processed_fragments'],
                                    y=performance_metrics['memory_consumption'],
                                    color=performance_metrics['preset_str'], symbol=performance_metrics['sample'],
                                    title='Maximum memory consumption over fragments processed',
                                    labels={'x': 'fragments processed / count', 'y': 'max. consumed memory / MiB'},
                                    color_discrete_sequence=color_sequence)
    # above: use 'lowess' for LOESS interpolated curve; this is linear trend line
    memory_scatter_fig.update_layout(width=1200, height=750, template='simple_white',
                                     font={'family': 'Ubuntu', 'size': 22,
                                           'color': 'rgb(20, 20, 20)'},
                                     title={'text': 'Memory footprint',
                                            'font': {'family': 'Ubuntu', 'size': 36,
                                                     'color': 'rgb(20, 20, 20)'},
                                            'xanchor': 'center', 'yanchor': 'middle', 'x': 0.5},
                                     legend=dict(orientation='h', xanchor='center', yanchor='top',
                                                 x=0.5, y=-0.18, title='', font=dict(family='Ubuntu', size=18)))
    memory_scatter_fig.update_xaxes(showgrid=True, gridcolor='lightgrey', zerolinecolor='rgb(20,20,20)')
    memory_scatter_fig.update_yaxes(showgrid=True, gridcolor='lightgrey', zerolinecolor='rgb(20,20,20)')
    memory_scatter_fig.update_traces(marker=dict(size=16, line=dict(width=2)))  # , opacity=0.5
    #memory_scatter_fig.show()
    memory_scatter_fig.write_image(out_path / f'GCparagon_memory_consumption_presets_SSD-ref.png')
    memory_scatter_fig = px.scatter(x=performance_metrics['processed_fragments'],
                                    y=performance_metrics['memory_consumption'],  trendline='ols',
                                    labels={'x': 'fragments processed / count', 'y': 'max. consumed memory / MiB'},
                                    color_discrete_sequence=color_sequence)
    # above: use 'lowess' for LOESS interpolated curve; this is linear trend line
    memory_scatter_fig.update_layout(width=1125, height=900, template='simple_white',
                                     annotations=[{"x": int((max(performance_metrics['processed_fragments']) -
                                                   min(performance_metrics['processed_fragments'])) * 0.68),
                                                   "y": min(performance_metrics['memory_consumption']) +
                                                   (max(performance_metrics['memory_consumption']) -
                                                    min(performance_metrics['memory_consumption'])) * 0.95,
                                                   "text": sci_notation(
                                                       p_value_memory, prefix=f"Pearson's R={r_value_memory_str} (p=",
                                                       decimal_digits=2, postfix=r')'),
                                                   "font": {"size": 18}, "showarrow": False, 'bgcolor': 'white'}],
                                     font={'family': 'Ubuntu', 'size': 22,
                                           'color': 'rgb(20, 20, 20)'},
                                     title={'text': 'Maximum memory consumption over fragments processed',
                                            'font': {'family': 'Ubuntu', 'size': 36,
                                                     'color': 'rgb(20, 20, 20)'},
                                            'xanchor': 'center', 'yanchor': 'middle', 'x': 0.5},
                                     legend=dict(orientation='h', xanchor='center', yanchor='top',
                                                 x=0.5, y=-0.18, title='', font=dict(family='Ubuntu', size=18)))
    memory_scatter_fig.update_xaxes(showgrid=True, gridcolor='lightgrey', zerolinecolor='rgb(20,20,20)')
    memory_scatter_fig.update_yaxes(showgrid=True, gridcolor='lightgrey', zerolinecolor='rgb(20,20,20)')
    memory_scatter_fig.update_traces(marker=dict(size=16, line=dict(width=2)))  # , opacity=0.5
    #memory_scatter_fig.show()
    memory_scatter_fig.write_image(out_path / f'GCparagon_memory_consumption_linRegress_SSD-ref.png')
