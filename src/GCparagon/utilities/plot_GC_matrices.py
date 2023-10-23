#!/usr/bin/env python3

import scipy
import logging
import numpy as np
import plotly.express as px
from pathlib import Path
from plotly import graph_objs as go
from pandas import DataFrame as pd_DF
from scipy.ndimage import gaussian_filter
from typing import Optional
# project imports
from GCparagon.utilities.gc_logging import log

# DEFINE COLOR SCALE FOR MATRIX PLOTS HERE!!!
USE_CONTINUOUS_COLOR_SCALE = 'turbo'  # also nice: 'thermal', 'agsunset', 'rainbow', 'plotly3', or 'blackbody'
# (better for grayscale conversion) earlier default was 'viridis'; others: 'matter', 'greys'


def plot_statistic_matrices(frq_data: dict, data_id_to_show: str, in_file: str, parent_logger: Optional[str] = None,
                            output_dir=None, percentile_plot=None, sample_id=None, y_tick_label_offset=0,
                            x_tick_label_offset=0, fig_width=1800, fig_height=2000, fig_fontsize=50, show_figure=False,
                            image_formats=('png',)):
    try:
        matrix_type = list(filter(lambda a: a is not None,
                                  [mat_type if mat_type in data_id_to_show else None
                                   for mat_type in ('W_gc', 'O_gc', 'S_gc', 'Mask', 'D_gc')]))[0]  # take first found
    except IndexError:  # list index out of range
        log(message=f"received value for 'data_id_to_show' (= '{data_id_to_show}') not allowed. Must contain a string "
                    f"from 'W_gc', 'O_gc', 'S_gc', 'Mask', or 'D_gc'.", logger_name=parent_logger,
            log_level=logging.ERROR, flush=True, close_handlers=True)
        raise AttributeError
    plot_data = frq_data[data_id_to_show]
    if plot_data.empty:  # AttributeError: 'numpy.ndarray' object has no attribute 'empty'
        log(message=f"no data available for category '{data_id_to_show}'. Returning without plotting ..",
            logger_name=parent_logger, log_level=logging.WARNING)
        return
    # distill annotations and data_type
    matrix_type_mapping = {'w_gc': 'correction weights',
                           's_gc': 'expected fragments',
                           'o_gc': 'observed fragments',
                           'mask': 'weighted fragments mask',
                           'd_gc': 'difference'}
    annotation = data_id_to_show[data_id_to_show.index(matrix_type) + len(matrix_type):].strip(
        '-').strip('_').strip('-')
    if 'Mask' == matrix_type:
        figure_title = f"{sample_id if sample_id else ''} " \
                       f"{matrix_type_mapping[matrix_type.lower()]}" \
                       f"{' (' + annotation + ')' if annotation else ''}"
        fig = px.imshow(frq_data[data_id_to_show],
                        title=figure_title,
                        template="simple_white", width=fig_width, height=fig_height,
                        labels={'x': 'GC bases / bp',
                                'y': 'fragment length / bp'},
                        color_continuous_scale=['white', 'black'])
    else:
        figure_title = f"{sample_id if sample_id else ''} " \
                              f"{matrix_type_mapping[matrix_type.lower()]}" \
                              f"{' (' + annotation + ')' if annotation else ''}"
        fig = px.imshow(frq_data[data_id_to_show],
                        title=figure_title,
                        template="simple_white", width=fig_width, height=fig_height,
                        labels={'x': 'GC bases / bp',
                                'y': 'fragment length / bp'},
                        color_continuous_scale=USE_CONTINUOUS_COLOR_SCALE)
    fig.update_layout(font_family="Ubuntu",
                      font_size=fig_fontsize,
                      legend=dict(title=data_id_to_show,
                                  orientation="h", y=1, yanchor="bottom", x=0.5, xanchor="center"),
                      title={'text': figure_title,
                             'font': {'family': 'Ubuntu', 'size': fig_fontsize+4, 'color': 'rgb(20, 20, 20)'},
                             'xanchor': 'center', 'yanchor': 'middle', 'x': 0.5},
                      margin=go.layout.Margin(
                          l=5,   # left margin
                          r=5,   # right margin
                          b=5,   # bottom margin
                          t=60)  # top margin
                      )
    for idx in range(len(fig.data)):
        fig.data[idx].y = fig.data[idx].y + y_tick_label_offset
        fig.data[idx].x = fig.data[idx].x + x_tick_label_offset
    if show_figure:
        fig.show()
    in_file_path = Path(in_file)
    if output_dir is None:
        output_dir = in_file_path.parent
    for image_format in image_formats:
        out_file = Path(output_dir) / \
            ((f"{'.'.join(in_file_path.name.split('.')[:-1])}" if sample_id is None else sample_id) +
             f".{data_id_to_show}.heatmap.{image_format}")
        fig.write_image(out_file)  # requires: requests
    # below percentile plot
    if percentile_plot is not None:
        if data_id_to_show == 'R_gc':
            replacement = 1.
        else:
            replacement = 0.
        cleaned_data = remove_above_percentile(data=frq_data[data_id_to_show], percentile=percentile_plot,
                                               replacement_value=replacement)
        fig = px.imshow(cleaned_data, template="simple_white", width=fig_width, height=fig_height,
                        title=f"{sample_id if sample_id else ''} "
                              f"{matrix_type_mapping[matrix_type.lower()]}"
                              f"{' (' + annotation + ')' if annotation else ''}"
                              f"(below {percentile_plot}th percentile)",
                        color_continuous_scale=USE_CONTINUOUS_COLOR_SCALE)
        fig.update_layout(font_family="Ubuntu",
                          font_size=fig_fontsize,
                          legend=dict(title=f"{data_id_to_show} below {percentile_plot}th percentile)",
                                      orientation="h", y=1, yanchor="bottom", x=0.5, xanchor="center"))
        if show_figure:
            fig.show()
        for image_format in image_formats:
            out_file = Path(output_dir) / \
                (f"{'.'.join(in_file_path.name.split('.')[:-1])}" if sample_id is None else sample_id) + \
                f".{data_id_to_show}.{percentile_plot}ltPerc.heatmap.{image_format}"
            fig.write_image(out_file)


# post-processing functions
def remove_above_percentile(data: pd_DF, percentile: int, replacement_value=1.) -> pd_DF:
    data_percentile = np.percentile(data, percentile)
    new_data = data[:]
    new_data = new_data.where(data < data_percentile, replacement_value)
    return new_data


def limit_extreme_outliers(outliers_matrix: np.array, outliers_factor=8, parent_logger: Optional[str] = None,
                           detection_method='IQR') -> np.array:
    relevant_masked_weights = outliers_matrix[(outliers_matrix != 0.) * (outliers_matrix != 1.)].flatten()
    outliers_replaced_matrix = outliers_matrix.copy()
    if detection_method == 'IQR':  # gives lower threshold in general (more weights get limited)
        q3 = np.quantile(relevant_masked_weights, 0.75, method='linear')
        iqr = q3 - np.quantile(relevant_masked_weights, 0.25, method='linear')
        outliers_threshold = q3 + outliers_factor * iqr  # Q3 + 3SD is def. of extreme outliers; here: lenient def.
        outliers_mask = outliers_matrix > outliers_threshold
        outliers_replaced_matrix[outliers_mask] = outliers_threshold
        number_of_outlier_weights = outliers_mask.sum()
    elif detection_method == 'SD':
        outliers_threshold = float(np.median(relevant_masked_weights)) + \
                             outliers_factor * np.std(relevant_masked_weights)
        outliers_mask = outliers_matrix > outliers_threshold
        outliers_replaced_matrix[outliers_mask] = outliers_threshold
        number_of_outlier_weights = outliers_mask.sum()
    else:
        if parent_logger is None:
            print(f"ERROR - unknown outliers detection method '{detection_method}'. Must be one of 'SD', 'IQR'.")
        else:
            log(message=f"unknown outliers detection method '{detection_method}'. Must be one of 'SD', 'IQR'.",
                logger_name=parent_logger, log_level=logging.ERROR, flush=True, close_handlers=True)
        raise AttributeError
    # give user feedback
    log(message=f" i :  {outliers_factor} x{detection_method.upper()} outlier threshold was {outliers_threshold:.3} " +
                (f"(median of weights != {{0.0; 1.0}} + {outliers_factor} standard deviations). "
                 if detection_method == 'SD' else
                 f"(Q3 of weights != {{0.0; 1.0}} + {outliers_factor} * IQR). ") +
                f"Number of correction weights above {outliers_threshold:.3}: {number_of_outlier_weights:,} (= "
                f"{number_of_outlier_weights / len(relevant_masked_weights):.4%} of all weights != {{0.0; 1.0}})",
        logger_name=parent_logger, log_level=logging.INFO)
    return outliers_replaced_matrix


# postprocessing functions
def smooth_2d_gc_weights(smooth_matrix: np.array, min_flen: int, parent_logger: str, default_matrix_value=1.,
                         smoothing_kernel='gauss', smoothing_intensity=1):
    available_smoothing_kernels = ('constant', 'gauss')
    if smoothing_intensity < 1:
        log(message=f"intensity of smoothing must be positive (not zero) but was {smoothing_intensity}",
            logger_name=parent_logger, log_level=logging.ERROR, flush=True, close_handlers=True)
        raise AttributeError
    if smoothing_kernel not in available_smoothing_kernels:
        log(message=f"""unknown smoothing method '{smoothing_kernel}'. Must be one of: """
                    f"""'{"', '".join(available_smoothing_kernels)}'!""", logger_name=parent_logger,
            log_level=logging.ERROR, flush=True, close_handlers=True)
        raise AttributeError
    if smoothing_intensity < 1:
        log(message=f"intensity of smoothing must be a positive integer but was {smoothing_intensity}",
            logger_name=parent_logger, log_level=logging.ERROR, flush=True, close_handlers=True)
        raise AttributeError
    # replace non-existent attribute combinations (zeros) with ones:
    smooth_matrix[smooth_matrix == 0.] = 1.  # WARNING: might replace low-precision data points which values approach
    # zero with ones! -> use higher float precision in these cases! [ DEFAULT should be 6 digits after the comma ]
    # smooth according to input method:
    if smoothing_kernel == 'constant':  # equals a moving average operation
        window = np.ones((1 + 2*smoothing_intensity, 1 + 2*smoothing_intensity)) / \
                 ((1. + 2*smoothing_intensity) ** 2)  # kernel sum should be 1
        smoothed_matrix = scipy.signal.convolve2d(smooth_matrix, window, 'same', boundary='fill',
                                                  fillvalue=default_matrix_value)
    elif smoothing_kernel == 'gauss':
        #   sigma, standard deviations of the Gaussian filter are given for each axis (1., 1.)
        #   order=0, derivatives of a Gaussian:  An order of 0 corresponds to convolution with a Gaussian kernel.
        #   mode='reflect', The mode parameter determines how the array borders are handled, where cval is the value
        #                   when mode is equal to ‘constant’. Default is ‘reflect’.
        #   cval=0.0, Value to fill past edges of input if mode is ‘constant’. Default is 0.0.
        #   truncate=4.0, Truncate the filter at this many standard deviations. Default is 4.0
        gauss_sigma = 0.2 * smoothing_intensity + 0.3
        smoothed_matrix = gaussian_filter(input=smooth_matrix, sigma=(gauss_sigma, gauss_sigma), mode='constant',
                                          order=0, cval=1., truncate=4.)
    else:
        log(message=f"ended up in an impossible else clause due to passed unknown kernel attribute {smoothing_kernel}",
            logger_name=parent_logger, log_level=logging.ERROR, flush=True, close_handlers=True)
        raise AttributeError  # should never occur; basically just for linter
    # set back non-existing values to zero
    for f_len in range(min_flen, smooth_matrix.shape[1]):
        smoothed_matrix[f_len-min_flen][f_len+1:] = 0.
    return smoothed_matrix


if __name__ == '__main__':
    pass  # create your tests here!
