#!/usr/bin/env python3

import sys
import pathlib
import logging


def set_up_logging(logfile_path: str, dev=False, verbose=False, logger_name=None) -> logging.Logger:
    logging.raiseExceptions = dev
    # concerning level:
    # "Logging messages which are less severe than level will be ignored; logging messages which have severity level or
    #  higher will be emitted by whichever handler or handlers service this logger, unless a handler’s level has been
    #  set to a higher severity level than level."
    log_path = pathlib.Path(logfile_path)
    root_logger = logging.getLogger(name='GCparagon_dev' if logger_name is None else logger_name)
    root_logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    formatter = logging.Formatter('| %(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler = logging.FileHandler(logfile_path, encoding='utf-8')
    file_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    file_handler.setFormatter(formatter)
    file_handler.set_name('dd_hndl')
    root_logger.addHandler(file_handler)
    stdout_handler = logging.StreamHandler(stream=sys.stdout)
    stdout_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    stdout_handler.setStream(sys.stdout)
    stdout_handler.setFormatter(formatter)
    stdout_handler.set_name('zln_hndl')
    root_logger.addHandler(stdout_handler)
    error_handler = logging.FileHandler(log_path.parent / (log_path.stem + '.err'), encoding='utf-8')
    error_handler.setLevel(logging.ERROR)
    error_handler.setFormatter(formatter)
    error_handler.set_name('fhl_hndl')
    root_logger.addHandler(error_handler)
    return root_logger


def gib_cmd_logger() -> logging.Logger:
    any_logger = logging.getLogger(name='GCparagon_dev')
    any_logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('| %(asctime)s - %(name)s - %(levelname)s - %(message)s')
    stdout_handler = logging.StreamHandler(stream=sys.stdout)
    stdout_handler.setLevel(logging.DEBUG)
    stdout_handler.setStream(sys.stdout)
    stdout_handler.setFormatter(formatter)
    stdout_handler.set_name('zln_hndl')
    any_logger.addHandler(stdout_handler)
    return any_logger


def log(message: str, log_level: int, i_log_with: logging.Logger, flush=True, close_handlers=False):
    """
    :param message:
    :param log_level:
    :param i_log_with:
    :param flush:
    :param close_handlers:
    :return:
    """
    match log_level:
        case logging.NOTSET:
            i_log_with.info(message)
        case logging.DEBUG:
            i_log_with.debug(message)
        case logging.INFO:
            i_log_with.info(message)
        case logging.WARNING:
            i_log_with.warning(message)
        case logging.ERROR:
            i_log_with.error(message)
        case logging.CRITICAL:
            i_log_with.critical(message)
    if flush and isinstance(i_log_with, logging.Logger):  # do nothing if i_log_with is undefined
        for hdlr in i_log_with.handlers:
            hdlr.flush()
    if close_handlers:
        for hdlr in i_log_with.handlers:
            hdlr.flush()
            hdlr.close()


def main():
    print("WARNING: no tests defined for gc_logging.py!")  # run your logging tests here


if __name__ == '__main__':
    sys.exit(main())
