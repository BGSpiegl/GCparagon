#!/usr/bin/env python3

import sys
from pathlib import Path
import logging
from typing import Union as OneOf


def set_up_logging(logfile_path: OneOf[str, Path], dev=False, verbose=False, logger_name=None) -> str:
    """
    Use ONLY if no logger is available! For resetting log file paths, use the 'set_new_log_paths()' function below
    :param logfile_path:
    :param dev:
    :param verbose:
    :param logger_name:
    :return:
    """
    logging.raiseExceptions = dev
    # concerning level:
    # "Logging messages which are less severe than level will be ignored; logging messages which have severity level or
    #  higher will be emitted by whichever handler or handlers service this logger, unless a handler’s level has been
    #  set to a higher severity level than level."
    log_path = Path(logfile_path)
    logger_name = 'GCparagon-DEFAULT' if logger_name is None else logger_name
    root_logger = logging.getLogger(name=logger_name)
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
    return logger_name


def set_new_log_paths(logfile_path: OneOf[str, Path], logger_name: str, dev=False, verbose=False):
    logging.raiseExceptions = dev
    logfile_path = Path(logfile_path)
    target_logger = logging.getLogger(logger_name)
    # flush and remove old handlers
    for handler in target_logger.handlers:
        handler.flush()
        handler.close()
        target_logger.removeHandler(handler)
    # create new log paths and new commandline handler
    formatter = logging.Formatter('| %(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler = logging.FileHandler(logfile_path, encoding='utf-8')
    file_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    file_handler.setFormatter(formatter)
    file_handler.set_name('dd_hndl')
    target_logger.addHandler(file_handler)
    stdout_handler = logging.StreamHandler(stream=sys.stdout)
    stdout_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    stdout_handler.setStream(sys.stdout)
    stdout_handler.setFormatter(formatter)
    stdout_handler.set_name('zln_hndl')
    target_logger.addHandler(stdout_handler)
    error_handler = logging.FileHandler(logfile_path.parent / (logfile_path.stem + '.err'), encoding='utf-8')
    error_handler.setLevel(logging.ERROR)
    error_handler.setFormatter(formatter)
    error_handler.set_name('fhl_hndl')
    target_logger.addHandler(error_handler)


def gib_cmd_logger() -> str:
    logger_name = 'GCparagon-CMDlogger'
    any_logger = logging.getLogger(name=logger_name)
    any_logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('| %(asctime)s - %(name)s - %(levelname)s - %(message)s')
    stdout_handler = logging.StreamHandler(stream=sys.stdout)
    stdout_handler.setLevel(logging.DEBUG)
    stdout_handler.setStream(sys.stdout)
    stdout_handler.setFormatter(formatter)
    stdout_handler.set_name('zln_hndl')
    any_logger.addHandler(stdout_handler)
    return logger_name


def log(message: str, log_level: int, logger_name: str, flush=True, close_handlers=False):
    """
    :param message:
    :param log_level:
    :param logger_name: name of logger; will be used with logging.getLogger()
    :param flush:
    :param close_handlers:
    :return:
    """
    current_logger = logging.getLogger(logger_name)
    match log_level:
        case logging.NOTSET:
            current_logger.info(message)
        case logging.DEBUG:
            current_logger.debug(message)
        case logging.INFO:
            current_logger.info(message)
        case logging.WARNING:
            current_logger.warning(message)
        case logging.ERROR:
            current_logger.error(message)
        case logging.CRITICAL:
            current_logger.critical(message)
    if flush and isinstance(current_logger, logging.Logger):  # do nothing if logger_name is undefined
        for hdlr in current_logger.handlers:
            hdlr.flush()
    if close_handlers:
        for hdlr in current_logger.handlers:
            hdlr.flush()
            hdlr.close()


def main():
    print("WARNING: no tests defined for gc_logging.py!")  # run your logging tests here


if __name__ == '__main__':
    sys.exit(main())
