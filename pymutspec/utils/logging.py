import os
import logging
import logging.config

import yaml

DEFAULT_PATH_TO_LOGCONF = os.path.join(os.path.dirname(__file__), "configs/log_settings.yaml")


def load_logger(path=None, stream_level: str = None, filename=None):
    path = path or DEFAULT_PATH_TO_LOGCONF
    with open(path, "r") as file:
        log_config = yaml.safe_load(file)
        if stream_level is not None:
            log_config["handlers"]["stream_handler"]["level"] = stream_level
        if filename is not None:
            log_config["handlers"]["full_file_handler"]["filename"] = filename
        logging.config.dictConfig(log_config)
    logger = logging.getLogger('MutSpecCalc')
    return logger


def basic_logger():
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s', '%d-%m-%y %H:%M:%S')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger
