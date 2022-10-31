import os
import logging
import logging.config

import yaml

DEFAULT_PATH_TO_LOGCONF = os.path.join(os.path.dirname(__file__), "configs/log_settings.yaml")


def load_logger(path=DEFAULT_PATH_TO_LOGCONF, stream_level: str = None, filename=None):
    with open(path, "r") as file:
        log_config = yaml.safe_load(file)
        if stream_level is not None:
            log_config["handlers"]["stream_handler"]["level"] = stream_level
        if filename is not None:
            log_config["handlers"]["full_file_handler"]["filename"] = filename
        logging.config.dictConfig(log_config)
    logger = logging.getLogger('MutSpecCalc')
    return logger
