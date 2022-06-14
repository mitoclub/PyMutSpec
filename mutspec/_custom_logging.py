import os
import logging
import logging.config

import yaml

DEFAULT_PATH_TO_LOGCONF = "./configs/log_settings.yaml"

with open(os.environ.get("LOG_CONFIG", DEFAULT_PATH_TO_LOGCONF), "rb") as file:
    log_config = yaml.safe_load(file)
    logging.config.dictConfig(log_config)
    del log_config

logger = logging.getLogger('MutSpecCalc')
