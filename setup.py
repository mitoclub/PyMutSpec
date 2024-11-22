import os
from setuptools import setup, find_packages

scripts = []
if os.path.exists("scripts"):
    for file in os.listdir("scripts"):
        if file.endswith(".py") and not file.startswith("_"):
            fp = os.path.join("scripts", file)
            scripts.append(fp)

setup(
    packages=find_packages(),
    scripts=scripts,
    package_data={'pymutspec': ['utils/configs/log_settings.yaml']},
)
