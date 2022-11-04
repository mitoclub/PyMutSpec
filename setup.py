import os
from setuptools import setup, find_packages

extra_index_urls = []
packages = []
scripts = []

with open("requirements.txt", encoding="utf-8") as file:
    for line in map(str.strip, file):
        if line:
            if line.startswith("-f"):
                extra_index_urls.append(line.split()[1])
            elif line.startswith("#"):
                continue
            else:
                packages.append(line)

if os.path.exists("scripts"):
    for file in os.listdir("scripts"):
        if file.endswith(".py") and not file.startswith("_"):
            fp = os.path.join("scripts", file)
            with open(fp) as fin:
                if "python" in fin.readline():
                    scripts.append(fp)

setup(
    name="mutspec_utils",
    version="0.0.3",
    author="kpotoh",
    description="Utilities for advanced analysis of mutational spectra",
    url="https://github.com/mitoclub/mutspec-utils",
    author_email="None",
    license="MIT",
    install_requires=packages,
    dependency_links=extra_index_urls,
    scripts=scripts,
    packages=find_packages(),
    # include_package_data=True,
    package_data={'mutspec_utils': ['utils/configs/log_settings.yaml']},
    python_requires=">=3.8",
)
