import os
import re
from setuptools import setup, find_packages

path_to_pyproject = 'pyproject.toml'


def read_specs(path: str):
    specs = dict()
    with open(path, 'r') as fin:
        for line in fin:
            if line.startswith('name'):
                specs['name'] = re.match('name = "(.+)"', line).group(1)
            elif line.startswith('version'):
                specs['version'] = re.match('version = "(.+)"', line).group(1)
            elif line.startswith('description'):
                specs['description'] = re.match('description = "(.+)"', line).group(1)
            elif line.startswith('requires-python'):
                specs['python_requires'] = re.match('requires-python = "(.+)"', line).group(1)
                break
    return specs


def read_requirements(path: str):
    with open(path, encoding="utf-8") as file:
        requirements = []
        extra_index_urls = []
        for line in file:
            line = line.strip()
            if line:
                if line.startswith('-f'):
                    extra_index_urls.append(line.split()[1])
                elif line.startswith('#'):
                    continue
                else:
                    requirements.append(line)
        return requirements, extra_index_urls


requirements, _ = read_requirements("requirements.txt")

scripts = []
if os.path.exists("scripts"):
    for file in os.listdir("scripts"):
        if file.endswith(".py") and not file.startswith("_"):
            fp = os.path.join("scripts", file)
            with open(fp) as fin:
                if "python" in fin.readline():
                    scripts.append(fp)

specs = read_specs(path_to_pyproject)

setup(
    author="kpotoh",
    author_email="axepon@mail.ru",
    url="https://github.com/mitoclub/PyMutSpec",
    license="MIT",
    install_requires=requirements,
    scripts=scripts,
    packages=find_packages(),
    package_data={'pymutspec': ['utils/configs/log_settings.yaml']},
    **specs
)
