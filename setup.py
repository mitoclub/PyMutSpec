from setuptools import setup, find_packages

extra_index_urls = []
packages = []

with open("requirements.txt", encoding="utf-8") as file:
    for line in map(str.strip, file):
        if line:
            if line.startswith("-f"):
                extra_index_urls.append(line.split()[1])
            else:
                packages.append(line)

setup(
    name="mutspec_utils",
    version="0.0.2",
    author="kpotoh, genkvg",
    author_email="None",
    license="MIT",
    install_requires=packages,
    dependency_links=extra_index_urls,
    packages=find_packages(),
    python_requires=">=3.8",
)
