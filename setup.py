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
    name="mutspec",
    version="0.0.1",
    author="kpotoh, genkvg",
    author_email="None",
    license="MIT",
    install_requires=packages,
    dependency_links=extra_index_urls,
    # package_data={"": ["data/*", "data/**/*"]},
    packages=find_packages(include="core*"),
    python_requires=">=3.8",
)
