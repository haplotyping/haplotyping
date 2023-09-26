import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
    
with open("haplotyping/_version.py", "r") as fh:
    exec(fh.read())

setuptools.setup(
    name="haplotyping", 
    version=__version__,
    author="Matthijs Brouwer",
    author_email="matthijs.brouwer@wur.nl",
    description="Haplotyping",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/haplotyping/haplotyping",
    packages=setuptools.find_packages(exclude=["tests","tests.*"]),
    install_requires=[
      "graphviz >= 0.20.1",
      "pandas >= 2.1.0",
      "numpy >= 1.25.2",
      "networkit >= 10.1",
      "h5py >= 3.9.0",
      "tables >= 3.8.0",
      "pyahocorasick >= 2.0.0",
      "Flask >= 2.3.3",
      "flask-restx >= 1.1.0",
      "waitress >= 2.1.2",
      "Flask-caching >= 2.0.2",
      "openpyxl >= 3.1.2",
      "frictionless == 5.10.5",
      "psutil >= 5.9.5",
      "metis >= 0.2a5",
      "networkit >= 10.1",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.11",
    include_package_data=True,
    package_data={"": ["data/schema/*.json"]},
)
