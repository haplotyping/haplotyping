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
    packages=setuptools.find_packages(),
    install_requires=[
      "graphviz >= 0.20.1",
      "pandas >= 1.5.2",
      "networkit >= 10.0",
      "h5py >= 3.7.0",
      "tables >= 3.7.0",
      "pyahocorasick >= 1.4.4",
      "flask >= 2.2.2",
      "flask-restx >= 1.0.3",
      "waitress >= 2.1.2",
      "flask-caching >= 2.0.1",
      "openpyxl >= 3.0.10",
      "frictionless >= 4.40.8"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    test_suite="tests",
    include_package_data=True,
    package_data={"": ["data/schema/*.json"]},
)
