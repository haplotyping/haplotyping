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
      "pandas >= 2.0.3",
      "numpy >= 1.25.1",
      "networkit >= 10.0",
      "h5py >= 3.9.0",
      "tables >= 3.8.0",
      "pyahocorasick >= 2.0.0",
      "flask >= 2.2.2",
      "flask-restx >= 1.0.3",
      "waitress >= 2.1.2",
      "flask-caching >= 2.0.1",
      "openpyxl >= 3.0.10",
      "frictionless == 5.10.5",
      "psutil >= 5.9.5",
      "metis >= 0.2a5"
    ],
    test_requires=[
      "pytest >= 7.4.0",
      "pytest-ordering >= 0.6"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
    test_suite="tests",
    include_package_data=True,
    package_data={"": ["data/schema/*.json"]},
)
