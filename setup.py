import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="haplotyping", 
    version="0.0.1",
    author="Matthijs Brouwer",
    author_email="matthijs.brouwer@wur.nl",
    description="Haplotyping",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/haplotyping/haplotyping",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    test_suite="tests",
)
