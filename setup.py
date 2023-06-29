import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="alf",
    version="3.1.3",
    author="Ryan Lee Hayes",
    author_email="rlhayes100@gmail.com",
    description="Adaptive Landscape Flattening scripts for Multisite lambda Dynamics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ryanleehayes/alf",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    # package_dir={"": "src"},
    # packages=setuptools.find_packages(where="src"),
    packages=['alf'],
    install_requires=['numpy','scipy','MDAnalysis'],
    python_requires=">=2.7",
    # To add external library like wham:
    # https://stackoverflow.com/questions/47360113/compile-c-library-on-pip-install
    # https://stackoverflow.com/questions/42585210/extending-setuptools-extension-to-use-cmake-in-setup-py/48015772#48015772
)
