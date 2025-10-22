import setuptools
import platform


def install_torch(package: str, version: str = '', cpu: bool=True):
    """This is needed to make setup compatible with M1 chip."""
    cuda = "arm" not in platform.platform()
    python_version = ''.join(platform.python_version().split('.')[:2])

    if cpu:
        # https://download.pytorch.org/whl/cpu/torch-2.1.2%2Bcpu-cp311-cp311-linux_x86_64.whl
        # https://download.pytorch.org/whl/cpu/torchvision-0.16.2%2Bcpu-cp311-cp311-linux_x86_64.whl
        return ''.join([
            f'{package} @ https://download.pytorch.org/whl/',
            f'cpu/',
            f'{package}',
            f'-{version}' if version else '',
            '%2Bcpu',
            f'-cp{python_version}-cp{python_version}',
            'm' if int(python_version) <= 37 else '',
            '-linux_x86_64.whl',
        ])
    else:
        return ''.join([
            f'{package} @ https://download.pytorch.org/whl/',
            f'cu111/' if cuda else '',
            f'{package}',
            f'-{version}' if version else '',
            '%2Bcu111' if cuda else '',
            f'-cp{python_version}-cp{python_version}',
            'm' if int(python_version) <= 37 else '',
            '-linux_x86_64.whl',
        ])

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="molannotator",
    version="0.0.1",
    author="Lucas Morin",
    author_email="lum@zurich.ibm.com",
    description="A Python library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.ibm.com/LUM/graph-recognition/",
    packages=setuptools.find_packages(exclude=["tests.*", "tests"]),
    install_requires=[
        "mol-depict @ git+ssh://git@github.ibm.com/LUM/molecule-depictor.git",
        install_torch('torch', '2.1.2', cpu=True),
        install_torch('torchvision', '0.16.2', cpu=True),
        "rdkit-pypi",
        "CairoSVG",
        "streamlit",
        "streamlit-ketcher",
        "matplotlib"
    ],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        "License :: Other/Proprietary License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
    ],
    python_requires='>=3.11',
    package_data={"": ["*.json"]},
)
