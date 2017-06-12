from distutils.core import setup

setup(
    name='ipas',
    version='0.1dev',
    packages=['ipas'],
    install_requires=[
        "shapely >= 1.6",
        "numpy >= 1.13",
    ]
)
