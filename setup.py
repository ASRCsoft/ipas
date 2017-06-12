from distutils.core import setup

setup(
    name='ipas',
    version='0.1dev',
    packages=['ipas'],
    install_requires=[
        "shapely >= 1.5",
        "numpy >= 1.12",
    ]
)
