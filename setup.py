import setuptools

with open('README.md', 'r') as f:
    long_description = f.read()
with open('requirements.txt', 'r') as f:
    requirements = f.read().strip('\n').split('\n')

entry_points = {
    'console_scripts': [
        'spice_stew=spice_stew.spice_stew:main',
        ]
    }

setuptools.setup(
    name='spice_stew',
    version='2022.05.19',
    author='Gabriel Pelouze',
    author_email='gabriel.pelouze@universite-paris-saclay.fr',
    description='Correct the pointing of SOLO/SPICE with SPICE kernels',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/gpelouze/spice_stew',
    entry_points=entry_points,
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
    install_requires=requirements,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
)
