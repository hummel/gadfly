import sys
try:
    from setuptools import setup
    have_setuptools = True
except ImportError:
    from distutils.core import setup
    have_setuptools = False

def setup_pyGadget():

    setup(
        name='pyGadget',
        version='0.0',
        description='Toolkit for analyzing Gadget SPH using pandas dataframes',
        author='Jacob Hummel',
        author_email='jhummel@astro.as.utexas.edu',
        url='http://www.as.utexas.edu/~jhummel/',
        license='MIT',
        classifiers=[
            'Development Status :: 1 - Planning',
            'Intended Audience :: Science/Research',
            'Intended Audience :: Developers',
            'License :: MIT License',
            'Programming Language :: Python :: 2',
            'Topic :: Scientific/Engineering :: Astronomy',
            'Topic :: Scientific/Engineering :: Physics'
        ],
        # zip_safe=False,
        packages=['pyGadget'], # ['pyGadget', 'pyGadget.more'],
        package_dir={
            'pyGadget' : 'pyGadget',
            # 'pyGadget.more' : 'pyGadget/more',
        },
    )
    return

if __name__ == '__main__':
    setup_pyGadget()
