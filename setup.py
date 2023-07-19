from setuptools import setup, find_packages

setup(
    name='colloid-model',
    version='1.0.0',
    description='colloid-model is used to explore theoretical colloid performance.',
    url='https://github.com/lvenneri/Colloid-model',
    author='Lorenzo Venneri',
    author_email='lorenzo.venneri@gmail.com',
    license='MIT',
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'colloid-model': [
                    ],
    },
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
        'matplotlib',
        "openpyxl==3.1.0",
        "xlsxwriter",
        'dill',
        'plotly',
        'streamlit',
        'latex',


    ],
    entry_points={},

    classifiers=[
        'Development Status :: Concept',
        'Intended Audience :: Science/Research',
        'License :: NONE :: NONE',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.9',
    ],
)
