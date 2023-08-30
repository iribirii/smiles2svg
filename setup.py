from setuptools import setup, find_packages
version = "0.1.1"

with open('./README.md', 'r') as f:
    long_description = f.read()

setup(
        name="smiles2svg",
        packages=find_packages(exclude=["tests"]),
        version=version,
        license="MIT",
        description="Generate 2D molecular graphs from SMILES codes",
        long_description=long_description,
        long_description_content_type="text/markdown",
        author="Inigo Iribarren Aguirre",
        author_email="innigo.iribarren@gmail.com",
        keywords=[
            "computational chemistr",
            "smiles",
            "SMILES",
            "molecular graphs",
            "molecules",
            "rdkit",
            "2d molecules"
            ],
        url='https://github.com/iribirii/smiles2svg',
        classifiers=[
            "Development Status :: 3 - Beta",
            "Intended Audience :: Chemists",
            "Topic :: Molecular graphs",
            "License :: MIT License",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: 3.11"
            ],
        install_requires=[
            "CairoSVG",
            "numpy",
            "rdkit",
            "svgwrite"
            ],
        python_requires=">=3.0",
        include_package_data=True
        )
