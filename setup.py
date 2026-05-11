from setuptools import setup, find_packages

setup(
    name="HERACLES",
    version="0.1.0",
    description="Hera Engine for small RNA Analysis with CLustering and Expression Signatures",
    long_description=open("README.md").read() if __import__("os").path.exists("README.md") else "",
    long_description_content_type="text/markdown",
    author="HERANOVA",
    author_email="Allen.Hu@yizhenbio.com",
    url="https://github.com/Heranova-Lifesciences/HERACLES",
    packages=find_packages(),
    install_requires=[
        "pandas>=1.3.0",
        "numpy>=1.21.0",
        "biopython>=1.79",
        "pydeseq2>=0.4.0",
        "networkx>=2.6",
        "tqdm>=4.62.0",
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "multiqc>=1.12",
        "gseapy>=1.0.0",
    ],
    entry_points={
        "console_scripts": [
            "heracles=main:main",
        ],
    },
    python_requires=">=3.10",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
)
