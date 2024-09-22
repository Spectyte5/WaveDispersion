from setuptools import setup

setup(
    name="WaveDispersion",
    version="0.1.0",
    description="desc",
    author="Rafal Wrobel",
    author_email="rwrobel@student.agh.edu.pl",
    python_requires=">=3.8",
    install_requires=[
        "scipy>=1.11.4",
        "matplotlib>=3.8.2"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=["wavedispersion"],
)