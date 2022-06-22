from setuptools import setup

setup(
    name="illumina2vcf",
    version="1.0.0",
    author="Adam Faulconbridge",
    author_email="adam@sanogenetics.com",
    packages=["illumina2vcf"],
    description="Converter to turn Illumina final report into a VCF file.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/sanogenetics/illumina2vcf",
    install_requires=open("requirements.txt").readlines(),
    extras_require={
        "dev": open("requirements-dev.txt").readlines(),
    },
)
