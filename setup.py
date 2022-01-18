from setuptools import setup, find_packages



with open("requirements.txt", "rt", encoding="utf-8") as fh:
requirements = [line.strip() for line in fh.readlines()]



setup(
	name="tools",
	version='0.0.1',
	author="os334",
	author_email="os4@sanger.ac.uk",
	description="tools",
	long_description_content_type="text/markdown",
	url="https://github.com/os334/tools/",
	packages=find_packages(),
	setup_requires=["setuptools>=45", "wheel"],
	install_requires=requirements,
	include_package_data=True,
	classifiers=[
		"Programming Language :: Python :: 3.7",
		"Programming Language :: Python :: 3.8",
		"Programming Language :: Python :: 3.9",
	],
	zip_safe=False
)