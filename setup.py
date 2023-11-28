from setuptools import setup, find_packages
import os
import shutil


#remove the dist folder first if exists
if os.path.exists("dist"):
	shutil.rmtree("dist")

def readme():
	this_directory = os.path.abspath(os.path.dirname(__file__))
	with open(os.path.join(this_directory, 'README.md'), encoding='latin-1') as f:
		long_description = f.read()
		return(long_description)
	# with open('README.rst') as f:
	# 	return(f.read())

VERSION = '1.1.1'

def write_version_py(filename='SigProfilerClusters/version.py'):
	# Copied from numpy setup.py
	cnt = """
# THIS FILE IS GENERATED FROM SIGPROFILECLUSTERS SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
	
	"""
	fh = open(filename, 'w')
	fh.write(cnt % {'version': VERSION,})
	fh.close()

write_version_py()

setup(name='SigProfilerClusters',
		version=VERSION,
		description='SigProfilerClusters tool',
		long_description= readme(),
		long_description_content_type='text/markdown',
		url='',
		author='Erik Bergstrom',
		author_email='ebergstr@eng.ucsd.edu',
		license='UCSD',
		packages=find_packages(),#['SigProfilerMatrixGenerator'],
		python_requires=">=3.8",
		install_requires =[
			"matplotlib>=3.3.0,<3.5",
			"sigProfilerPlotting>=1.3.18",
			"SigProfilerMatrixGenerator>=1.2.19",
			"SigProfilerSimulator>=1.1.5",
			"statsmodels>=0.9.0",
			"scipy>=1.1.0",
			"pandas>=0.23.4",
			"numpy>= 1.17.3",
			"SigProfilerExtractor>=1.0.17", 
			"seaborn>=0.13.0"],
		include_package_data=True,
		zip_safe=False)
