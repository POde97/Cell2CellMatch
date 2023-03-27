from setuptools import setup

setup(
    name='Cell2CellMatch',
    version='1.0',    
    description='functional network of cells from scRna-seq ',
    url='https://github.com/POde97/Cell2CellMatch',
    author='Paolo Odello',
    author_email='paoloodeo.o@gmail.com',
    license='MIT license',
    packages=['Cell2CellMatch'],
    install_requires=['igraph==0.10.2',
			'leidenalg',
			'networkx==2.6.3',
		      	'CellIDpy'
			]
)
