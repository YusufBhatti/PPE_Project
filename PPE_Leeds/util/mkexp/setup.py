from distutils.core import setup

name = 'mkexp'
version = '0.2.1dev'

setup(
    name = name,
    version = version,
    description = 'Run-script generation for earth system models',
    long_description = open('README.txt').read(),
    author = 'Karl-Hermann Wieners',
    author_email = 'karl-hermann.wieners@mpimet.mpg.de',
    url = 'http://code.zmaw.de/projects/esmenv',
    py_modules = ['configobj', 'validate', 'feedback', 'expconfig'],
    scripts = ['mkexp', 'getexp', 'rmexp', 'diffexp', 'diffpath'],
    platforms = ['Posix'],
    license = 'LICENSE.txt',
    requires = ['Jinja2(>= 2.6)']
)
