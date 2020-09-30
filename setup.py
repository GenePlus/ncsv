import os
import sys

from setuptools import setup, Extension
from Cython.Build import cythonize

basepath = os.path.dirname(os.path.join(os.getcwd(), __file__))
sys.path.insert(0, os.path.join(basepath, 'lib')) # for import config object
os.chdir(basepath) # into the basepath, setup must be run there

from ncsv.config import GlobalConfig as SvConfig

setup(
    name=SvConfig.PACKAGE_NAME,
    version=SvConfig.VERSION,
    description='NoahCare sv calling',
    package_dir={'': 'lib'},
    packages=[SvConfig.PACKAGE_NAME],
    install_requires=['pysam>=0.14', # use the header method, before 0.14 header is only a dict
                      'ConcurrentLogHandler==0.9.1',
                      'ujson==1.33', # use given version, high(2.0.2) show "is not JSON serializable" error
                      'intervaltree==2.1.0',
                     ],
    zip_safe=False,
    ext_modules=cythonize(Extension(
        "sv_helper", # the extension name
        sources=["extension/sv_helper/src/c_sv_helper.pyx", "extension/sv_helper/src/sv_helper.cpp",
                 "extension/sv_helper/src/svabaASQG.cpp", "extension/sv_helper/src/svabaOverlapAlgorithm.cpp"],
        language="c++",
        library_dirs=['extension/sv_helper/prelib'],
        include_dirs=['extension/sv_helper/SeqLib', 'extension/sv_helper/SeqLib/htslib',
                      'extension/sv_helper/SGA/SGA', 'extension/sv_helper/SGA/StringGraph',
                      'extension/sv_helper/SGA/Algorithm', 'extension/sv_helper/SGA/SuffixTools',
                      'extension/sv_helper/SGA/Bigraph', 'extension/sv_helper/SGA/Util',
                      'extension/sv_helper/SGA/SQG'],
        libraries=['seqlib', 'bwa', 'hts',
                   'sga', 'stringgraph', 'algorithm', 'suffixtools', 'bigraph', 'util', 'sqg',
                   'z', 'pthread', 'bz2', 'lzma',],
        extra_compile_args=['-std=c++11'],
    ))
)
