from distutils.core import setup, Extension
from Cython.Build import cythonize



library_dirs = ['prelib']

include_dirs = ['SeqLib','SeqLib/htslib',
                'SGA/SGA', 'SGA/StringGraph', 'SGA/Algorithm', 'SGA/SuffixTools', 'SGA/Bigraph', 'SGA/Util', 'SGA/SQG']

setup(ext_modules=cythonize(Extension(
   "sv_helper", # the extension name
   sources=["src/c_sv_helper.pyx", "src/sv_helper.cpp",
            "src/svabaASQG.cpp", "src/svabaOverlapAlgorithm.cpp"],
   language="c++",
   library_dirs=library_dirs,
   include_dirs=include_dirs,
   libraries=['seqlib','bwa', 'hts',
              'sga', 'stringgraph', 'algorithm', 'suffixtools', 'bigraph', 'util', 'sqg',
              'z', 'pthread', 'bz2', 'lzma',],
)))