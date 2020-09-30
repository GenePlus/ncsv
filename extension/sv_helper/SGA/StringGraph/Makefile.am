noinst_LIBRARIES = libstringgraph.a
AM_CXXFLAGS = -std=c++11
libstringgraph_a_CPPFLAGS = \
	-I$(top_srcdir)/SGA/Bigraph \
	-I$(top_srcdir)/SGA/Util \
	-I$(top_srcdir)/SGA/Thirdparty \
	-I$(top_srcdir)/SGA/Algorithm \
	-I$(top_srcdir)/SGA/SQG \
	-I$(top_srcdir)/SeqLib

libstringgraph_a_SOURCES = \
        SGUtil.cpp SGUtil.h \
        SGAlgorithms.cpp SGAlgorithms.h \
        SGVisitors.h SGVisitors.cpp \
        CompleteOverlapSet.h CompleteOverlapSet.cpp \
        RemovalAlgorithm.h RemovalAlgorithm.cpp \
		SGSearch.h SGSearch.cpp \
		GraphSearchTree.h \
		SGWalk.h SGWalk.cpp

