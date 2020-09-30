AUTOMAKE_OPTIONS = foreign
SUBDIRS = SeqLib SGA/SQG SGA/Bigraph SGA/Algorithm SGA/Util SGA/StringGraph SGA/SGA SGA/SuffixTools

install:
	mkdir -p prelib && cp -f SGA/*/*.a prelib && cp -f SeqLib/*/*.a prelib && cp config.h src # copy the config to src