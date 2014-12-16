AUTOMAKE_OPTIONS = foreign
SUBDIRS = lib liee

buildlib:
	cd lib && make buildlib
.PHONY: buildlib

all: buildlib

.DEFAULT_GOAL := all
