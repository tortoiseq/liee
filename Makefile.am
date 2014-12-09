AUTOMAKE_OPTIONS = foreign
SUBDIRS = liee

buildlib:
	cd liee/lib && make buildlib
.PHONY: buildlib

all: buildlib

.DEFAULT_GOAL := all
