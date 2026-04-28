
PREFIX     := /usr/local

SRC_DIR    := src
BUILD_DIR  := build
SCRIPT_DIR := scripts
SUB_DIR    := submodules
BIN_DIR    := $(BUILD_DIR)/bin
LIB_DIR    := $(BUILD_DIR)/lib

ECHO       := echo
PYTHON     := $(filter /%,$(shell /bin/sh -c 'type python'))
INSTALL    := $(filter /%,$(shell /bin/sh -c 'type install'))
MKDIR      := $(filter /%,$(shell /bin/sh -c 'type mkdir'))
AWK        := $(filter /%,$(shell /bin/sh -c 'type awk'))
CAT        := $(filter /%,$(shell /bin/sh -c 'type cat'))
SED        := $(filter /%,$(shell /bin/sh -c 'type sed'))
CP         := $(filter /%,$(shell /bin/sh -c 'type cp'))
RM         := $(filter /%,$(shell /bin/sh -c 'type rm'))
GIT        := $(filter /%,$(shell /bin/sh -c 'type git'))

GIT_SUBUPDATE = $(GIT) submodule update --init --recursive
GIT_CHECKOUT  = $(GIT) checkout

CP_R        = $(CP) -R
RM_R        = $(RM) -r

PYTHON_VER := $(shell $(PYTHON) --version 2>&1 | awk '{if (/Python/) {split($$2,v,".");print "python"v[1]"."v[2]}}')
INSTALL_REG = $(INSTALL) -m 644 -p
INSTALL_DIR = $(INSTALL) -m 755 -d
INSTALL_EXE = $(INSTALL) -m 755 -p
INSTALL_LIB = $(CP_R) -a
MKDIR_P     = $(MKDIR) -p


PACKAGE    := OGTK
VERSION    := $(shell $(GIT) describe --long --tags --always)
CONTACT    := https:\/\/github.com\/JGI-Bioinformatics\/ogtk
LICENSE    := LICENSE

BIN_TARGETS = \
	$(BIN_DIR)/add-singleton-orthogroups \
	$(BIN_DIR)/assign-BLG \
	$(BIN_DIR)/assign-CLG \
	$(BIN_DIR)/assign-hog-orthogroups \
	$(BIN_DIR)/assign-JVLG \
	$(BIN_DIR)/assign-locus-labels \
	$(BIN_DIR)/cluster-orthogroups \
	$(BIN_DIR)/count-orthogroups \
	$(BIN_DIR)/filter-orthogroups \
	$(BIN_DIR)/freq-chr-pairs-ava \
	$(BIN_DIR)/freq-chr-pairs-two \
	$(BIN_DIR)/join-orthogroups \
	$(BIN_DIR)/parse-orthogroups \
	$(BIN_DIR)/plot-color-legend \
	$(BIN_DIR)/plot-freq-chr-pairs \
	$(BIN_DIR)/plot-genomic-intervals \
	$(BIN_DIR)/plot-label-density \
	$(BIN_DIR)/propagate-LGs \
	$(BIN_DIR)/rename-assembly-report \
	$(BIN_DIR)/score-orthogroups \
	$(BIN_DIR)/select-orthogroups \
	$(BIN_DIR)/update-training-orthogroups \
	$(BIN_DIR)/upset-orthogroups

LIB_TARGETS = \
	$(LIB_DIR)/og/__init__.py \
	$(LIB_DIR)/og/constants.py \
	$(LIB_DIR)/og/core/members.py \
	$(LIB_DIR)/og/core/parsers/assembly_report.py \
	$(LIB_DIR)/og/core/parsers/bed.py \
	$(LIB_DIR)/og/core/parsers/config.py \
	$(LIB_DIR)/og/core/parsers/newick.py \
	$(LIB_DIR)/og/core/parsers/orthogroups.py \
	$(LIB_DIR)/og/core/parsers/tsv.py \
	$(LIB_DIR)/og/core/trees.py \
	$(LIB_DIR)/og/core/utils.py

SUB_TARGETS = \
	$(LIB_DIR)/og/core/compression


.SUFFIXES: .py .sh .R

.PHONY: install activate clean

all: build activate

$(BIN_DIR): 
	@$(MKDIR_P) $@

$(LIB_DIR): 
	@$(MKDIR_P) $@

$(BIN_DIR)/%: $(SCRIPT_DIR)/%.py
	@$(AWK) 'BEGIN{print "#!/usr/bin/env python3"} {print "#",$$0}' $(LICENSE) | $(CAT) - $< | \
		$(SED) "s/__PACKAGE_NAME__/$(PACKAGE)/;s/__PACKAGE_VERSION__/$(VERSION)/;s/__PACKAGE_CONTACT__/$(CONTACT)/" >$@

$(BIN_DIR)/%: $(SCRIPT_DIR)/%.sh
	@$(AWK) 'BEGIN{print "#!/usr/bin/env bash"} {print "#",$$0}' $(LICENSE) | $(CAT) - $< | \
		$(SED) "s/__PACKAGE_NAME__/$(PACKAGE)/;s/__PACKAGE_VERSION__/$(VERSION)/;s/__PACKAGE_CONTACT__/$(CONTACT)/" >$@

$(BIN_DIR)/%: $(SCRIPT_DIR)/%.R
	@$(AWK) 'BEGIN{print "#!/usr/bin/env Rscript"} {print "#",$$0}' $(LICENSE) | $(CAT) - $< | \
		$(SED) "s/__PACKAGE_NAME__/$(PACKAGE)/;s/__PACKAGE_VERSION__/$(VERSION)/;s/__PACKAGE_CONTACT__/$(CONTACT)/" >$@

$(LIB_DIR)/%: $(SRC_DIR)/%
	@$(MKDIR_P) $(@D)
	@$(AWK) 'BEGIN{print "#!/usr/bin/env python3"} {print "#",$$0}' $(LICENSE) | $(CAT) - $< | \
		$(SED) "s/__PACKAGE_NAME__/$(PACKAGE)/;s/__PACKAGE_VERSION__/$(VERSION)/;s/__PACKAGE_CONTACT__/$(CONTACT)/" >$@

$(LIB_DIR)/og/core/%: $(SUB_DIR)/%/src/%
	@$(CP_R) $(SUB_DIR)/$*/src/$* $(LIB_DIR)/og/core

$(SUB_DIR)/%/src/%:
	$(GIT_SUBUPDATE) $<

build: $(LIB_DIR) $(LIB_TARGETS) $(SUB_TARGETS) $(BIN_DIR) $(BIN_TARGETS)


activate:
	@$(ECHO) 'export PYTHONPATH="$(PREFIX)/lib/$(PYTHON_VER)/site-packages:$$PYTHONPATH";' >activate
	@$(ECHO) 'export PATH="$(PREFIX)/bin:$$PATH";' >>activate
	@$(ECHO) '#setenv PYTHONPATH "$(PREFIX)/lib/$(PYTHON_VER)/site-packages:$$PYTHONPATH";' >>activate
	@$(ECHO) '#setenv PATH "$(PREFIX)/bin:$$PATH";' >>activate


install: all
	$(INSTALL_DIR) $(PREFIX)/bin
	$(INSTALL_DIR) $(PREFIX)/lib/$(PYTHON_VER)/site-packages
	$(INSTALL_EXE) $(BIN_DIR)/* $(PREFIX)/bin
	$(INSTALL_LIB) $(LIB_DIR)/* $(PREFIX)/lib/$(PYTHON_VER)/site-packages
	$(INSTALL_REG) activate $(PREFIX)/


clean:
	-$(RM_R) $(BUILD_DIR)
