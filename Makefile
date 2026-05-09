
PREFIX     ?= /usr/local
INSTALL_PATH ?= $(PREFIX)/lib/$(PYTHON_VERSION)/site-packages

CURR_PATH   = $(shell pwd)

SRC_DIR    := $(CURR_PATH)/src
BUILD_DIR  := $(CURR_PATH)/build
SCRIPT_DIR := $(CURR_PATH)/scripts
SUB_DIR    := $(CURR_PATH)/submodules
BIN_DIR    := $(BUILD_DIR)/bin
LIB_DIR    := $(BUILD_DIR)/lib

ECHO       := $(shell which echo 2>/dev/null)
PYTHON     := $(shell which python 2>/dev/null)
INSTALL    := $(shell which install 2>/dev/null)
MKDIR      := $(shell which mkdir 2>/dev/null)
AWK        := $(shell which awk 2>/dev/null)
CAT        := $(shell which cat 2>/dev/null)
SED        := $(shell which sed 2>/dev/null)
CP         := $(shell which cp 2>/dev/null)
CP_R        = $(CP) -R
RM         := $(shell which rm 2>/dev/null)
RM_R        = $(RM) -r
GIT        := $(shell which git 2>/dev/null)

INSTALL_DIR = $(INSTALL) -m 755 -d
INSTALL_EXE = $(INSTALL) -m 755 -p
INSTALL_LIB = $(CP_R) -a
INSTALL_REG = $(INSTALL) -m 644 -p
MKDIR_P     = $(MKDIR) -p

GIT_SUBUPDATE = $(GIT) submodule update --init --recursive
GIT_CHECKOUT  = $(GIT) checkout

PROJECT    := OGTK
LIBRARY    := og
VERSION    := $(shell $(GIT) describe --long --tags --always)
CONTACT    := https:\/\/github.com\/JGI-Bioinformatics\/ogtk
LICENSE    := LICENSE


ifneq ($(shell which python3),)
PYTHON     := $(shell which python3)
else ifneq ($(shell which python),)
PYTHON     := $(shell which python)
else
$(error "Python interpreter not found. Please install Python and ensure it is accessible via PATH.")
endif

PYTHON_VERSION := $(shell $(PYTHON) --version 2>&1 | awk '{if (/Python/) {split($$2,v,".");print "python"v[1]"."v[2]}}')



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
	$(LIB_DIR)/$(LIBRARY)/__init__.py \
	$(LIB_DIR)/$(LIBRARY)/constants.py \
	$(LIB_DIR)/$(LIBRARY)/core/members.py \
	$(LIB_DIR)/$(LIBRARY)/core/parsers/assembly_report.py \
	$(LIB_DIR)/$(LIBRARY)/core/parsers/bed.py \
	$(LIB_DIR)/$(LIBRARY)/core/parsers/config.py \
	$(LIB_DIR)/$(LIBRARY)/core/parsers/newick.py \
	$(LIB_DIR)/$(LIBRARY)/core/parsers/orthogroups.py \
	$(LIB_DIR)/$(LIBRARY)/core/parsers/tsv.py \
	$(LIB_DIR)/$(LIBRARY)/core/trees.py \
	$(LIB_DIR)/$(LIBRARY)/core/utils.py

SUB_TARGETS = \
	$(LIB_DIR)/bgzip.py \
	$(LIB_DIR)/$(LIBRARY)/core/strand.py \
	$(LIB_DIR)/$(LIBRARY)/core/compression \
	$(LIB_DIR)/$(LIBRARY)/core/intervals \


.SUFFIXES:
.SUFFIXES: .py .sh .R

.PHONY: all install activate clean

all: $(LIB_DIR) $(BIN_DIR) $(LIB_TARGETS) $(SUB_TARGETS) $(BIN_TARGETS) activate

$(BUILD_DIR):
	@$(MKDIR_P) $@

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


$(LIB_DIR)/$(LIBRARY)/core/%: $(LIB_DIR)/$(LIBRARY) $(SUB_DIR)/%/src/%
	make -C $(SUB_DIR)/$(subst .py,,$*) install INSTALL_PATH=$(@D)

$(LIB_DIR)/bgzip.py: $(SUB_DIR)/compression/src/bgzip.py
	make -C $(SUB_DIR)/compression install-bgzip INSTALL_PATH=$(@D)

$(SUB_DIR)/%/src/%:
	$(GIT_SUBUPDATE) $<


activate:
	@$(ECHO) 'export PYTHONPATH="$(INSTALL_PATH):$$PYTHONPATH";' >activate
	@$(ECHO) 'export PATH="$(PREFIX)/bin:$$PATH";' >>activate
	@$(ECHO) '#setenv PYTHONPATH "$(INSTALL_PATH):$$PYTHONPATH";' >>activate
	@$(ECHO) '#setenv PATH "$(PREFIX)/bin:$$PATH";' >>activate


install: all
	$(INSTALL_DIR) $(PREFIX)/bin
	$(INSTALL_DIR) $(INSTALL_PATH)
	$(INSTALL_EXE) $(BIN_DIR)/* $(PREFIX)/bin
	$(INSTALL_LIB) $(LIB_DIR)/* $(INSTALL_PATH)
	$(INSTALL_REG) activate $(PREFIX)/


clean:
	-$(RM_R) $(BUILD_DIR)
