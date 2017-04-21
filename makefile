########################################################################
# Makefile for LiDAR Viewer, a visualization and analysis application
# for large 3D point cloud data.
# Copyright (c) 2004-2017 Oliver Kreylos
#
# This file is part of the WhyTools Build Environment.
# 
# The WhyTools Build Environment is free software; you can redistribute
# it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
# 
# The WhyTools Build Environment is distributed in the hope that it will
# be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with the WhyTools Build Environment; if not, write to the Free
# Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# 02111-1307 USA
########################################################################

# Directory containing the Vrui build system. The directory below
# matches the default Vrui installation; if Vrui's installation
# directory was changed during Vrui's installation, the directory below
# must be adapted.
VRUI_MAKEDIR := /usr/local/share/Vrui-4.2/make
ifdef DEBUG
  VRUI_MAKEDIR = $(VRUI_MAKEDIR)/debug
endif

# Base installation directory for LiDAR Viewer. If this is set to the
# default of $(PWD), LiDAR Viewer does not have to be installed to be
# run. Created executables and resources will be installed in the bin
# and share directories under the given base directory, respectively.
# Important note: Do not use ~ as an abbreviation for the user's home
# directory here; use $(HOME) instead.
INSTALLDIR := $(shell pwd)

########################################################################
# Everything below here should not have to be changed
########################################################################

# Version number for installation subdirectories. This is used to keep
# subsequent release versions of LiDAR Viewer from clobbering each
# other. The value should be identical to the major.minor version
# number found in VERSION in the root package directory.
VERSION = 2.13

# Set up resource directories: */
CONFIGDIR = etc/LidarViewer-$(VERSION)

# Include definitions for the system environment and system-provided
# packages
include $(VRUI_MAKEDIR)/SystemDefinitions
include $(VRUI_MAKEDIR)/Packages.System
include $(VRUI_MAKEDIR)/Configuration.Vrui
include $(VRUI_MAKEDIR)/Packages.Vrui

# Set installation directory structure:
EXECUTABLEINSTALLDIR = $(INSTALLDIR)/$(EXEDIR)
ETCINSTALLDIR = $(INSTALLDIR)/$(CONFIGDIR)

########################################################################
# List common packages used by all components of this project
# (Supported packages can be found in $(VRUI_MAKEDIR)/Packages.*)
########################################################################

PACKAGES = MYGEOMETRY MYMATH MYIO MYTHREADS MYMISC

########################################################################
# Specify all final targets
########################################################################

ALL = $(EXEDIR)/CalcLasRange \
      $(EXEDIR)/LidarPreprocessor \
      $(EXEDIR)/LidarSpotRemover \
      $(EXEDIR)/LidarSubtractor \
      $(EXEDIR)/LidarIlluminator \
      $(EXEDIR)/LidarColorMapper \
      $(EXEDIR)/PaulBunyan \
      $(EXEDIR)/LidarExporter \
      $(EXEDIR)/LidarGridder \
      $(EXEDIR)/LidarViewer \
      $(EXEDIR)/PointSetSimilarity \
      $(EXEDIR)/PrintPrimitiveFile
.PHONY: all
all: $(ALL)

########################################################################
# Specify other actions to be performed on a `make clean'
########################################################################

.PHONY: extraclean
extraclean:

.PHONY: extrasqueakyclean
extrasqueakyclean:

# Include basic makefile
include $(VRUI_MAKEDIR)/BasicMakefile

########################################################################
# Specify build rules for executables
########################################################################

CALCLASRANGE_SOURCES = CalcLasRange.cpp

$(EXEDIR)/CalcLasRange: PACKAGES += MYCOMM
$(EXEDIR)/CalcLasRange: $(OBJDIR)/CalcLasRange.o
.PHONY: CalcLasRange
CalcLasRange: $(EXEDIR)/CalcLasRange

LIDARPREPROCESSOR_SOURCES = SplitPoints.cpp \
                            TempOctree.cpp \
                            PointAccumulator.cpp \
                            LidarProcessOctree.cpp \
                            LidarOctreeCreator.cpp \
                            ReadPlyFile.cpp \
                            LidarPreprocessor.cpp

$(OBJDIR)/LidarPreprocessor.o: CFLAGS += -DLIDARVIEWER_CONFIGFILENAME='"$(ETCINSTALLDIR)/LidarViewer.cfg"'

$(EXEDIR)/LidarPreprocessor: PACKAGES += MYCOMM
$(EXEDIR)/LidarPreprocessor: $(LIDARPREPROCESSOR_SOURCES:%.cpp=$(OBJDIR)/%.o)
.PHONY: LidarPreprocessor
LidarPreprocessor: $(EXEDIR)/LidarPreprocessor

LIDARSPOTREMOVER_SOURCES = LidarProcessOctree.cpp \
                           LidarSpotRemover.cpp
$(EXEDIR)/LidarSpotRemover: $(LIDARSPOTREMOVER_SOURCES:%.cpp=$(OBJDIR)/%.o)
.PHONY: LidarSpotRemover
LidarSpotRemover: $(EXEDIR)/LidarSpotRemover

LIDARSUBTRACTOR_SOURCES = LidarProcessOctree.cpp \
                          SplitPoints.cpp \
                          TempOctree.cpp \
                          LidarOctreeCreator.cpp \
                          PointAccumulator.cpp \
                          SubtractorHelper.cpp \
                          LidarSubtractor.cpp

$(OBJDIR)/LidarSubtractor.o: CFLAGS += -DLIDARVIEWER_CONFIGFILENAME='"$(ETCINSTALLDIR)/LidarViewer.cfg"'

$(EXEDIR)/LidarSubtractor: $(LIDARSUBTRACTOR_SOURCES:%.cpp=$(OBJDIR)/%.o)
.PHONY: LidarSubtractor
LidarSubtractor: $(EXEDIR)/LidarSubtractor

LIDARILLUMINATOR_SOURCES = LidarProcessOctree.cpp \
                           NormalCalculator.cpp \
                           LidarIlluminator.cpp

$(EXEDIR)/LidarIlluminator: $(LIDARILLUMINATOR_SOURCES:%.cpp=$(OBJDIR)/%.o)
.PHONY: LidarIlluminator
LidarIlluminator: $(EXEDIR)/LidarIlluminator

LIDARCOLORMAPPER_SOURCES = LidarProcessOctree.cpp \
                           LidarColorMapper.cpp

$(EXEDIR)/LidarColorMapper: PACKAGES += MYIMAGES
$(EXEDIR)/LidarColorMapper: $(LIDARCOLORMAPPER_SOURCES:%.cpp=$(OBJDIR)/%.o)
.PHONY: LidarColorMapper
LidarColorMapper: $(EXEDIR)/LidarColorMapper

PAULBUNYAN_SOURCES = LidarProcessOctree.cpp \
                     PaulBunyan.cpp

$(EXEDIR)/PaulBunyan: $(PAULBUNYAN_SOURCES:%.cpp=$(OBJDIR)/%.o)
.PHONY: PaulBunyan
PaulBunyan: $(EXEDIR)/PaulBunyan

LIDAREXPORTER_SOURCES = LidarProcessOctree.cpp \
                        LidarExporter.cpp

$(EXEDIR)/LidarExporter: $(LIDAREXPORTER_SOURCES:%.cpp=$(OBJDIR)/%.o)
.PHONY: LidarExporter
LidarExporter: $(EXEDIR)/LidarExporter

LIDARGRIDDER_SOURCES = LidarProcessOctree.cpp \
                       LidarGridder.cpp

$(EXEDIR)/LidarGridder: PACKAGES += MYVRUI
$(EXEDIR)/LidarGridder: $(LIDARGRIDDER_SOURCES:%.cpp=$(OBJDIR)/%.o)
.PHONY: LidarGridder
LidarGridder: $(EXEDIR)/LidarGridder

LIDARVIEWER_SOURCES = LidarOctree.cpp \
                      PointBasedLightingShader.cpp \
                      LidarTool.cpp \
                      PlanePrimitive.cpp \
                      BruntonPrimitive.cpp \
                      LinePrimitive.cpp \
                      PointPrimitive.cpp \
                      SpherePrimitive.cpp \
                      CylinderPrimitive.cpp \
                      ProfileExtractor.cpp \
                      ProfileTool.cpp \
                      SceneGraph.cpp \
                      LidarProcessOctree.cpp \
                      LoadPointSet.cpp \
                      LidarViewer.cpp

$(OBJDIR)/LidarViewer.o: CFLAGS += -DLIDARVIEWER_CONFIGFILENAME='"$(ETCINSTALLDIR)/LidarViewer.cfg"'

$(EXEDIR)/LidarViewer: PACKAGES += MYVRUI
# $(EXEDIR)/LidarViewer: CFLAGS += -DVISUALIZE_WATER
$(EXEDIR)/LidarViewer: $(LIDARVIEWER_SOURCES:%.cpp=$(OBJDIR)/%.o)
.PHONY: LidarViewer
LidarViewer: $(EXEDIR)/LidarViewer

$(EXEDIR)/PointSetSimilarity: $(OBJDIR)/PointSetSimilarity.o
.PHONY: PointSetSimilarity
PointSetSimilarity: $(EXEDIR)/PointSetSimilarity

$(EXEDIR)/PrintPrimitiveFile: $(OBJDIR)/PrintPrimitiveFile.o
.PHONY: PrintPrimitiveFile
PrintPrimitiveFile: $(EXEDIR)/PrintPrimitiveFile

install: $(ALL)
	@echo Installing LiDAR Viewer in $(INSTALLDIR)...
	@install -d $(INSTALLDIR)
	@install -d $(EXECUTABLEINSTALLDIR)
	@install $(ALL) $(EXECUTABLEINSTALLDIR)
	@install -d $(ETCINSTALLDIR)
	@install -m u=rw,go=r $(CONFIGDIR)/LidarViewer.cfg $(ETCINSTALLDIR)
