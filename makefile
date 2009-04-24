########################################################################
# Makefile for LiDAR Viewer, a visualization and analysis application
# for large 3D point cloud data.
# Copyright (c) 2004-2009 Oliver Kreylos
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

# Root directory of the Vrui software installation. This must match the
# same setting in Vrui's makefile. By default the directories match; if
# the installation directory was adjusted during Vrui's installation, it
# must be adjusted here as well.
VRUIDIR = $(HOME)/Vrui-1.0

# Base installation directory for LiDAR Viewer and its configuration
# file.If this is set to the default of $(PWD), LiDAR Viewer does not
# have to be installed to be run. LiDAR Viewer's executables, and
# configuration file will be installed in the bin and etc directories
# underneath the given base directory, respectively.
INSTALLDIR = $(shell pwd)

# Set up additional flags for the C++ compiler:
CFLAGS = 

# Set up destination directories for compilation products:
OBJDIRBASE = o
BINDIRBASE = bin

# Create debug or fully optimized versions of the software:
ifdef DEBUG
  # Include the debug version of the Vrui application makefile fragment:
  include $(VRUIDIR)/etc/Vrui.debug.makeinclude
  # Enable debugging and disable optimization:
  CFLAGS += -g3 -O0
  # Set destination directories for created objects:
  OBJDIR = $(OBJDIRBASE)/debug
  BINDIR = $(BINDIRBASE)/debug
else
  # Include the release version of the Vrui application makefile fragment:
  include $(VRUIDIR)/etc/Vrui.makeinclude
  # Disable debugging and enable optimization:
  CFLAGS += -g0 -O3 -DNDEBUG
  # Set destination directories for created objects:
  OBJDIR = $(OBJDIRBASE)
  BINDIR = $(BINDIRBASE)
endif

# Pattern rule to compile C++ sources:
$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(OBJDIR)/$(*D)
	@echo Compiling $<...
	@g++ -c -o $@ $(VRUI_CFLAGS) $(CFLAGS) $<

# Rule to build all LiDAR Viewer components:
ALL = $(BINDIR)/CalcLasRange \
      $(BINDIR)/LidarPreprocessor \
      $(BINDIR)/LidarSubtractor \
      $(BINDIR)/LidarIlluminator \
      $(BINDIR)/LidarViewer \
      $(BINDIR)/PointSetSimilarity \
      $(BINDIR)/PrintPrimitiveFile
.PHONY: all
all: $(ALL)

# Rule to remove build results:
clean:
	-rm -f $(OBJDIR)/*.o
	-rm -f $(ALL)
	-rmdir $(BINDIR)

# Rule to clean the source directory for packaging:
distclean:
	-rm -rf $(OBJDIRBASE)
	-rm -rf $(BINDIRBASE)

CALCLASRANGE_SOURCES = CalcLasRange.cpp

$(BINDIR)/CalcLasRange: $(CALCLASRANGE_SOURCES:%.cpp=$(OBJDIR)/%.o)
	@mkdir -p $(BINDIR)
	@echo Linking $@...
	@g++ -o $@ $^ $(VRUI_LINKFLAGS)
.PHONY: CalcLasRange
CalcLasRange: $(BINDIR)/CalcLasRange

LIDARPREPROCESSOR_SOURCES = SplitPoints.cpp \
                            TempOctree.cpp \
                            PointAccumulator.cpp \
                            LidarProcessOctree.cpp \
                            LidarOctreeCreator.cpp \
                            LidarPreprocessor.cpp

$(OBJDIR)/LidarPreprocessor.o: CFLAGS += -DLIDARVIEWER_CONFIGFILENAME='"$(INSTALLDIR)/etc/LidarViewer.cfg"'

$(BINDIR)/LidarPreprocessor: $(LIDARPREPROCESSOR_SOURCES:%.cpp=$(OBJDIR)/%.o)
	@mkdir -p $(BINDIR)
	@echo Linking $@...
	@g++ -o $@ $^ $(VRUI_LINKFLAGS)
.PHONY: LidarPreprocessor
LidarPreprocessor: $(BINDIR)/LidarPreprocessor

LIDARSUBTRACTOR_SOURCES = LidarProcessOctree.cpp \
                          SplitPoints.cpp \
                          TempOctree.cpp \
                          LidarOctreeCreator.cpp \
                          PointAccumulator.cpp \
                          LidarSubtractor.cpp

$(BINDIR)/LidarSubtractor: $(LIDARSUBTRACTOR_SOURCES:%.cpp=$(OBJDIR)/%.o)
	@mkdir -p $(BINDIR)
	@echo Linking $@...
	@g++ -o $@ $^ $(VRUI_LINKFLAGS)
.PHONY: LidarSubtractor
LidarSubtractor: $(BINDIR)/LidarSubtractor

LIDARILLUMINATOR_SOURCES = LidarProcessOctree.cpp \
                           NormalCalculator.cpp \
                           LidarIlluminator.cpp

$(BINDIR)/LidarIlluminator: $(LIDARILLUMINATOR_SOURCES:%.cpp=$(OBJDIR)/%.o)
	@mkdir -p $(BINDIR)
	@echo Linking $@...
	@g++ -o $@ $^ $(VRUI_LINKFLAGS)
.PHONY: LidarIlluminator
LidarIlluminator: $(BINDIR)/LidarIlluminator

LIDARVIEWER_SOURCES = LidarOctree.cpp \
                      LidarTool.cpp \
                      PlanePrimitive.cpp \
                      LinePrimitive.cpp \
                      PointPrimitive.cpp \
                      SpherePrimitive.cpp \
                      CylinderPrimitive.cpp \
                      PointBasedLightingShader.cpp \
                      LidarViewer.cpp

$(OBJDIR)/LidarViewer.o: CFLAGS += -DLIDARVIEWER_CONFIGFILENAME='"$(INSTALLDIR)/etc/LidarViewer.cfg"'

$(BINDIR)/LidarViewer: $(LIDARVIEWER_SOURCES:%.cpp=$(OBJDIR)/%.o)
	@mkdir -p $(BINDIR)
	@echo Linking $@...
	@g++ -o $@ $^ $(VRUI_LINKFLAGS)
.PHONY: LidarViewer
LidarViewer: $(BINDIR)/LidarViewer

$(BINDIR)/PointSetSimilarity: $(OBJDIR)/PointSetSimilarity.o
	@mkdir -p $(BINDIR)
	@echo Linking $@...
	@g++ -o $@ $^ $(VRUI_LINKFLAGS)
.PHONY: PointSetSimilarity
PointSetSimilarity: $(BINDIR)/PointSetSimilarity

$(BINDIR)/PrintPrimitiveFile: $(OBJDIR)/PrintPrimitiveFile.o
	@mkdir -p $(BINDIR)
	@echo Linking $@...
	@g++ -o $@ $^ $(VRUI_LINKFLAGS)
.PHONY: PrintPrimitiveFile
PrintPrimitiveFile: $(BINDIR)/PrintPrimitiveFile

install: $(ALL)
	@echo Installing LiDAR Viewer in $(INSTALLDIR)...
	@install -d $(INSTALLDIR)
	@install -d $(INSTALLDIR)/bin
	@install $(BINDIR)/CalcLasRange $(INSTALLDIR)/bin
	@install $(BINDIR)/LidarPreprocessor $(INSTALLDIR)/bin
	@install $(BINDIR)/LidarSubtractor $(INSTALLDIR)/bin
	@install $(BINDIR)/LidarIlluminator $(INSTALLDIR)/bin
	@install $(BINDIR)/LidarViewer $(INSTALLDIR)/bin
	@install $(BINDIR)/PointSetSimilarity $(INSTALLDIR)/bin
	@install $(BINDIR)/PrintPrimitiveFile $(INSTALLDIR)/bin
	@install -d $(INSTALLDIR)/etc
	@install etc/LidarViewer.cfg $(INSTALLDIR)/etc
