all : 

PROJ_NAME := ring

DEV_DIR := .
BUILD_DIR := $(DEV_DIR)/build
INSTALL_DIR := $(DEV_DIR)/install
GTEST_DIR := $(DEV_DIR)/third_party/gtest

SHELL := /bin/sh
CXX := g++-4.7
LD := g++-4.7
CP := cp -r
RSYNC := rsync -iCau --exclude='\.*' --delete
AR := ar
RM := rm -f
ZIP := zip
MKDIR := mkdir -p

INCLUDE_FLAGS := -I$(INSTALL_DIR)/include \
                 -I$(DEV_DIR)/include \
								 -I$(GTEST_DIR)/include \
								 -I$(GTEST_DIR) \
                 -I/usr/local/Cellar/include/glew/1.10.0/include \
                 -I/usr/local/Cellar/include/GLFW 

FLAGS := -m64 -Wall -Wextra -Wshadow -Werror -pedantic
CXXFLAGS := -std=c++11 -Weffc++ $(FLAGS)
LDFLAGS := -L$(INSTALL_DIR)/lib -lm -lpthread /usr/lib/libgmp.10.dylib -lflint -lssl -lcrypto

DEBUG_FLAGS := -g -O2 -D _DEBUG
RELEASE_FLAGS := -O2 -D NDEBUG

ARTIFACTS := 

include make/lib.mk
include make/demo.mk
#include make/another_demo.mk
include make/test.mk

BUILD_DIRS := $(sort $(BUILD_DIR_PROJ) \
		$(BUILD_DIR_PROJ_DEMO) \
		$(BUILD_DIR_PROJ_BINARY) \
		$(BUILD_DIR_PROJ_TEST)) 

TARGETS := $(sort $(TARGET_PROJ) \
                  $(TARGET_PROJ_DEMO) \
                  $(TARGET_PROJ_BINARY) \
                  $(TARGET_PROJ_TEST)) 

all : $(TARGETS)
lib : $(TARGET_PROJ)
demo : $(TARGET_PROJ_DEMO)
another-demo: $(TARGET_PROJ_BINARY)
test : $(TARGET_PROJ_TEST)

COMPILE_CXX = $(CXX) -c -o $@ -fopenmp $< $(INCLUDE_FLAGS) $(CXXFLAGS) $(DEBUG_FLAGS) 

release : COMPILE_CXX = $(CXX) -c -o $@ -fopenmp $< $(INCLUDE_FLAGS) $(CXXFLAGS) $(RELEASE_FLAGS)
release : all

%-test.x : | $(INSTALL_DIR)/test $(INSTALL_DIR)/test/resources
	$(LD) -o $@ -fopenmp $^ $(LDFLAGS)
	$@

%.a : | $(INSTALL_DIR)/lib $(INSTALL_DIR)/include
	$(AR) rcs $@ $^

%.x : | $(INSTALL_DIR)/bin $(INSTALL_DIR)/bin/resources
	$(LD) -o $@ -fopenmp $^ $(LDFLAGS)

$(TARGET_PROJ) : $(OBJECTS_PROJ)
$(TARGET_PROJ_DEMO) : $(OBJECTS_PROJ_DEMO)
$(TARGET_PROJ_BINARY) : $(OBJECTS_PROJ_BINARY)
$(TARGET_PROJ_TEST) : $(OBJECTS_PROJ_TEST)

$(OBJECTS_PROJ) $(patsubst %.o,%.d, $(OBJECTS_PROJ)) : | $(BUILD_DIR_PROJ)
$(OBJECTS_PROJ_DEMO) $(patsubst %.o,%.d, $(OBJECTS_PROJ_DEMO)) : | $(BUILD_DIR_PROJ_DEMO)
$(OBJECTS_PROJ_BINARY) $(patsubst %.o,%.d, $(OBJECTS_PROJ_BINARY)) : | $(BUILD_DIR_PROJ_BINARY)
$(OBJECTS_PROJ_TEST) $(patsubst %.o,%.d, $(OBJECTS_PROJ_TEST)) : | $(BUILD_DIR_PROJ_TEST)

OBJECTS := $(sort $(OBJECTS_PROJ) \
		$(OBJECTS_PROJ_DEMO) \
		$(OBJECTS_PROJ_BINARY) \
		$(OBJECTS_PROJ_TEST))

$(OBJECTS) : %.o : %.d


$(addsuffix /%.d, $(BUILD_DIRS)) : %.cpp
	$(SHELL) -ec "$(CXX) -std=c++11 $(INCLUDE_FLAGS) -I$(DEV_DIR)/include -MM $< \
	| sed 's|$(notdir $*)\.o[ ]*:|$*\.o $@ :|g' > $@; \
	[ -s $@ ] || $(RM) $@"

$(addsuffix /%.d, $(BUILD_DIRS)) : %.cc
	$(SHELL) -ec "$(CXX) -std=c++11 $(INCLUDE_FLAGS) -I$(DEV_DIR)/include -MM $< \
	| sed 's|$(notdir $*)\.o[ ]*:|$*\.o $@ :|g' > $@; \
	[ -s $@ ] || $(RM) $@"


$(BUILD_DIR_PROJ)/%.o : %.cpp ; $(COMPILE_CXX)
$(BUILD_DIR_PROJ_DEMO)/%.o : %.cpp ; $(COMPILE_CXX)
$(BUILD_DIR_PROJ_BINARY)/%.o : %.cpp ; $(COMPILE_CXX)
$(BUILD_DIR_PROJ_TEST)/%.o : %.cpp ; $(COMPILE_CXX)
$(BUILD_DIR_PROJ_TEST)/%.o : %.cc ; $(COMPILE_CXX)

clean:
	$(RM) -r $(TARGETS) $(BUILD_DIRS) $(ARTIFACTS)

-include $(patsubst %.o,%.d,$(OBJECTS))

$(BUILD_DIRS) $(INSTALL_DIR)/bin \
		$(INSTALL_DIR)/bin/resources \
		$(INSTALL_DIR)/include \
		$(INSTALL_DIR)/lib \
		$(INSTALL_DIR)/test/resources \
		$(INSTALL_DIR)/test :
	$(MKDIR) $@

.PHONY : all clean release lib demo test another-demo

.PRECIOUS : %-test.x

