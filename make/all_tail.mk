## Here we set up the directories that all of the 
## object files will be ending up in.

BUILD_DIRS := $(sort $(BUILD_DIR_PROJ) $(BUILD_DIR_PROJ_DEMO) $(BUILD_DIR_PROJ_BINARY) $(BUILD_DIR_PROJ_TEST))

TARGETS := $(sort $(TARGET_PROJ) \
                  $(TARGET_PROJ_DEMO) \
                  $(TARGET_PROJ_BINARY) \
                  $(TARGET_PROJ_TEST)) 

## UPPER LEVEL (PHONY) TARGETS

all : $(TARGETS)
nbody : $(TARGET_PROJ)
nbody-demo : $(TARGET_PROJ_DEMO)
nbody-binary : $(TARGET_PROJ_BINARY)
nbody-test : $(TARGET_PROJ_TEST)

## Set compile flags for debug mode by default
COMPILE_CXX = $(CXX) -c -o $@ $< $(INCLUDE_FLAGS) $(CXXFLAGS) $(DEBUG_FLAGS) 

## Running 'make release' will set compile flags to release mode
release : COMPILE_CXX = $(CXX) -c -o $@ $< $(INCLUDE_FLAGS) $(CXXFLAGS) $(RELEASE_FLAGS)
release : all

## Rule for making gtest executables:
## First we link %-test.x, then we execute it to run the unit tests.
%-test.x : | $(INSTALL_DIR)/test $(INSTALL_DIR)/test/resources
	$(LD) -o $@ $^ $(LDFLAGS)
	$@

## Rule for making any static libraries. If you want to find out
## more, see http://linux.die.net/man/1/ar
%.a : | $(INSTALL_DIR)/lib $(INSTALL_DIR)/include
	$(AR) rcs $@ $^

## Rule for making any regular executables - link the object files.
%.x : | $(INSTALL_DIR)/bin $(INSTALL_DIR)/bin/resources
	$(LD) -o $@ $^ $(LDFLAGS)

## The targets depend on each of their objects. If any object
## files are newer, we'll rebuild the target.
$(TARGET_PROJ) : $(OBJECTS_PROJ)
$(TARGET_PROJ_DEMO) : $(OBJECTS_PROJ_DEMO)
$(TARGET_PROJ_BINARY) : $(OBJECTS_PROJ_BINARY)
$(TARGET_PROJ_TEST) : $(OBJECTS_PROJ_TEST)

## We need BUILD_DIR_PROJ around before we try to generate the object
## files or dependency files. See the link:
## http://www.gnu.org/software/make/manual/make.html#Prerequisite-Types
$(OBJECTS_PROJ) $(patsubst %.o,%.d, $(OBJECTS_PROJ)) : | $(BUILD_DIR_PROJ)
$(OBJECTS_PROJ_DEMO) $(patsubst %.o,%.d, $(OBJECTS_PROJ_DEMO)) : | $(BUILD_DIR_PROJ_DEMO)
$(OBJECTS_PROJ_BINARY) $(patsubst %.o,%.d, $(OBJECTS_PROJ_BINARY)) : | $(BUILD_DIR_PROJ_BINARY)
$(OBJECTS_PROJ_TEST) $(patsubst %.o,%.d, $(OBJECTS_PROJ_TEST)) : | $(BUILD_DIR_PROJ_TEST)
## patsubst = "pattern substitution". For this and addsuffix below, see
## http://www.gnu.org/software/make/manual/make.html#Text-Functions

OBJECTS := $(sort $(OBJECTS_PROJ) $(OBJECTS_PROJ_DEMO) $(OBJECTS_PROJ_BINARY) $(OBJECTS_PROJ_TEST))
$(OBJECTS) : %.o : %.d

## For each required dependency file in BUILD_DIRS, we use sed and
## g++ -MM to generate a list of the included headers (dependencies)
$(addsuffix /%.d, $(BUILD_DIRS)) : %.cpp
	$(SHELL) -ec "$(CXX) -std=c++11 $(INCLUDE_FLAGS) -I$(DEV_DIR)/include -MM $< \
	| sed 's|$(notdir $*)\.o[ ]*:|$*\.o $@ :|g' > $@; \
	[ -s $@ ] || $(RM) $@"

## You should name all of your files so they end in .cpp, but gtest
## gtest will have some files that end in .cc we need to take care of.
$(addsuffix /%.d, $(BUILD_DIRS)) : %.cc
	$(SHELL) -ec "$(CXX) -std=c++11 $(INCLUDE_FLAGS) -I$(DEV_DIR)/include -MM $< \
	| sed 's|$(notdir $*)\.o[ ]*:|$*\.o $@ :|g' > $@; \
	[ -s $@ ] || $(RM) $@"


## Build all of the objects from our source files and do it
## in our build/projectname directory.
$(BUILD_DIR_PROJ)/%.o : %.cpp ; $(COMPILE_CXX)
$(BUILD_DIR_PROJ_DEMO)/%.o : %.cpp ; $(COMPILE_CXX)
$(BUILD_DIR_PROJ_BINARY)/%.o : %.cpp ; $(COMPILE_CXX)
$(BUILD_DIR_PROJ_TEST)/%.o : %.cpp ; $(COMPILE_CXX)
$(BUILD_DIR_PROJ_TEST)/%.o : %.cc ; $(COMPILE_CXX)

## Command to remove all of the contents of the build/ directory
## related to our project and all of the generated executables
## and copied header files found in install/
clean:
	$(RM) -r $(TARGETS) $(BUILD_DIRS) $(ARTIFACTS)

## Include the '.d' (dependencies) files; the leading '-' dash
## indicates that we don't want make to warn us if the files
## don't already exist.
-include $(patsubst %.o,%.d,$(OBJECTS))

## Make all of the build and install directories if they do not already exist
$(BUILD_DIRS) $(INSTALL_DIR)/bin \
						  $(INSTALL_DIR)/bin/resources \
						  $(INSTALL_DIR)/include \
							$(INSTALL_DIR)/lib \
							$(INSTALL_DIR)/test/resources \
							$(INSTALL_DIR)/test :
	$(MKDIR) $@

# These names don't represent real targets
.PHONY : all clean release $(PROJ_NAME) $(PROJ_NAME)-demo $(PROJ_NAME)-test

# Don't delete our unit test executable even if make is killed or 
# interrupted while it's being run.
.PRECIOUS : %-test.x
