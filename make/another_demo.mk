## Look in src/demo for the demo main files
vpath %.cpp $(DEV_DIR)/src/binary

## The result of our nbody-demo compile: an executable
## named list-demo.x which will run our program
TARGET_PROJ_BINARY := $(INSTALL_DIR)/bin/$(PROJ_NAME)-binary.x
BUILD_DIR_PROJ_BINARY := $(BUILD_DIR)/$(PROJ_NAME)-binary

## nbody-demo.x depends on libnbody.a having been created first
$(TARGET_PROJ_BINARY) : | $(INSTALL_DIR)/lib/lib$(PROJ_NAME).a \

## Link our nbody-demo.x executable with libnbody.a
## (here is where you add any extra third-party library linkings needed)
$(TARGET_PROJ_BINARY) : LDFLAGS += $(INSTALL_DIR)/lib/lib$(PROJ_NAME).a \
                                 -framework Cocoa \
                                 -framework OpenGL \
                                 -framework IOKit \
                                 -framework CoreVideo \
                                 -lglew \
                                 /usr/local/Cellar/glew/1.10.0/lib/libGLEW.a \
                                 /usr/local/Cellar/glfw3/3.0.4/lib/libglfw3.3.dylib \

## Add any more files to this list
OBJECTS_PROJ_BINARY := \
	$(BUILD_DIR_PROJ_BINARY)/binary.o \
  $(BUILD_DIR_PROJ_DEMO)/GlutWrapper.o \
  $(BUILD_DIR_PROJ_DEMO)/Shaders.o \
	#$(BUILD_DIR_PROJ_BINARY)/otherFiles.o \
	#$(BUILD_DIR_PROJ_BINARY)/moreOtherFiles.o \
## Add any more files you like!

## Add the public include files to the list of artifacts
## to be cleaned up on a 'make clean'
ARTIFACTS += $(INSTALL_DIR)/include/$(PROJ_NAME)-binary
ARTIFACTS += $(ISNTALL_DIR)/bin/resources/$(PROJ_NAME)-binary


