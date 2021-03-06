vpath %.cpp $(DEV_DIR)/src

TARGET_PROJ := $(INSTALL_DIR)/lib/lib$(PROJ_NAME).a

BUILD_DIR_PROJ := $(BUILD_DIR)/$(PROJ_NAME)

OBJECTS_PROJ := \
	$(BUILD_DIR_PROJ)/GRingArray.o \
	$(BUILD_DIR_PROJ)/IBE.o \
	$(BUILD_DIR_PROJ)/ABE.o \
	$(BUILD_DIR_PROJ)/Trapdoor.o \

$(OBJECTS_PROJ) : | $(INSTALL_DIR)/include/$(PROJ_NAME) \
		$(INSTALL_DIR)/bin/resources/$(PROJ_NAME) \
		$(INSTALL_DIR)/test/resources/$(PROJ_NAME)

$(INSTALL_DIR)/include/$(PROJ_NAME) ::
	$(MKDIR) $@
	$(RSYNC) $(DEV_DIR)/include/$(PROJ_NAME)/ $@/

$(INSTALL_DIR)/bin/resources/$(PROJ_NAME) ::
	$(MKDIR) $@
	$(RSYNC) $(DEV_DIR)/resources/$(PROJ_NAME)/ $@/

$(INSTALL_DIR)/test/resources/$(PROJ_NAME) ::
	$(MKDIR) $@
	$(RSYNC) $(DEV_DIR)/resources/$(PROJ_NAME)/ $@/

ARTIFACTS += $(INSTALL_DIR)/include/$(PROJ_NAME)
ARTIFACTS += $(INSTALL_DIR)/bin/resources/$(PROJ_NAME)
ARTIFACTS += $(INSTALL_DIR)/test/resources/$(PROJ_NAME)

