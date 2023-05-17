# ================= Libraries and Includes ====================================

DATASTRUCT_INCLUDE_DIR = ./DataStruct/include
INOUT_INCLUDE_DIR = ./InOut/include
METRIC_INCLUDE_DIR = ./Metrics/include
CURV_INCLUDE_DIR = ./Curvature/include
OT1_INCLUDE_DIR = ./OT1/include
PH0_INCLUDE_DIR = ./PH0/include
CLUSTERING_INCLUDE_DIR = ./Clustering/include
NANOGUI_INCLUDE_DIR = ./nanogui/include

BLAS_LIBS = -lblas -llapack
THREAD_LIBS = -lpthread

# ================= Project Directories ====================================

INC_DIR = ./project/include
SRC_DIR = ./project/src
OBJ_DIR = ./project/src
BIN_DIR = ./bin

# ================= Project Name ===========================================

EXT=
NAME1=Comp2DShapes
NAMEOBJ1=$(OBJ_DIR)/$(NAME1).o
NAMEBIN1=$(BIN_DIR)/$(NAME1)$(EXT)
NAME2=2DShape
NAMEOBJ2=$(OBJ_DIR)/$(NAME2).o
NAMEBIN2=$(BIN_DIR)/$(NAME2)$(EXT)
NAME3=Clustering
NAMEOBJ3=$(OBJ_DIR)/$(NAME3).o
NAMEBIN3=$(BIN_DIR)/$(NAME3)$(EXT)
NAME4=GUI
NAMEOBJ4=$(OBJ_DIR)/$(NAME4).o
NAMEBIN4=$(BIN_DIR)/$(NAME4)$(EXT)
NAME5=test
NAMEOBJ5=$(OBJ_DIR)/$(NAME5).o
NAMEBIN5=$(BIN_DIR)/$(NAME5)$(EXT)

# ================= Compilers and Flags ====================================

CC 		:= gcc
CFLAGS		:= -c -O3 -ansi -Wall -Werror -pedantic
CPP		:= g++
CPPFLAGS	:= -c -std=c++17 -O3 -Wdeprecated-declarations -Wuninitialized -ansi -Wall -std=c++0x

INCLUDE_DIRS = -I$(INC_DIR) -I$(DATASTRUCT_INCLUDE_DIR) -I$(INOUT_INCLUDE_DIR)\
		-I$(METRIC_INCLUDE_DIR) -I$(CURV_INCLUDE_DIR) -I$(OT1_INCLUDE_DIR) -I$(PH0_INCLUDE_DIR) -I$(CLUSTERING_INCLUDE_DIR) -I$(NANOGUI_INCLUDE_DIR)

LIB_DIRS =

LIBS = $(BLAS_LIBS) $(THREAD_LIBS) -lstdc++

LOAD_LIB_PATH =

LD_FLAGS = -O3

# ================= Pattern rules ==========================================

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDE_DIRS) $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDE_DIRS) $< -o $@

$(F_OBJ_DIR)/%.o: $(F_SRC_DIR)/%.f
	$(FC) $(FFLAGS) $< -o $@

# ================= Compile source code ====================================

OBJECTS1 = \
$(NAMEOBJ1)

OBJECTS2 = \
$(NAMEOBJ2)

OBJECTS3 = \
$(NAMEOBJ3)

OBJECTS4 = \
$(NAMEOBJ4)

OBJECTS5 = \
$(NAMEOBJ5)



# ================= Generate Executable ====================================

$(NAMEBIN1) : $(OBJECTS1)
	$(CPP) -o $(NAMEBIN1) $(LD_FLAGS) $(OBJECTS1) $(LIB_DIRS) $(LIBS) $(LOAD_LIB_PATH)

$(NAMEBIN2) : $(OBJECTS2)
	$(CPP) -o $(NAMEBIN2) $(LD_FLAGS) $(OBJECTS2) $(LIB_DIRS) $(LIBS) $(LOAD_LIB_PATH)

$(NAMEBIN3) : $(OBJECTS3)
	$(CPP) -o $(NAMEBIN3) $(LD_FLAGS) $(OBJECTS3) $(LIB_DIRS) $(LIBS) $(LOAD_LIB_PATH)

$(NAMEBIN4) : $(OBJECTS4)
	$(CPP) -o $(NAMEBIN4) $(LD_FLAGS) $(OBJECTS4) $(LIB_DIRS) $(LIBS) $(LOAD_LIB_PATH)

$(NAMEBIN5) : $(OBJECTS5)
	$(CPP) -o $(NAMEBIN5) $(LD_FLAGS) $(OBJECTS5) $(LIB_DIRS) $(LIBS) $(LOAD_LIB_PATH)


all: $(OBJECTS1) $(OBJECTS2) $(OBJECTS3) #$(OBJECTS4)
	$(CPP) -o $(NAMEBIN1) $(LD_FLAGS) $(OBJECTS1) $(LIB_DIRS) $(LIBS) $(LOAD_LIB_PATH)
	$(CPP) -o $(NAMEBIN2) $(LD_FLAGS) $(OBJECTS2) $(LIB_DIRS) $(LIBS) $(LOAD_LIB_PATH)
	$(CPP) -o $(NAMEBIN3) $(LD_FLAGS) $(OBJECTS3) $(LIB_DIRS) $(LIBS) $(LOAD_LIB_PATH)
#$(CPP) -o $(NAMEBIN4) $(LD_FLAGS) $(OBJECTS4) $(LIB_DIRS) $(LIBS) $(LOAD_LIB_PATH)

clustering: $(OBJECTS3)
	$(CPP) -o $(NAMEBIN3) $(LD_FLAGS) $(OBJECTS3) $(LIB_DIRS) $(LIBS) $(LOAD_LIB_PATH)

GUI: $(OBJECTS4)
	$(CPP) -o $(NAMEBIN4) $(LD_FLAGS) $(OBJECTS4) $(LIB_DIRS) $(LIBS) $(LOAD_LIB_PATH)

test: $(OBJECTS5)
	$(CPP) -o $(NAMEBIN5) $(LD_FLAGS) $(OBJECTS5) $(LIB_DIRS) $(LIBS) $(LOAD_LIB_PATH)

clean:
	touch $(OBJ_DIR)/junk.o; rm -f $(OBJ_DIR)/*.o $(NAMEBIN1) $(NAMEBIN2) $(NAMEBIN3)

$(OBJECTS1):
$(OBJECTS2):
$(OBJECTS3):
$(OBJECTS4):
$(OBJECTS5):

