.SUFFIXES :
.SUFFIXES : .cpp .o 

# c ============================================
CC := g++
# CPP_FLAGS := -O2 -m64
CPP_FLAGS := -g -m64
CPP_LDFLAGS := -lm -lnsl -m64 # -lpthread
LINK := $(CC) -fPIC
# end c ============================================

# commands
MAKE = make
AR = ar cr
RM = -rm -rf

# includes
INCLUDES := -I.

# cpp
SRC_CPP=$(wildcard src/*.cpp)
OBJ_CPP=$(SRC_CPP:%.cpp=%.o)

TARGET=CoreFAnaly

$(TARGET):$(OBJ_CPP)
	$(LINK) -o $@ $^ $(CPP_LDFLAGS) $(LDFLAGS)

$(OBJ_CPP):%.o:%.cpp
	$(CC) -c $< -o $@ $(INCLUDES) $(CFLAGS) $(CPP_FLAGS)

clean:
	$(RM) $(TARGET) $(OBJ_CPP) $(LIB)

clean-all:
	$(RM) $(TARGET) $(OBJ_CPP) $(LIB)
	$(RM) bin/*
	$(RM) lib/*
