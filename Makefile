CXX = g++
CXXFLAGS = -std=c++20 -O3 -march=native -ftree-vectorize -ffast-math -Wall #-g

SRC = src
BIN = bin
LIB = lib
BUILD = bin

CXXFLAGS += -Ilib/num_array/include

TARGET = $(BIN)/libsppe.a
OBJECTS = $(addprefix $(BUILD)/, sppe.o particle.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(AR) -crs $(TARGET) $(OBJECTS)

#$(BUILD)/%.o: $(SRC)/%.cc $(SRC)/%.h
#	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILD)/sppe.o: $(addprefix $(SRC)/, sppe.cc sppe.h spatial_map.h types.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILD)/particle.o: $(addprefix $(SRC)/, particle.cc particle.h types.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	$(RM) $(OBJECTS)
