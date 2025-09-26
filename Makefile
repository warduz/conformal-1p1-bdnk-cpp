CXX = g++
CXXFLAGS = -O3 -march=native # -fopenmp
# CXXFLAGS = -O0 -g -fstack-protector-all -Wall -Wextra # debug

BUILD_DIR = build
CONFIG = "config.toml"

SRC = main.cpp Config.cpp Evolver.cpp initialize.cpp kt.cpp Metric.cpp
OBJS = $(addprefix $(BUILD_DIR)/, $(SRC:.cpp=.o))

TARGET = $(BUILD_DIR)/solver

all: $(TARGET)

$(TARGET): $(OBJS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)

run: $(TARGET)
	./$(TARGET) $(CONFIG)