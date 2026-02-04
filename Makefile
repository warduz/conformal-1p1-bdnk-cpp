CXX = g++
CXXFLAGS = -O3 -march=native # -fopenmp
# CXXFLAGS = -O0 -g -fstack-protector-all -Wall -Wextra # debug

BUILD_DIR = build
CONFIG = "config/config.toml"

SRC = src/main.cpp src/Config.cpp src/Evolver.cpp src/initialize.cpp src/kt.cpp src/Metric.cpp
OBJS = $(addprefix $(BUILD_DIR)/, $(notdir $(SRC:.cpp=.o)))

TARGET = $(BUILD_DIR)/solver

all: $(TARGET)

$(TARGET): $(OBJS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

$(BUILD_DIR)/%.o: src/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)

run: $(TARGET)
	./$(TARGET) $(CONFIG)