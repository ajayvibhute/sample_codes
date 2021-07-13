################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/autodbnn.cpp \
../src/mpi_autodbnn.cpp \
../src/mpi_autodbnn_12_March_2013.cpp \
../src/mpi_autodbnn_15_feb_2013.cpp \
../src/mpi_autodbnn_18_feb_2013.cpp \
../src/mpi_autodbnn_20_feb_2013.cpp \
../src/mpi_autodbnn_25_feb_2013.cpp \
../src/mpi_autodbnn_4_march_2013.cpp 

OBJS += \
./src/autodbnn.o \
./src/mpi_autodbnn.o \
./src/mpi_autodbnn_12_March_2013.o \
./src/mpi_autodbnn_15_feb_2013.o \
./src/mpi_autodbnn_18_feb_2013.o \
./src/mpi_autodbnn_20_feb_2013.o \
./src/mpi_autodbnn_25_feb_2013.o \
./src/mpi_autodbnn_4_march_2013.o 

CPP_DEPS += \
./src/autodbnn.d \
./src/mpi_autodbnn.d \
./src/mpi_autodbnn_12_March_2013.d \
./src/mpi_autodbnn_15_feb_2013.d \
./src/mpi_autodbnn_18_feb_2013.d \
./src/mpi_autodbnn_20_feb_2013.d \
./src/mpi_autodbnn_25_feb_2013.d \
./src/mpi_autodbnn_4_march_2013.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/opt/openmpi/1.4.2/include/ -includempi.h -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


