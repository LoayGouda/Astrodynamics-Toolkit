REM Directory where modules are located
SET MODULES_DIR=modules

REM List of source files
SET SRC_FILES=%MODULES_DIR%\numbers_module.f90 %MODULES_DIR%\vector_module.f90 %MODULES_DIR%\newton_module.f90 ^
%MODULES_DIR%\lambert_module.f90 %MODULES_DIR%\kepler_module.f90 astrodynamics_lib.f90


REM Compiler and flags
SET FC=gfortran
SET FLAGS=-O2

REM Compile all source files
%FC% %FLAGS% -o astrodynamics_lib.exe %SRC_FILES%

REM Execute the compiled program
astrodynamics_lib.exe
