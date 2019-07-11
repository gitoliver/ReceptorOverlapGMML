TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    src/main.cpp \
    src/io.cpp \
    src/residue_name_conversion_map.cpp \
    src/bead_residues.cpp

HEADERS += \
    includes/io.h \
    includes/residue_name_conversion_map.h \
    includes/bead_residues.h
