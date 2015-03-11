############################################################################
# CMakeLists.txt file for building ROOT io/hdfs package
############################################################################

ROOT_USE_PACKAGE(io/io)
ROOT_USE_PACKAGE(net/net)
include_directories(${HDFS_INCLUDE_DIRS})
add_definitions(-pthread)

ROOT_GENERATE_DICTIONARY(G__HDFS *.h MODULE HDFS LINKDEF LinkDef.h)

ROOT_LINKER_LIBRARY(HDFS *.cxx G__HDFS.cxx LIBRARIES Core ${HDFS_LIBRARIES} DEPENDENCIES Net RIO)
ROOT_INSTALL_HEADERS()
