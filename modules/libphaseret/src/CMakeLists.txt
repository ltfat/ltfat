include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

SET(sources
    dgtrealwrapper.c gla.c legla.c pghi.c rtisila.c rtpghi.c spsi.c utils.c
    gsrtisila.c gsrtisilapghi.c)

SET(sources_typeconstant
    dgtrealwrapper_typeconstant.c legla_typeconstant.c pghi_typeconstant.c)

if (USECPP)
    SET_SOURCE_FILES_PROPERTIES( ${sources} ${sources_typeconstant} PROPERTIES LANGUAGE CXX)
endif (USECPP)

add_library(commondouble OBJECT ${sources_typeconstant})
add_library(commonsingle OBJECT ${sources_typeconstant})

add_library(double OBJECT ${sources})
add_library(single OBJECT ${sources})

target_compile_definitions(commondouble PUBLIC LTFAT_DOUBLE)
target_compile_definitions(commonsingle PUBLIC LTFAT_SINGLE)
target_compile_definitions(double PUBLIC LTFAT_DOUBLE)
target_compile_definitions(single PUBLIC LTFAT_SINGLE)

if (MSVC)
    target_compile_definitions(commondouble PUBLIC PHASERET_BUILD_SHARED)
    target_compile_definitions(commonsingle PUBLIC PHASERET_BUILD_SHARED)
    target_compile_definitions(double PUBLIC PHASERET_BUILD_SHARED)
    target_compile_definitions(single PUBLIC PHASERET_BUILD_SHARED)
endif(MSVC)

add_library(phaseret_static STATIC
    $<TARGET_OBJECTS:commondouble> $<TARGET_OBJECTS:double> $<TARGET_OBJECTS:single>)
set_target_properties(phaseret_static PROPERTIES OUTPUT_NAME phaseret)

add_library(phaseretd_static STATIC
    $<TARGET_OBJECTS:commondouble> $<TARGET_OBJECTS:double>)
set_target_properties(phaseretd_static PROPERTIES OUTPUT_NAME phaseretd)

add_library(phaseretf_static STATIC
    $<TARGET_OBJECTS:commonsingle> $<TARGET_OBJECTS:single>)
set_target_properties(phaseretf_static PROPERTIES OUTPUT_NAME phaseretf)

target_link_libraries(phaseret_static)
target_link_libraries(phaseretf_static)
target_link_libraries(phaseretd_static)

add_library(phaseret SHARED
    $<TARGET_OBJECTS:commondouble> $<TARGET_OBJECTS:double> $<TARGET_OBJECTS:single>)


add_library(phaseretd SHARED
    $<TARGET_OBJECTS:commondouble> $<TARGET_OBJECTS:double>)

add_library(phaseretf SHARED
    $<TARGET_OBJECTS:commonsingle> $<TARGET_OBJECTS:single>)


if (MSVC)
    set_target_properties(phaseret PROPERTIES IMPORT_SUFFIX ".dll.lib")
    set_target_properties(phaseretf PROPERTIES IMPORT_SUFFIX ".dll.lib")
    set_target_properties(phaseretd PROPERTIES IMPORT_SUFFIX ".dll.lib")
endif(MSVC)



target_link_libraries(phaseret  ${LTFAT_LIB} ${LIBS})
target_link_libraries(phaseretf ${LTFATF_LIB} ${LIBS})
target_link_libraries(phaseretd ${LTFATD_LIB} ${LIBS})


set(CMAKE_BUILD_TYPE Release)
