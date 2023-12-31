###############################################################################
# find swig
#find_package(SWIG REQUIRED)
#include(${SWIG_USE_FILE})

###############################################################################
# Build swig

option(ITK_USE_SYSTEM_SWIG "Use system swig. If OFF, swig is built as an external project." OFF)

set(swig_version 2.0.4)
set(swig_md5     4319c503ee3a13d2a53be9d828c3adc0)
set(swigwin_md5  4ab8064b1a8894c8577ef9d0fb2523c8)

if(WIN32)
  set(swig_ep "${CMAKE_CURRENT_BINARY_DIR}/swigwin-${swig_version}/swig.exe")
else()
  set(swig_ep "${CMAKE_CURRENT_BINARY_DIR}/swig/bin/swig")
endif()

if(ITK_USE_SYSTEM_SWIG)
  # the path set for the EP build prevents find_package to do its job
  if("${SWIG_EXECUTABLE}" STREQUAL "${swig_ep}")
    unset(SWIG_DIR CACHE)
    unset(SWIG_EXECUTABLE CACHE)
    unset(SWIG_VERSION CACHE)
  endif()
  find_package(SWIG REQUIRED)
  # set(SWIG_EXECUTABLE ${swig} CACHE FILEPATH "swig executable." FORCE)
else()
  include(ExternalProject)
  if(WIN32)
    ExternalProject_Add(swig
      URL http://prdownloads.sourceforge.net/swig/swigwin-${swig_version}.zip
      URL_MD5 ${swigwin_md5}
      SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/swigwin-${swig_version}
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ""
      STEP_TARGETS install
      )
    set(SWIG_DIR ${CMAKE_CURRENT_BINARY_DIR}/swigwin-${swig_version})
    set(SWIG_VERSION "${swig_version}")
  else()
    # build swig as an external project
    ExternalProject_add(PCRE
    URL http://downloads.sourceforge.net/project/pcre/pcre/8.12/pcre-8.12.tar.gz
    URL_MD5 fa69e4c5d8971544acd71d1f10d59193
    PREFIX PCRE
    CONFIGURE_COMMAND ../PCRE/configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/PCRE --enable-shared=no
    )
    ExternalProject_add(swig
    URL http://prdownloads.sourceforge.net/swig/swig-${swig_version}.tar.gz
    URL_MD5 ${swig_md5}
    PREFIX swig
    CONFIGURE_COMMAND ../swig/configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/swig
    --with-pcre-prefix=${CMAKE_CURRENT_BINARY_DIR}/PCRE
    DEPENDS PCRE
    STEP_TARGETS install
    )
    set(SWIG_DIR ${CMAKE_CURRENT_BINARY_DIR}/swig/share/swig/${swig_version} CACHE FILEPATH "swig directory." FORCE)
  endif()
  set(SWIG_EXECUTABLE ${swig_ep} CACHE FILEPATH "swig executable." FORCE)
endif()


###############################################################################
# find python
find_package(PythonInterp REQUIRED)

set(CABLE_INDEX ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/idx.py CACHE INTERNAL "cableidx path")

###############################################################################
# install the files requires for swiginterface
if(NOT EXTERNAL_WRAP_ITK_PROJECT)
  WRAP_ITK_INSTALL(/Configuration/Generators/SwigInterface CMakeLists.txt)
  WRAP_ITK_INSTALL(/Configuration/Generators/SwigInterface Master.mdx.in)
  WRAP_ITK_INSTALL(/Configuration/Generators/SwigInterface empty.in)
  WRAP_ITK_INSTALL(/Configuration/Generators/SwigInterface module.i.in)
  WRAP_ITK_INSTALL(/Configuration/Generators/SwigInterface module.includes.in)

  install(DIRECTORY pygccxml-1.0.0
    DESTINATION "${WRAP_ITK_INSTALL_PREFIX}/Configuration/Generators/SwigInterface"
    PATTERN ".svn" EXCLUDE
    PATTERN "CVS" EXCLUDE
    PATTERN "docs" EXCLUDE
    PATTERN "unittests" EXCLUDE
  )

  install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/idx.py"
                "${CMAKE_CURRENT_SOURCE_DIR}/igenerator.py"
    DESTINATION "${WRAP_ITK_INSTALL_PREFIX}/Configuration/Generators/SwigInterface"
  )

endif(NOT EXTERNAL_WRAP_ITK_PROJECT)


###############################################################################
# store the current dir, so it can be reused later
set(ITK_WRAP_SWIGINTERFACE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE INTERNAL "swig interface source dir")
set(ITK_WRAP_SWIGINTERFACE_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}" CACHE INTERNAL "swig interface binary dir")

set(WRAPPER_MASTER_INDEX_OUTPUT_DIR "${PROJECT_BINARY_DIR}/Typedefs" CACHE INTERNAL "typedefs dir")
set(IGENERATOR  "${CMAKE_CURRENT_SOURCE_DIR}/igenerator.py" CACHE INTERNAL "igenerator.py path" FORCE)

macro(itk_wrap_module_swig_interface library_name)
  # store the content of the mdx file
  set(SWIG_INTERFACE_MDX_CONTENT )
  # store the content of the .i file for the module - a set of import of all the .i files generated for the module
  set(SWIG_INTERFACE_MODULE_CONTENT )
  # the list of .idx files generated for the module
  set(SWIG_INTERFACE_IDX_FILES )
  # build a list of modules to create the ignerator custom command in
  # itk_end_wrap_module_swig_interface
  set(SWIG_INTERFACE_MODULES )
endmacro(itk_wrap_module_swig_interface)

macro(itk_end_wrap_module_swig_interface)

  # Loop over the extra swig input files and copy them to the Typedefs directory
  foreach(source ${WRAPPER_LIBRARY_SWIG_INPUTS})
    get_filename_component(basename ${source} NAME)
    set(dest "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${basename}")
    exec_program(${CMAKE_COMMAND} ARGS -E copy_if_different "\"${source}\"" "\"${dest}\"")
    WRAP_ITK_INSTALL("/Configuration/Typedefs" "${dest}")
#    set(SWIG_INTERFACE_MODULE_CONTENT "${SWIG_INTERFACE_MODULE_CONTENT}%import ${basename}\n")
  endforeach(source)

  # the list of .i files generated for the module
  set(SWIG_INTERFACE_FILES )

  # prepare dependencies
  set(DEPS )
  foreach(dep ${WRAPPER_LIBRARY_DEPENDS})
    set(DEPS ${DEPS} ${${dep}IdxFiles} ${${dep}SwigFiles})
    set(SWIG_INTERFACE_MDX_CONTENT "${dep}.mdx\n${SWIG_INTERFACE_MDX_CONTENT}")
  endforeach(dep)

  # add some libs required by this module
  set(swig_libs )
  foreach(swig_lib ${WRAPPER_SWIG_LIBRARY_FILES})
    get_filename_component(basename ${swig_lib} NAME)
    set(swig_libs ${swig_libs} --swig-include ${basename})
    set(dest "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${basename}")
    exec_program(${CMAKE_COMMAND} ARGS -E copy_if_different "\"${swig_lib}\"" "\"${dest}\"")
    WRAP_ITK_INSTALL("/Configuration/Typedefs" "${dest}")
  endforeach(swig_lib WRAPPER_SWIG_LIBRARY_FILES)

  # accumulate the idx files already generated before this module to generate usable depenencies
  set(idx_files )
  foreach(module ${SWIG_INTERFACE_MODULES})
    # create the swig interface
    set(interface_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${module}.i")
    set(xml_file "${WRAPPER_LIBRARY_OUTPUT_DIR}/${module}.xml")
    set(idx_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${module}.idx")
    set(idx_files ${idx_files} ${idx_file})
#    set(includes_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${WRAPPER_LIBRARY_NAME}.includes")
    set(module_interface_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${WRAPPER_LIBRARY_NAME}.i")
    set(typedef_in_file "${WRAPPER_LIBRARY_OUTPUT_DIR}/${module}SwigInterface.h.in")
    set(typedef_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${module}SwigInterface.h")
    # prepare the options
    set(opts )
    foreach(dep ${WRAPPER_LIBRARY_DEPENDS})
      set(opts ${opts} --mdx "${WRAP_ITK_TYPEDEFS_DIRECTORY}/${dep}.mdx")
#      set(opts ${opts} --include "${dep}.includes")
#      set(opts ${opts} --import "${dep}.i")
    endforeach(dep)
    # import the interface files previously defined instead of importing all the files defined
    # in the module to avoid many warnings when running swig
#    foreach(i ${SWIG_INTERFACE_FILES})
#      get_filename_component(bi "${i}" NAME)
#      set(opts ${opts} --import "${bi}")
#    endforeach(i)

  # message("${opts}")

    # configure the test driver, to set the python path to the pygccxml dir
  #   find_program(ITK_TEST_DRIVER itkTestDriver)
  #   set(PYTHON_PYGCCXML_DRIVER "${ITK_TEST_DRIVER}"
  #     --add-before-env PYTHONPATH "${WRAP_ITK_CMAKE_DIR}/pygccxml-1.0.0/"
  #     "${PYTHON_EXECUTABLE}"
  #   )

    if(ITK_WRAP_EXPLICIT)
      set(opts ${opts} --include "${WRAPPER_LIBRARY_NAME}Explicit.h")
    endif(ITK_WRAP_EXPLICIT)

    add_custom_command(
      OUTPUT ${interface_file} ${typedef_file}
      COMMAND ${PYTHON_EXECUTABLE} ${IGENERATOR}
        ${opts}
        --swig-include itk.i
        ${swig_libs}
        --mdx ${mdx_file}
#        --include ${WRAPPER_LIBRARY_NAME}.includes
#        --import ${module_interface_file}
        --swig-include ${module}_ext.i
        -w1 -w3 -w51 -w52 -w53 -w54 #--warning-error
        -A protected -A private
        --typedef-input ${typedef_in_file}
        --typedef-output ${typedef_file}
        --include ${module}SwigInterface.h
        # --verbose
        ${xml_file}
        ${interface_file}
      DEPENDS ${DEPS} ${idx_files} ${IGENERATOR} # ${SWIG_INTERFACE_IDX_FILES} ${SWIG_INTERFACE_FILES} ${typedef_in_file} # ${includes_file}
    )
  #   add_custom_target(${module}Swig DEPENDS ${interface_file})
  #   add_dependencies(${module}Swig ${WRAPPER_LIBRARY_NAME}Idx)
  #   add_dependencies(${WRAPPER_LIBRARY_NAME}Swig ${module}Swig)

    set(SWIG_INTERFACE_FILES ${SWIG_INTERFACE_FILES} ${interface_file})

    WRAP_ITK_INSTALL("/Configuration/Typedefs" "${interface_file}")

  endforeach(module)


  # the mdx file
  set(mdx_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${WRAPPER_LIBRARY_NAME}.mdx")
  set(CONFIG_INDEX_FILE_CONTENT "${SWIG_INTERFACE_MDX_CONTENT}")
  configure_file("${ITK_WRAP_SWIGINTERFACE_SOURCE_DIR}/Master.mdx.in" "${mdx_file}"
     @ONLY IMMEDIATE )
  WRAP_ITK_INSTALL("/Configuration/Typedefs" "${mdx_file}")

  set(module_interface_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${WRAPPER_LIBRARY_NAME}.i")
  set(module_interface_ext_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${WRAPPER_LIBRARY_NAME}_ext.i")
  set(deps_imports )
  set(deps_includes )
  foreach(dep ${WRAPPER_LIBRARY_DEPENDS})
    set(deps_imports "${deps_imports}%import ${dep}.i\n")
#    set(deps_includes "${deps_includes}#include \"${dep}.includes\"\n")
  endforeach(dep)

#  set(SWIG_INTERFACE_INCLUDES "${deps_includes}#include \"${WRAPPER_LIBRARY_NAME}.includes\"")
  set(CONFIG_MODULE_INTERFACE_CONTENT ) #"${deps_imports}${SWIG_INTERFACE_MODULE_CONTENT}")
  configure_file("${ITK_WRAP_SWIGINTERFACE_SOURCE_DIR}/module.i.in" "${module_interface_file}"
    @ONLY IMMEDIATE )
  WRAP_ITK_INSTALL("/Configuration/Typedefs/" "${module_interface_file}")

#  set(WRAP_ITK_FILE_CONTENT )
#  configure_file("${WRAP_ITK_CONFIG_DIR}/empty.in" "${module_interface_ext_file}"
#    @ONLY IMMEDIATE )
#  install(FILES "${module_interface_ext_file}"
#    DESTINATION "${WRAP_ITK_INSTALL_PREFIX}/Configuration/Typedefs/"
#  )

  add_custom_target(${WRAPPER_LIBRARY_NAME}Idx DEPENDS ${SWIG_INTERFACE_IDX_FILES})
  add_dependencies(${WRAPPER_LIBRARY_NAME}Idx ${WRAPPER_LIBRARY_NAME}GccXML)
  set(${WRAPPER_LIBRARY_NAME}IdxFiles ${SWIG_INTERFACE_IDX_FILES} CACHE INTERNAL "Internal ${WRAPPER_LIBRARY_NAME}Idx file list.")

  add_custom_target(${WRAPPER_LIBRARY_NAME}Swig DEPENDS ${SWIG_INTERFACE_FILES})
  set(${WRAPPER_LIBRARY_NAME}SwigFiles ${SWIG_INTERFACE_FILES} CACHE INTERNAL "Internal ${WRAPPER_LIBRARY_NAME}Swig file list.")
  add_dependencies(${WRAPPER_LIBRARY_NAME}Swig ${WRAPPER_LIBRARY_NAME}Idx)
  if(NOT EXTERNAL_WRAP_ITK_PROJECT)
    # don't depends on the targets from wrapitk in external projects
    foreach(dep ${WRAPPER_LIBRARY_DEPENDS})
      add_dependencies(${WRAPPER_LIBRARY_NAME}Swig ${dep}Swig ${dep}Idx)
    endforeach(dep)
  endif(NOT EXTERNAL_WRAP_ITK_PROJECT)

endmacro(itk_end_wrap_module_swig_interface)


macro(itk_wrap_include_swig_interface include_file)
  list(APPEND SWIG_INTERFACE_INCLUDES ${include_file})
endmacro(itk_wrap_include_swig_interface)


macro(itk_wrap_submodule_swig_interface module)
  # store the content of the SwigInterface.h files - a set of #includes for that module
  set(SWIG_INTERFACE_INCLUDES )
  # typedefs for swig
  set(SWIG_INTERFACE_TYPEDEFS )
endmacro(itk_wrap_submodule_swig_interface)

macro(itk_end_wrap_submodule_swig_interface module)
  # variables used:
  # WRAPPER_LIBRARY_NAME
  # WRAPPER_LIBRARY_OUTPUT_DIR
  # WRAPPER_LIBRARY_CABLESWIG_INPUTS
  # WRAPPER_LIBRARY_DEPENDS
  # WRAPPER_MASTER_INDEX_OUTPUT_DIR
  # MODULE_INCLUDES


  set(SWIG_INTERFACE_INCLUDES_CONTENT )
#  foreach(dep ${WRAPPER_LIBRARY_DEPENDS})
#    set(SWIG_INTERFACE_INCLUDES_CONTENT "${SWIG_INTERFACE_INCLUDES_CONTENT}#include \"${dep}.includes\"\n")
#  endforeach(dep)
  if(SWIG_INTERFACE_INCLUDES)
    list(REMOVE_DUPLICATES SWIG_INTERFACE_INCLUDES)
    foreach(include_file ${SWIG_INTERFACE_INCLUDES})
      if("${include_file}" MATCHES "<.*>")
        set(SWIG_INTERFACE_INCLUDES_CONTENT "${SWIG_INTERFACE_INCLUDES_CONTENT}#include ${include_file}\n")
      else()
        set(SWIG_INTERFACE_INCLUDES_CONTENT "${SWIG_INTERFACE_INCLUDES_CONTENT}#include \"${include_file}\"\n")
      endif()
    endforeach()
  endif()
  # create the file which store all the includes
  set(includes_file "${WRAPPER_LIBRARY_OUTPUT_DIR}/${module}SwigInterface.h.in")
  configure_file("${ITK_WRAP_SWIGINTERFACE_SOURCE_DIR}/module.includes.in"
     ${includes_file}
     @ONLY IMMEDIATE )
  WRAP_ITK_INSTALL("/Configuration/Typedefs/" "${includes_file}")

  # the master idx file
  set(mdx_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${WRAPPER_LIBRARY_NAME}.mdx")

  # generate the idx file for all the groups of the module
  set(idx_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${module}.idx")
  add_custom_command(
    OUTPUT ${idx_file}
    COMMAND ${CABLE_INDEX}
          ${xml_file} ${idx_file}
    DEPENDS ${CABLE_INDEX} ${xml_file}
  )
  WRAP_ITK_INSTALL("/Configuration/Typedefs/" "${idx_file}")
#   add_custom_target(${module}Idx DEPENDS ${idx_file})

  # store the path of the idx file to store it in the mdx file
  set(SWIG_INTERFACE_MDX_CONTENT "${SWIG_INTERFACE_MDX_CONTENT}${module}.idx\n")
  set(SWIG_INTERFACE_IDX_FILES ${SWIG_INTERFACE_IDX_FILES} ${idx_file})

  set(SWIG_INTERFACE_MODULE_CONTENT "${SWIG_INTERFACE_MODULE_CONTENT}%import ${module}.i\n")

  set(SWIG_INTERFACE_MODULES ${SWIG_INTERFACE_MODULES} ${module})

#  set(interface_ext_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${module}_ext.i")
#  set(WRAP_ITK_FILE_CONTENT )
#  configure_file("${WRAP_ITK_CONFIG_DIR}/empty.in" "${interface_ext_file}"
#    @ONLY IMMEDIATE )
#  install(FILES "${interface_ext_file}"
#    DESTINATION "${WRAP_ITK_INSTALL_PREFIX}/Configuration/Typedefs/"
#  )

endmacro(itk_end_wrap_submodule_swig_interface)


macro(itk_wrap_one_type_swig_interface wrap_method wrap_class swig_name)
  # Add one  typedef to WRAPPER_TYPEDEFS
  # 'wrap_method' is the one of the valid WRAPPER_WRAP_METHODS from itk_wrap_class,
  # 'wrap_class' is the fully-qualified C++ name of the class
  # 'swig_name' is what the swigged class should be called
  # The optional last argument is the template parameters that should go between
  # the < > brackets in the C++ template definition.
  # Only pass 3 parameters to wrap a non-templated class
  #
  # Global vars used: none
  # Global vars modified: WRAPPER_TYPEDEFS

  # get the base C++ class name (no namespaces) from wrap_class:
  string(REGEX REPLACE "(.*::)" "" base_name "${wrap_class}")

  set(wrap_pointer 0)
  set(template_parameters "${ARGV3}")
  if(template_parameters)
    set(full_class_name "${wrap_class}< ${template_parameters} >")
  else(template_parameters)
    set(full_class_name "${wrap_class}")
  endif(template_parameters)

  # itk_wrap_one_type_all_generators("${wrap_method}" "${wrap_class}" "${swig_name}" "${ARGV3}")

  # Add a typedef for the class. We have this funny looking full_name::base_name
  # thing (it expands to, for example "typedef itk::Foo<baz, 2>::Foo") used
  # for gccxml typedefs

  if("${wrap_method}" MATCHES "2_SUPERCLASSES")
    itk_wrap_simple_type_swig_interface("${full_class_name}::Superclass::Superclass" "${swig_name}_Superclass_Superclass")
    itk_wrap_simple_type_swig_interface("${full_class_name}::Superclass::Superclass::Pointer" "${swig_name}_Superclass_Superclass_Pointer")
  endif("${wrap_method}" MATCHES "2_SUPERCLASSES")

  if("${wrap_method}" MATCHES "SUPERCLASS")
    itk_wrap_simple_type_swig_interface("${full_class_name}::Superclass" "${swig_name}_Superclass")
    itk_wrap_simple_type_swig_interface("${full_class_name}::Superclass::Pointer" "${swig_name}_Superclass_Pointer")
  endif("${wrap_method}" MATCHES "SUPERCLASS")

  # the same output with or without FORCE_INSTANTIATE
  itk_wrap_simple_type_swig_interface("${full_class_name}" "${swig_name}")

  if("${wrap_method}" MATCHES "POINTER")
    if("${wrap_method}" STREQUAL "AUTOPOINTER")
      # add a pointer typedef if we are so asked
      itk_wrap_simple_type_swig_interface("${full_class_name}::SelfAutoPointer" "${swig_name}_AutoPointer")
    else("${wrap_method}" STREQUAL "AUTOPOINTER")
      # add a pointer typedef if we are so asked
      itk_wrap_simple_type_swig_interface("${full_class_name}::Pointer" "${swig_name}_Pointer")
    endif("${wrap_method}" STREQUAL "AUTOPOINTER")
  endif("${wrap_method}" MATCHES "POINTER")

endmacro(itk_wrap_one_type_swig_interface)


macro(itk_wrap_simple_type_swig_interface wrap_class swig_name)
  set(SWIG_INTERFACE_TYPEDEFS "${SWIG_INTERFACE_TYPEDEFS}typedef ${wrap_class} ${swig_name};\n")
endmacro(itk_wrap_simple_type_swig_interface)

include_directories("${WRAPPER_MASTER_INDEX_OUTPUT_DIR}")
