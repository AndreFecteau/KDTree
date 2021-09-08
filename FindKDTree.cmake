# Distributed under the MPL 2 License.  See accompanying
# file LICENSE.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
FindKDTree
-------

Finds the KDTree library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``KDTree::KDTree``
  The KDTree library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``KDTree_FOUND``
  True if the system has the KDTree library.
``KDTree_VERSION``
  The version of the KDTree library which was found.
``KDTree_INCLUDE_DIRS``
  Include directories needed to use KDTree.
``KDTree_LIBRARIES``
  Libraries needed to link to KDTree.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``KDTree_INCLUDE_DIR``
  The directory containing ``KDTree.h``.
``KDTree_LIBRARY``
  The path to the KDTree library.

#]=======================================================================]

include(FindPackageHandleStandardArgs)
set(KDTree_VERSION 1.0.0)

find_path(KDTree_INCLUDE_DIR
        NAMES KDTree.h
        PATHS KDTree
        PATH_SUFFIXES include
        NO_DEFAULT_PATH
        )

find_package_handle_standard_args(KDTree DEFAULT_MSG
        FOUND_VAR KDTree_FOUND
        REQUIRED_VARS KDTree_INCLUDE_DIR
        VERSION_VAR KDTree_VERSION
        )

if (KDTree_FOUND AND NOT TARGET KDTree::KDTree)
    add_library(KDTree::KDTree INTERFACE IMPORTED)
    set_target_properties(KDTree::KDTree PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${KDTree_INCLUDE_DIR}"
            )
endif ()

mark_as_advanced(
        KDTree_INCLUDE_DIR
)