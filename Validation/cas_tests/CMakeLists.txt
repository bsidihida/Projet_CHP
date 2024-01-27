# -------------------------------------------------------------------------
#
#  NEOS
#
#  -------------------------------------------------------------------------
#  License
#  This file is part of Neos.
#
#  Neos is free software: you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License v3 (LGPL)
#  as published by the Free Software Foundation.
#
#  Neos is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with Neos. If not, see <http://www.gnu.org/licenses/>.
#
#---------------------------------------------------------------------------

set(NEOS_TEST_SRCS "")

if (BUILD_INTERPOLATOR)
  list(APPEND NEOS_TEST_SRCS
    interpolator.cpp
    )
endif()

if (BUILD_LAPLACIAN)
  list(APPEND NEOS_TEST_SRCS
    laplacian.cpp
    )
endif()

if (BUILD_CORE)
  list(APPEND NEOS_TEST_SRCS
    grid.cpp
    )
endif()

if (BUILD_GEOMETRY AND BUILD_CORE)
  list(APPEND NEOS_TEST_SRCS
    neos.cpp
    )
endif()

if (BUILD_GEOMETRY AND BUILD_CORE AND BUILD_INTERPOLATOR)
  list(APPEND NEOS_TEST_SRCS
    transport.cpp
    )
endif()

add_executable(
  neos_tests
  main.cpp
  ${NEOS_TEST_SRCS}
  )

add_dependencies(neos_tests googletest)

target_link_libraries(
  neos_tests
  googletest
  neos
  pthread
  )
target_include_directories(neos_tests PRIVATE ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT})
target_include_directories(neos_tests PRIVATE ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}/include)

add_test(interpolator ${PROJECT_BINARY_DIR}/tests/neos_tests)
