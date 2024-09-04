lib_sources = [  # All source files (excluding headers)
    'hello.cpp',
    'verlet.cpp',
    'verlet_2.cpp',
    'lj_direction_summation.cpp',
    'xyz.cpp',
    'neighbors.cpp',
    'lj_neighborlist.cpp',
    'ducastelle.cpp',
    'domain.cpp',
    'ih.cpp'
]

lib = static_library(
    'lib',  # Unique name
    lib_sources,
    dependencies : [eigen, mpi]
)

lib_incdirs = include_directories('.')
