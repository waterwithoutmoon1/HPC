executable(
	'milestone05',
	'main.cpp',
	include_directories : [lib_incdirs],
	link_with : [lib],
	dependencies : [eigen, mpi]

)

# Define the source and destination paths
source_file = 'lj54.xyz'
destination_file = 'lj54.xyz'  # This is the name it will have in the build directory

# Use configure_file to copy the file
configure_file(
	input : source_file,
	output : destination_file,
	copy: true
)