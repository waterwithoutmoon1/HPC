executable(
	'milestone09',
	'main.cpp',
	include_directories : [lib_incdirs],
	link_with : [lib],
	dependencies : [eigen, mpi]

)
# Define the source and destination paths
source_file = 'whisker_small.xyz'
destination_file = 'whisker_small.xyz'  # This is the name it will have in the build directory

# Use configure_file to copy the file
configure_file(
	input : source_file,
	output : destination_file,
	copy: true
)

# Define the source and destination paths
source_file_1 = 'whisker_large.xyz'
destination_file_1 = 'whisker_large.xyz'  # This is the name it will have in the build directory

# Use configure_file to copy the file
configure_file(
	input : source_file_1,
	output : destination_file_1,
	copy: true
)