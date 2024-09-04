executable(
	'milestone07',
	'main.cpp',
	include_directories : [lib_incdirs],
	link_with : [lib],
	dependencies : [eigen, mpi]

)

# Define the source and destination paths
source_file = 'cluster_923.xyz'
destination_file = 'cluster_923.xyz'  # This is the name it will have in the build directory

# Use configure_file to copy the file
configure_file(
	input : source_file,
	output : destination_file,
	copy: true
)

# Define the source and destination paths
source_file_1 = 'cluster_3871.xyz'
destination_file_1 = 'cluster_3871.xyz'  # This is the name it will have in the build directory

# Use configure_file to copy the file
configure_file(
	input : source_file_1,
	output : destination_file_1,
	copy: true
)