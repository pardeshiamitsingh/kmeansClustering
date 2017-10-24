// Amit pardeshi
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>



// calculate distance between two data points
float calculate_distance(float *point1, float *point2, int no_of_data_points) {
	float distance = 0.0;
	for (int i = 0; i<no_of_data_points; i++) {
		float tempDist = point1[i] - point2[i];
		tempDist = tempDist * tempDist;
		distance += tempDist;
	}
	return distance;
}


// data generation
float* generate_data(int no_of_elements) {
	float *rand_data = (float *)malloc(sizeof(float) * no_of_elements);
	for (int i = 0; i < no_of_elements; i++) {
		rand_data[i] = (rand()/ (float)(RAND_MAX)) * 100;
		//rand_data[i] = ((float)rand() / (float)(RAND_MAX)) * i;
	}
	return rand_data;
}

void print_centroids(float * centroids, int cluster_num, int dim) {

	printf("=======PRINT CENTROIDS==========\n");
	for (int i = 0; i < cluster_num * dim; i++) {
		printf("%f \t", centroids[i]);
	}
	printf("\n");

}

int main(int argc, char** argv) {

	int total_records = 12000;

	int dim = 8; // dimensions.


	MPI_Init(NULL, NULL);
	int my_rank, my_world;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &my_world);
	int data_points_per_process = total_records / my_world;
	int num_of_cluster = my_world - 1;// number 0f clusters.

	float* receive_buffer_non_root = malloc(dim * data_points_per_process * sizeof(float));

	float* sums = malloc(dim * num_of_cluster * sizeof(float));

	int* counts = malloc(num_of_cluster * sizeof(int));

	// array holding centroid info
	float* current_centroids = malloc(dim * num_of_cluster * sizeof(float));
	//  cluster on each process
	int* final_tag = malloc(data_points_per_process * sizeof(int));


	float* dataset = NULL;

	float* sum_total = NULL;
	int* total_counts = NULL;
	int* tags_all = NULL;
	if (my_rank == 0) {
		dataset = generate_data(dim * data_points_per_process * my_world);
		// intial centroids
		for (int i = 0; i < dim * num_of_cluster; i++) {
			current_centroids[i] = dataset[i];
		}
		print_centroids(current_centroids, num_of_cluster, dim);
		sum_total = malloc(dim * num_of_cluster * sizeof(float));
		total_counts = malloc(num_of_cluster * sizeof(int));
		tags_all = malloc(my_world * data_points_per_process * sizeof(int));
	}
	// send to all process.
	MPI_Scatter(dataset, dim *data_points_per_process, MPI_FLOAT, receive_buffer_non_root,
		dim*data_points_per_process, MPI_FLOAT, 0, MPI_COMM_WORLD);


	float flag = 1.0;

	while (flag > 0.00001) { // Check for threshhold

		// Broadcast the current cluster centroids to all processes.
		MPI_Bcast(current_centroids, dim * num_of_cluster, MPI_FLOAT, 0, MPI_COMM_WORLD);


		for (int i = 0; i < num_of_cluster; i++) {
			counts[i] = 0;
		}

		for (int i = 0; i < dim * num_of_cluster; i++) {
			sums[i] = 0.0;
		}

		float* current_data_point = receive_buffer_non_root;
		//for each datpoint calcluate the nearest centroid and update the sum of the centroid
		for (int i = 0; i < data_points_per_process; i++) {

			//recalculate datapoint to new cluster
			int new_cluster_index = 0;
			float dist = calculate_distance(current_data_point, current_centroids, dim);
			float* nextCentroid = current_centroids + dim;
			for (int c = 1; c < num_of_cluster; c++, nextCentroid += dim) {
				float new_dist = calculate_distance(current_data_point, nextCentroid, dim);
				if (dist > new_dist) {
					new_cluster_index = c;
					dist = new_dist;
				}
			}

			counts[new_cluster_index]++;
			// update the sum vector for centroid
			float *sum_of_cluster = &sums[new_cluster_index*dim];
			for (int i = 0; i<dim; i++) {
				sum_of_cluster[i] += current_data_point[i];
			}

			current_data_point += dim;
		}

		MPI_Reduce(sums, sum_total, dim * num_of_cluster, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(counts, total_counts, num_of_cluster, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

		if (my_rank == 0) {
			// calculate the new centroids
			for (int i = 0; i<num_of_cluster; i++) {
				for (int j = 0; j<dim; j++) {
					int tmp = dim*i + j;
					sum_total[tmp] = sum_total[tmp] / total_counts[i];
				}
			}
			// check for centroid change flag
			flag = calculate_distance(sum_total, current_centroids, dim * num_of_cluster);
			printf("====THE CHANGE FLAG==== %f\n", flag);
			
			for (int i = 0; i<dim * num_of_cluster; i++) {
				current_centroids[i] = sum_total[i];
			}
			print_centroids(current_centroids, num_of_cluster, dim);
		}
		// broadcast flag to each process
		MPI_Bcast(&flag, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	}

	// centroid are fixed .
	float* data_point = receive_buffer_non_root;
	for (int i = 0; i < data_points_per_process; i++, data_point += dim) {
		int new_cluster_index = 0;
		float new_dist = calculate_distance(data_point, current_centroids, dim);
		float* nextCentroid = current_centroids + dim;
		for (int c = 1; c < num_of_cluster; c++, nextCentroid += dim) {
			float dist = calculate_distance(data_point, nextCentroid, dim);
			if (dist < new_dist) {
				new_cluster_index = c;
				new_dist = dist;
			}
		}
		final_tag[i] = new_cluster_index;
	}

	MPI_Gather(final_tag, data_points_per_process, MPI_INT,
		tags_all, data_points_per_process, MPI_INT, 0, MPI_COMM_WORLD);

	if ((my_rank == 0)) {
		float* data = dataset;
		for (int i = 0; i < my_world * data_points_per_process; i++) {
			for (int j = 0; j < dim; j++)
			{
				printf("%f \t", data[j]);
			}
			printf("cluster%d\n", tags_all[i] + 1);
			data += dim;
		}
	}

	MPI_Finalize();

}